# custom boundary conditions for 2D FermiHarmonics equations

# ======================================================================================================================
# Cache for BC Projectors
# ======================================================================================================================

# In our harmonic basis, boundary conditions require projecting the solution into incoming
# and outgoing characteristic modes (as defined by positive and negative eigenvalues of the
# flux Jacobians Ax and Ay). This is expensive, so it is important to precompute them.

"""
    BCProjectorCache

Precompute and store the projection matrices for all boundary faces/nodes. These are meant
to be thread-local to avoid repeated allocations
"""
mutable struct BCProjectorCache
    state_buffers::Vector{Vector{Float64}}
    target_buffers::Vector{Vector{Float64}}
    out_buffers::Vector{Vector{Float64}}
    projectors::Dict{Tuple{Int, Int}, SparseMatrixCSC{Float64, Int}}
    initialized::Bool
    nvars::Int
end
BCProjectorCache() = BCProjectorCache(
    [Float64[] for _ in 1:Threads.nthreads()],
    [Float64[] for _ in 1:Threads.nthreads()],
    [Float64[] for _ in 1:Threads.nthreads()],
    Dict{Tuple{Int, Int}, SparseMatrixCSC{Float64, Int}}(),
    false,
    0)


# ======================================================================================================================
# Custom Boundary Condition Types
# ======================================================================================================================

"""
    MaxwellWallBC(p_scatter; tol=0.0)

Wall boundary condition with diffuse/specular mixing. It blends between them with parameter p_scatter. The
physics of the diffuse case is that particles scatter from the wall at a random angle with an isotropic distribution,
while the specular case is that particles reflect with the same angle they arrived at. The blended target state is then projected onto the incoming characteristics to get the final BC state.

Parameters:
- `p_scatter`: diffuse fraction (`1.0` fully diffuse, `0.0` fully specular)
- `tol`: eigenvalue tolerance for incoming-mode projector construction

Returns:
- `MaxwellWallBC`.
"""
mutable struct MaxwellWallBC
    p_scatter::Float64 # 1 is diffuse, 0 is specular
    tol::Float64
    cache::BCProjectorCache
end
MaxwellWallBC(p_scatter::Real; tol::Real = 0.0) =
    MaxwellWallBC(Float64(p_scatter), Float64(tol), BCProjectorCache())


"""
    OhmicContactBC(bias; p_ohmic_absorb=1.0, tol=0.0)

This one is meant to model connection to an ohmic contact with fixed potential (or fixed a0, in our case.)
Ohmic contact boundary condition with fixed monopole value (`a0 = bias`).  The physics of the absorbing case is that particles are absorbed at the contact and re-emitted with a distribution corresponding to the imposed bias (fixed value
of only a0, with higher harmonics determined by diffuse scattering). The specular case is that particles reflect with
the same angle they arrived at, with no bias. The blended target state is then projected onto the incoming modes to
get the final BC state.  Note that this means in the case p_ohmic_absorb = 0, we no longer have any particles being
injected.

Parameters:
- `bias`: imposed contact value for the monopole mode
- `p_ohmic_absorb`: absorption fraction (`1.0` fully absorbing, `0.0` fully specular)
- `tol`: eigenvalue tolerance for incoming-mode projector construction

Returns:
- `OhmicContactBC`.
"""
mutable struct OhmicContactBC
    p_ohmic_absorb::Float64
    bias::Float64
    tol::Float64
    cache::BCProjectorCache
end
OhmicContactBC(bias::Real; p_ohmic_absorb::Real = 1.0, tol::Real = 0.0) =
    OhmicContactBC(Float64(p_ohmic_absorb), Float64(bias), Float64(tol), BCProjectorCache())


# helper functions to get a symbol name for the BC type for logging purposes
boundary_condition_name(bc) = nameof(typeof(bc))
boundary_condition_name(::MaxwellWallBC) = :maxwell_wall
boundary_condition_name(::OhmicContactBC) = :ohmic_contact


# ======================================================================================================================
# Maxwell and Ohmic BC implementations
# ======================================================================================================================

# these functions implement the types of boundary condition defined above

"""
    maxwell_wall!(out, state, unit_normal, P_in, p, target) -> out

Apply Maxwell blending boundary condition in-place: blend diffuse and specular targets, then project onto incoming characteristics.
"""
function maxwell_wall!(
    out::AbstractVector{Float64},
    state::AbstractVector{Float64},
    unit_normal::SVector{2, Float64},
    P_in::AbstractSparseMatrix,
    p::Real,
    target::AbstractVector{Float64}
)
    diffuse_target!(target, state, unit_normal)
    specular_target!(out, state, unit_normal)
    p_scatter = Float64(p)
    one_minus = 1.0 - p_scatter
    N = length(state)
    @inbounds for i in 1:N
        target[i] = p_scatter * target[i] + one_minus * out[i]
    end
    apply_projector!(out, state, target, P_in)
    return out
end


"""
    ohmic_contact!(out, state, unit_normal, P_in, p_ohmic_absorb, bias, target) -> out 
Apply ohmic contact boundary condition in-place: blend diffuse (with clamped monopole) and specular targets, then project onto incoming characteristics.
"""
function ohmic_contact!(out::AbstractVector{Float64},
                        state::AbstractVector{Float64},
                        unit_normal::SVector{2, Float64},
                        P_in::AbstractSparseMatrix,
                        p_ohmic_absorb::Real,
                        bias::Real,
                        target::AbstractVector{Float64})
    diffuse_target!(target, state, unit_normal)
    @inbounds target[1] = Float64(bias)
    specular_target!(out, state, unit_normal)
    p_absorb = Float64(p_ohmic_absorb)
    one_minus = 1.0 - p_absorb
    N = length(state)
    @inbounds for i in 1:N
        target[i] = p_absorb * target[i] + one_minus * out[i]
    end
    apply_projector!(out, state, target, P_in)
    return out
end


# ======================================================================================================================
# Diffuse and Specular Target Computation
# ======================================================================================================================

"""
    diffuse_target!(target, state, unit_normal) -> target

Compute diffuse scattering target state in-place.
"""
@inline function diffuse_target!(
    target::AbstractVector{Float64},
    state::AbstractVector{Float64},
    unit_normal::SVector{2, Float64},
)::AbstractVector{Float64}
    c1, s1 = normal_cos_sin(unit_normal)
    N = length(state)
    Mloc = (N - 1) ÷ 2
    @inbounds C = state[1]
    if Mloc >= 1
        cm = c1
        sm = s1
        @inbounds C += 0.5 * pi * (cm * state[cosine_index(1)] + sm * state[sine_index(1)])
        @inbounds for m in 2:Mloc
            cm_next = c1 * cm - s1 * sm
            sm_next = s1 * cm + c1 * sm
            cm = cm_next
            sm = sm_next
            if iseven(m)
                k = m ÷ 2
                sign = isodd(k) ? -1.0 : 1.0
                coeff = -2.0 * sign / (4.0 * k * k - 1.0)
                C += coeff * (cm * state[cosine_index(m)] + sm * state[sine_index(m)])
            end
        end
    end
    @inbounds for i in 1:N
        target[i] = 0.0
    end
    @inbounds target[1] = C
    return target
end

"""
    specular_target!(target, state, unit_normal) -> target

Compute specular reflection target state in-place.
"""
@inline function specular_target!(
    target::AbstractVector{Float64},
    state::AbstractVector{Float64},
    unit_normal::SVector{2, Float64},
)::AbstractVector{Float64}
    c1, s1 = normal_cos_sin(unit_normal)
    N = length(state)
    Mloc = (N - 1) ÷ 2
    @inbounds target[1] = state[1]
    cm = c1
    sm = s1
    @inbounds for m in 1:Mloc
        a = state[cosine_index(m)]
        b = state[sine_index(m)]
        a_rot = cm * a + sm * b
        b_rot = -sm * a + cm * b
        sign = isodd(m) ? -1.0 : 1.0
        a_rot_ref = sign * a_rot
        b_rot_ref = -sign * b_rot
        target[cosine_index(m)] = cm * a_rot_ref - sm * b_rot_ref
        target[sine_index(m)] = sm * a_rot_ref + cm * b_rot_ref
        if m < Mloc
            cm_next = c1 * cm - s1 * sm
            sm_next = s1 * cm + c1 * sm
            cm = cm_next
            sm = sm_next
        end
    end
    return target
end


# ======================================================================================================================
# Kinetic Boundary Conditions via Characteristic Projection
# ======================================================================================================================

# Here we implement the projector computation and application functions needed for the boundary conditions.  Notice that
# a discrete angle method may be more natural for these kinetic BCs.  However, in that case, we would need a large
# number of degrees of freedom even in the hydro, diffusive regimes.  So we work instead in this harmonic basis.

"""
    characteristic_spectral_data(Ax, Ay, unit_normal; tol=0.0)
        -> (; eigvectors, eigenvalues, incoming_mask, outgoing_mask, grazing_mask, Dx, Dxinv, tol_effective)

Compute weighted characteristic data for boundary projectors.

The weighted transform uses `A_s = Dxinv * A * Dx`, which is symmetric for this harmonic
basis, so its eigenvectors define robust incoming/outgoing/grazing subspaces.
"""
function characteristic_spectral_data(
    Ax::AbstractMatrix,
    Ay::AbstractMatrix,
    unit_normal::SVector{2, Float64},
    ;tol::Float64=0.0,
)
    tol >= 0 || throw(ArgumentError("tol must be >= 0"))

    # Weighted similarity transform that symmetrizes the truncated harmonic transport operator.
    # This lets us use an orthogonal eigendecomposition and build stable projectors.
    A = unit_normal[1] .* Ax .+ unit_normal[2] .* Ay
    n = size(A, 1)
    D = ones(Float64, n)
    D[1] = sqrt(2.0)
    Dinv = ones(Float64, n)
    Dinv[1] = 1 / sqrt(2.0)

    Dx = Diagonal(D)
    Dxinv = Diagonal(Dinv)
    A_s = Dxinv * A * Dx
    eig = eigen(Symmetric(A_s))

    # Roundoff-safe grazing band: exact zero-speed modes are physically tangential and should
    # remain unconstrained by incoming BC projection.
    tol_floor = 100.0 * eps(Float64) * max(1.0, opnorm(A_s, Inf))
    tol_effective = max(tol, tol_floor)

    incoming_mask = eig.values .< -tol_effective
    outgoing_mask = eig.values .> tol_effective
    grazing_mask = .!(incoming_mask .| outgoing_mask)

    return (
        eigvectors = eig.vectors,
        eigenvalues = eig.values,
        incoming_mask = incoming_mask,
        outgoing_mask = outgoing_mask,
        grazing_mask = grazing_mask,
        Dx = Dx,
        Dxinv = Dxinv,
        tol_effective = tol_effective,
    )
end

"""
    projector_from_mask(eigvectors, mask, Dx, Dxinv) -> SparseMatrixCSC{Float64, Int}

Construct one characteristic projector from a boolean eigenvector mask in weighted space.
"""
function projector_from_mask(
    eigvectors::AbstractMatrix,
    mask::AbstractVector{Bool},
    Dx::Diagonal,
    Dxinv::Diagonal,
)::SparseMatrixCSC{Float64, Int}
    n = size(eigvectors, 1)
    if !any(mask)
        return spzeros(Float64, n, n)
    end

    Q = eigvectors[:, mask]
    P_s = Q * transpose(Q)
    P_dense = Dx * P_s * Dxinv
    P_sparse = sparse(P_dense)
    droptol!(P_sparse, 1e-12)
    return P_sparse
end

"""
    incoming_projector(Ax, Ay, unit_normal; tol=0.0) -> SparseMatrixCSC{Float64}

Compute incoming projector for flux Jacobian A = Ax*nx + Ay*ny as sparse matrix.
Incoming modes are eigenvectors with eigenvalues strictly `< -tol`.
Grazing modes with `|lambda| <= tol` are left untouched by the boundary projector.
"""
function incoming_projector(
    Ax::AbstractMatrix,
    Ay::AbstractMatrix,
    unit_normal::SVector{2, Float64},
    ;tol::Float64=0.0,
)::SparseMatrixCSC{Float64, Int}
    spec = characteristic_spectral_data(Ax, Ay, unit_normal; tol=tol)
    return projector_from_mask(spec.eigvectors, spec.incoming_mask, spec.Dx, spec.Dxinv)
end

"""
    characteristic_projectors(Ax, Ay, unit_normal; tol=0.0)
        -> (; incoming, outgoing, grazing, eigenvalues, nincoming, noutgoing, ngrazing)

Compute incoming/outgoing/grazing characteristic projectors for
`A(unit_normal) = n_x A_x + n_y A_y`.

Sign convention:
- `lambda < -tol`: incoming-to-domain modes (constrained by BC targets),
- `lambda > +tol`: outgoing-from-domain modes (left as interior state),
- `|lambda| <= tol`: grazing/tangential modes (left untouched to avoid basis-dependent bias).

Numerical note:
- an internal roundoff floor is always applied to `tol` so machine-zero eigenvalues do not
  get misclassified as incoming/outgoing when users pass `tol=0`.
"""
function characteristic_projectors(
    Ax::AbstractMatrix,
    Ay::AbstractMatrix,
    unit_normal::SVector{2, Float64},
    ;tol::Float64=0.0,
)
    spec = characteristic_spectral_data(Ax, Ay, unit_normal; tol=tol)
    P_in = projector_from_mask(spec.eigvectors, spec.incoming_mask, spec.Dx, spec.Dxinv)
    P_out = projector_from_mask(spec.eigvectors, spec.outgoing_mask, spec.Dx, spec.Dxinv)
    P_g = projector_from_mask(spec.eigvectors, spec.grazing_mask, spec.Dx, spec.Dxinv)

    return (
        incoming = P_in,
        outgoing = P_out,
        grazing = P_g,
        eigenvalues = spec.eigenvalues,
        tol_effective = spec.tol_effective,
        nincoming = count(spec.incoming_mask),
        noutgoing = count(spec.outgoing_mask),
        ngrazing = count(spec.grazing_mask),
    )
end

"""
    audit_characteristic_projectors(Ax, Ay, unit_normal; tol=0.0)
        -> NamedTuple

Audit projector algebra and opposite-normal consistency. Useful for checking BC operators:
- idempotence (`P^2 = P`),
- decomposition (`P_in + P_out + P_g = I`),
- orthogonality between subspaces,
- reciprocity under `n -> -n` (`P_in(n) ≈ P_out(-n)` and `P_g(n) ≈ P_g(-n)`).
"""
function audit_characteristic_projectors(
    Ax::AbstractMatrix,
    Ay::AbstractMatrix,
    unit_normal::SVector{2, Float64},
    ;tol::Float64=0.0,
)
    proj = characteristic_projectors(Ax, Ay, unit_normal; tol=tol)
    proj_opposite = characteristic_projectors(Ax, Ay, -unit_normal; tol=tol)

    n = size(Ax, 1)
    Id = Matrix{Float64}(I, n, n)
    P_in = Matrix(proj.incoming)
    P_out = Matrix(proj.outgoing)
    P_g = Matrix(proj.grazing)
    P_in_opposite = Matrix(proj_opposite.incoming)
    P_out_opposite = Matrix(proj_opposite.outgoing)
    P_g_opposite = Matrix(proj_opposite.grazing)

    return (
        tol_effective = proj.tol_effective,
        nincoming = proj.nincoming,
        noutgoing = proj.noutgoing,
        ngrazing = proj.ngrazing,
        idempotence_in = opnorm(P_in * P_in - P_in, Inf),
        idempotence_out = opnorm(P_out * P_out - P_out, Inf),
        idempotence_grazing = opnorm(P_g * P_g - P_g, Inf),
        decomposition_error = opnorm(P_in + P_out + P_g - Id, Inf),
        orthogonality_error = max(
            opnorm(P_in * P_out, Inf),
            opnorm(P_in * P_g, Inf),
            opnorm(P_out * P_g, Inf),
        ),
        opposite_incoming_outgoing_error = opnorm(P_in - P_out_opposite, Inf),
        opposite_outgoing_incoming_error = opnorm(P_out - P_in_opposite, Inf),
        opposite_grazing_error = opnorm(P_g - P_g_opposite, Inf),
    )
end

"""
    apply_projector!(out, state, target, P_in) -> out

Apply incoming projector: out = state + P_in * (target - state). Sparse matrix implementation.
"""
@inline function apply_projector!(
    out::AbstractVector{Float64},
    state::AbstractVector{Float64},
    target::AbstractVector{Float64},
    P_in::AbstractSparseMatrix,
)::AbstractVector{Float64}
    N = length(state)
    @inbounds for i in 1:N
        out[i] = state[i]
    end
    colptr = P_in.colptr
    rowval = P_in.rowval
    nzval = P_in.nzval
    @inbounds for j in 1:N
        dj = target[j] - state[j]
        for idx in colptr[j]:(colptr[j + 1] - 1)
            out[rowval[idx]] += nzval[idx] * dj
        end
    end
    return out
end


# ======================================================================================================================
# Projector Cache Initialization
# ======================================================================================================================

# In this section we implement the caching of projector matrices per face for P4est meshes.
# The main function is `build_projectors` which constructs the sparse projector matrices
# for all faces assigned to a given boundary condition type. These are stored in the
# `BCProjectorCache` struct defined above.

"""
    build_projectors(equations, tol, mesh, dg, cache, boundary_indexing)

Build sparse projector matrices for each `(boundary face, node)` pair in a P4est mesh.
Returns `Dict{Tuple{Int, Int}, SparseMatrixCSC}` keyed by `(boundary_index, node_index)`.
"""
function build_projectors(
    equations::FermiHarmonics2D,
    tol::Float64,
    mesh::Trixi.P4estMesh{2},
    solver,
    cache,
    boundary_indexing::Vector{Int},
)::Dict{Tuple{Int, Int}, SparseMatrixCSC{Float64, Int}}
    n_nodes = Trixi.nnodes(solver)
    contravariant_vectors = cache.elements.contravariant_vectors
    boundaries = cache.boundaries
    projectors = Dict{Tuple{Int, Int}, SparseMatrixCSC{Float64, Int}}()

    # Build projectors at each boundary node so the cached projector uses the same local
    # normal as the runtime flux evaluation. Using one face-representative normal can
    # create directional bias on warped/curved faces.
    for global_idx in boundary_indexing
        element = boundaries.neighbor_ids[global_idx]
        node_indices = boundaries.node_indices[global_idx]
        direction = Trixi.indices2direction(node_indices)
        for node_index in 1:n_nodes
            i_index, j_index = boundary_node_ij(direction, node_index, n_nodes)
            normal_direction = Trixi.get_normal_direction(
                direction, contravariant_vectors, i_index, j_index, element
            )
            unit_n = unit_normal(
                SVector(Float64(normal_direction[1]), Float64(normal_direction[2]))
            )
            key = projector_cache_key(global_idx, node_index)
            projectors[key] = incoming_projector(
                equations.Ax, equations.Ay, unit_n; tol = tol
            )
        end
    end
    return projectors
end

"""
    init_projector_cache!(semi)

Initialize the BC projector cache for all boundary conditions in the semidiscretization.

The projector cache depends only on mesh geometry, solver degree (polydeg), and boundary
locations - NOT on physics parameters. This function checks if the cache is already 
initialized and reuses it across multiple solves, which is critical for parameter sweeps.

# Performance note
Building projectors is expensive (eigendecomposition of flux Jacobians). Reusing the cache
can save significant time when running multiple cases with the same mesh/solver but 
different physics parameters (gamma_mr, gamma_mc, etc.).
"""
function init_projector_cache!(
    semi::Trixi.SemidiscretizationHyperbolic{<:Any, <:FermiHarmonics2D},
)::Trixi.SemidiscretizationHyperbolic{<:Any, <:FermiHarmonics2D}
    boundary_conditions = semi.boundary_conditions
    nvars = Trixi.nvariables(semi.equations)
    
    # Reset caches if number of variables changed (e.g., adaptive harmonics in sweeps)
    for bc in boundary_conditions.boundary_condition_types
        if bc.cache.nvars != 0 && bc.cache.nvars != nvars
            empty!(bc.cache.projectors)
            bc.cache.initialized = false
        end
        bc.cache.nvars = nvars
    end

    # Check if cache already initialized - reuse if possible
    all_initialized = all(bc.cache.initialized for bc in boundary_conditions.boundary_condition_types)
    if all_initialized
        @debug "BC projector cache already initialized, reusing" mesh = typeof(semi.mesh)
        return semi
    end
    
    # Build projectors for any uninitialized boundary conditions
    for (bc, boundary_indexing) in zip(
        boundary_conditions.boundary_condition_types,
        boundary_conditions.boundary_indices,
    )
        if !bc.cache.initialized
            new_projectors = build_projectors(
                semi.equations, bc.tol, semi.mesh, semi.solver, semi.cache, boundary_indexing
            )
            merge!(bc.cache.projectors, new_projectors)
            bc.cache.initialized = true
        end
    end
    @debug "Initialized BC projector cache" mesh = typeof(semi.mesh)
    return semi
end

# ======================================================================================================================
# Helper Functions for Normal Vector Computation
# ======================================================================================================================

"""
    unit_normal(normal::SVector{2, Float64}) -> SVector{2, Float64}
    
Normalize a normal vector to unit length.
"""
@inline function unit_normal(normal::SVector{2, Float64})::SVector{2, Float64}
    nx, ny = normal
    nrm = hypot(nx, ny)
    if nrm < 1e-14
        return SVector(1.0, 0.0)
    else
        return SVector(nx / nrm, ny / nrm)
    end
end

"""
    normal_cos_sin(normal) -> (cosine, sine)

Compute cosine and sine of the angle defined by the normal vector.
"""
@inline function normal_cos_sin(normal::SVector{2, Float64})::Tuple{Float64, Float64}
    nx, ny = normal
    nrm = hypot(nx, ny)
    if nrm < 1e-14
        return 1.0, 0.0
    else
        return nx / nrm, ny / nrm
    end
end

"""
    projector_cache_key(boundary_index, node_index) -> (Int, Int)

Canonical dictionary key for cached boundary projectors.
"""
@inline projector_cache_key(boundary_index::Int, node_index::Int) = (boundary_index, node_index)

"""
    boundary_node_ij(direction, node_index, n_nodes) -> (i_index, j_index)

Map a boundary-local node number to `(i, j)` indices on a tensor-product P4est face.
Direction follows Trixi convention:
- `1,2`: x-normal faces (`i` fixed),
- `3,4`: y-normal faces (`j` fixed).
"""
@inline function boundary_node_ij(
    direction::Int,
    node_index::Int,
    n_nodes::Int,
)::Tuple{Int, Int}
    if direction == 1
        return 1, node_index
    elseif direction == 2
        return n_nodes, node_index
    elseif direction == 3
        return node_index, 1
    elseif direction == 4
        return node_index, n_nodes
    end
    throw(ArgumentError("Unsupported boundary direction $direction"))
end
