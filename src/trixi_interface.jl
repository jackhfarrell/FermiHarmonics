# This script implements required methods by the solver Trixi.jl for the custom FermiHarmonics2D equations type.
# TODO: clean up this file. Do we need all of these methods now that we only support the P4estMesh

# ======================================================================================================================
# Required Equation Methods
# ======================================================================================================================

function Trixi.varnames(::typeof(cons2cons), equations::FermiHarmonics2D)
    nvars = Trixi.nvariables(equations)
    names = Vector{String}(undef, nvars)
    names[1] = "a0"
    idx = 2
    max_harmonic = (nvars - 1) ÷ 2
    for m in 1:max_harmonic
        if idx <= nvars
            names[idx] = "a$(m)"
            idx += 1
        end
        if idx <= nvars
            names[idx] = "b$(m)"
            idx += 1
        end
    end
    return Tuple(names)
end

Trixi.varnames(::typeof(cons2prim), equations::FermiHarmonics2D) =
    Trixi.varnames(cons2cons, equations)

@inline Trixi.cons2prim(u, equations::FermiHarmonics2D) = u

@inline Trixi.cons2cons(u, equations::FermiHarmonics2D) = u

@inline function Trixi.flux(u, orientation::Integer, equations::FermiHarmonics2D{NVARS}) where {NVARS}
    normal = orientation == 1 ? SVector(1.0, 0.0) : SVector(0.0, 1.0)
    out = MVector{NVARS, Float64}(undef)
    harmonics_flux!(out, u, normal)
    return SVector{NVARS, Float64}(out)
end

@inline function Trixi.flux(u, normal_direction::AbstractVector, equations::FermiHarmonics2D{NVARS}) where {NVARS}
    normal = SVector(normal_direction[1], normal_direction[2])
    out = MVector{NVARS, Float64}(undef)
    harmonics_flux!(out, u, normal)
    return SVector{NVARS, Float64}(out)
end

@inline function (dissipation::Trixi.DissipationLocalLaxFriedrichs)(
    u_ll::AbstractVector{Float64}, u_rr::AbstractVector{Float64},
    orientation_or_normal_direction, equations::FermiHarmonics2D{NVARS}
) where {NVARS}
    λ = dissipation.max_abs_speed(u_ll, u_rr, orientation_or_normal_direction, equations)
    coeff = -0.5f0 * λ
    out = MVector{NVARS, Float64}(undef)
    @inbounds for i in 1:NVARS
        out[i] = coeff * (u_rr[i] - u_ll[i])
    end
    return SVector{NVARS, Float64}(out)
end

@inline function Trixi.max_abs_speed_naive(u_ll, u_rr, orientation::Integer,
                                          equations::FermiHarmonics2D)
    return equations.max_speed
end

@inline function Trixi.max_abs_speed_naive(u_ll, u_rr, normal_direction::AbstractVector,
                                          equations::FermiHarmonics2D)
    nrm = hypot(normal_direction[1], normal_direction[2])
    return equations.max_speed * nrm
end

@inline Trixi.have_constant_speed(::FermiHarmonics2D) = Trixi.True()

@inline Trixi.max_abs_speeds(u_or_eq::Union{FermiHarmonics2D, AbstractVector}, 
                            equations::FermiHarmonics2D) = (equations.max_speed, equations.max_speed)
@inline Trixi.max_abs_speeds(equations::FermiHarmonics2D) = (equations.max_speed, equations.max_speed)
# ======================================================================================================================
# Boundary Condition Interface
# ======================================================================================================================

@inline function get_thread_buffer!(buffers::Vector{Vector{Float64}}, nvars::Int)
    tid = Threads.threadid()
    @inbounds buf = buffers[tid]
    if length(buf) != nvars
        buf = Vector{Float64}(undef, nvars)
        buffers[tid] = buf
    end
    return buf
end

@inline function ensure_state_vector(u_inner, buffers::Vector{Vector{Float64}})
    u_inner isa AbstractVector{Float64} && return u_inner
    nvars = length(u_inner)
    buf = get_thread_buffer!(buffers, nvars)
    @inbounds for i in 1:nvars
        buf[i] = u_inner[i]
    end
    return buf
end

@inline function projector_direction_index(orientation::Integer, direction::Integer)
    return orientation == 1 ? (direction == 1 ? 1 : 2) : (direction == 3 ? 3 : 4)
end

@inline function apply_bc!(bc_type::Symbol, out, state, unit_n, P_in, bc, target)
    if bc_type === :maxwell
        maxwell_wall!(out, state, unit_n, P_in, bc.p_scatter, target)
    else
        ohmic_contact!(out, state, unit_n, P_in, bc.p_ohmic_absorb, bc.bias, target)
    end
    return out
end

@inline function bc_callable(bc, bc_type::Symbol, u_inner, normal_direction, x, t,
                             surface_flux_function, equations, 
                             boundary_index::Int=0, node_index::Int=0)
    normal = SVector(Float64(normal_direction[1]), Float64(normal_direction[2]))
    unit_n = unit_normal(normal)
    state = ensure_state_vector(u_inner, bc.cache.state_buffers)
    nvars = length(state)
    target = get_thread_buffer!(bc.cache.target_buffers, nvars)
    out = get_thread_buffer!(bc.cache.out_buffers, nvars)
    if bc.cache.initialized && boundary_index > 0 && haskey(bc.cache.projectors, boundary_index)
        @inbounds P_in = bc.cache.projectors[boundary_index]
    else
        P_in = incoming_projector(equations.Ax, equations.Ay, unit_n; tol = bc.tol)
    end
    apply_bc!(bc_type, out, state, unit_n, P_in, bc, target)
    return surface_flux_function(state, out, normal_direction, equations)
end

@inline bc_type(::MaxwellWallBC) = :maxwell
@inline bc_type(::OhmicContactBC) = :ohmic

@inline function (bc::Union{MaxwellWallBC, OhmicContactBC})(u_inner, normal_direction::AbstractVector, 
                                                             x, t, surface_flux_function, equations)
    bc_callable(bc, bc_type(bc), u_inner, normal_direction, x, t, surface_flux_function, equations)
end

@inline function (bc::Union{MaxwellWallBC, OhmicContactBC})(u_inner, normal_direction::AbstractVector,
                                                             x, t, surface_flux_function, equations, 
                                                             boundary_index::Integer, node_index::Integer)
    bc_callable(bc, bc_type(bc), u_inner, normal_direction, x, t, surface_flux_function, equations,
                boundary_index, node_index)
end

# ======================================================================================================================
# calc_boundary_flux! Methods
# ======================================================================================================================

@inline function calc_boundary_flux_generic!(surface_flux_values, flux, equations, 
                                            node_index, index2, index3)
    @inbounds for v in Trixi.eachvariable(equations)
        surface_flux_values[v, node_index, index2, index3] = flux[v]
    end
    return nothing
end

@inline function Trixi.calc_boundary_flux!(surface_flux_values, t,
                                           boundary_condition::Union{MaxwellWallBC, OhmicContactBC},
                                           mesh::Trixi.UnstructuredMesh2D,
                                           have_nonconservative_terms::Trixi.False,
                                           equations::FermiHarmonics2D,
                                           surface_integral, dg::Trixi.DG, cache,
                                           node_index, side_index, element_index,
                                           boundary_index)
    u_inner = Trixi.get_node_vars(cache.boundaries.u, equations, dg, node_index, boundary_index)
    outward_direction = Trixi.get_surface_normal(cache.elements.normal_directions, node_index, 
                                                 side_index, element_index)
    x = Trixi.get_node_coords(cache.boundaries.node_coordinates, equations, dg, node_index, boundary_index)
    flux = boundary_condition(u_inner, outward_direction, x, t, surface_integral.surface_flux, 
                              equations, boundary_index, node_index)
    calc_boundary_flux_generic!(surface_flux_values, flux, equations, node_index, side_index, element_index)
end

@inline function Trixi.calc_boundary_flux!(surface_flux_values, t,
                                           boundary_condition::Union{MaxwellWallBC, OhmicContactBC},
                                           mesh::Trixi.P4estMesh{2},
                                           have_nonconservative_terms::Trixi.False,
                                           equations::FermiHarmonics2D,
                                           surface_integral, dg::Trixi.DG, cache,
                                           i_index, j_index,
                                           node_index, direction_index, element_index,
                                           boundary_index)
    u_inner = Trixi.get_node_vars(cache.boundaries.u, equations, dg, node_index, boundary_index)
    normal_direction = Trixi.get_normal_direction(direction_index, cache.elements.contravariant_vectors,
                                                  i_index, j_index, element_index)
    x = Trixi.get_node_coords(cache.elements.node_coordinates, equations, dg, i_index, j_index, element_index)
    flux = boundary_condition(u_inner, normal_direction, x, t, surface_integral.surface_flux,
                              equations, boundary_index, node_index)
    calc_boundary_flux_generic!(surface_flux_values, flux, equations, node_index, direction_index, element_index)
end

@inline function Trixi.calc_boundary_flux!(surface_flux_values, t,
                                           boundary_condition::Union{MaxwellWallBC, OhmicContactBC},
                                           mesh::Trixi.P4estMeshView{2},
                                           have_nonconservative_terms::Trixi.False,
                                           equations::FermiHarmonics2D,
                                           surface_integral, dg::Trixi.DG, cache,
                                           i_index, j_index,
                                           node_index, direction_index, element_index,
                                           boundary_index)
    return Trixi.calc_boundary_flux!(surface_flux_values, t, boundary_condition,
                                     mesh.parent, have_nonconservative_terms, equations,
                                     surface_integral, dg, cache,
                                     i_index, j_index, node_index, direction_index,
                                     element_index, boundary_index)
end
# ======================================================================================================================
# Convergence Monitoring
# ======================================================================================================================

# Infinity norm for single-node state vectors.
@inline function Trixi.residual_steady_state(du::AbstractVector, equations::FermiHarmonics2D)
    return maximum(abs, du)
end

# Infinity norm for full DG arrays.
@inline function Trixi.residual_steady_state(du::AbstractArray{<:Any, 4}, equations::FermiHarmonics2D)
    return maximum(abs, du)
end

# ======================================================================================================================
# Semidiscretize Hook
# ======================================================================================================================

function Trixi.semidiscretize(
    semi::Trixi.SemidiscretizationHyperbolic{<:Any, <:FermiHarmonics2D},
    tspan;
    kwargs...,
)
    init_projector_cache!(semi)
    return invoke(
        Trixi.semidiscretize,
        Tuple{Trixi.AbstractSemidiscretization, Any},
        semi, tspan; kwargs...
    )
end

# Restart overload.
function Trixi.semidiscretize(
    semi::Trixi.SemidiscretizationHyperbolic{<:Any, <:FermiHarmonics2D},
    tspan,
    restart_file::AbstractString;
    kwargs...,
)
    init_projector_cache!(semi)
    return invoke(
        Trixi.semidiscretize,
        Tuple{Trixi.AbstractSemidiscretization, Any, AbstractString},
        semi, tspan, restart_file; kwargs...
    )
end

# ======================================================================================================================
# Custom Analysis Integrals
# ======================================================================================================================

struct ResidualSteadyStateIntegral end

"""
    residual_steady_state_integral

Analysis-integral object for Trixi `AnalysisCallback` that reports the steady-state residual.

Returns:
- scalar residual value when used through Trixi analysis hooks.
"""
const residual_steady_state_integral = ResidualSteadyStateIntegral()

# Residual analysis hook.
function Trixi.analyze(::ResidualSteadyStateIntegral, 
                       du, u, t, semi::Trixi.AbstractSemidiscretization)
    equations = semi.equations
    return Trixi.residual_steady_state(du, equations)
end

# Pretty printing for residual integral.
Trixi.pretty_form_utf(::ResidualSteadyStateIntegral) = "residual"
Trixi.pretty_form_ascii(::ResidualSteadyStateIntegral) = "residual"

# Disable default analysis integrals/errors for this equations type.
function Trixi.default_analysis_integrals(::FermiHarmonics2D)
    return ()
end

function Trixi.default_analysis_errors(::FermiHarmonics2D)
    return Symbol[]
end

# No-op entropy time derivative.
function Trixi.entropy_timederivative(u, equations::FermiHarmonics2D)
    return 0.0
end

# Return zero state integral for analysis initialization.
function Trixi.integrate(u_ode::AbstractVector, semi::Trixi.AbstractSemidiscretization)
    equations = semi.equations
    nvars = Trixi.nvariables(equations)
    return StaticArrays.SVector{nvars, Float64}(ntuple(i -> 0.0, nvars))
end

# Scalar placeholder integral for custom analysis object.
function Trixi.integrate(func, u, mesh::Trixi.P4estMesh{2}, 
                         equations::FermiHarmonics2D, 
                         dg::Trixi.DGSEM, cache; normalize=false)
    return 0.0
end
