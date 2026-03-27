"""
I/O helpers for FermiHarmonics outputs and restart compatibility.
"""

using HDF5
using NLsolve
using Trixi

# ======================================================================================================================
# Analysis Output
# ======================================================================================================================
# TODO: make save_for_analysis more flexible, we should be able to give it a list of variable names and it saves those
"""
    save_for_analysis(sol, semi, filename; nvisnodes=400, observables=nothing)

Save observables on a uniform Cartesian grid in a lightweight HDF5 format for post-processing
and analysis. The grid is determined from the simulation domain bounds.

Arguments:
- `sol`: time integration solution.
- `semi`: semidiscretization object.
- `filename`: output HDF5 path.
- `nvisnodes`: number of nodes per axis of uniform sampling grid.
- `observables`: optional list of observable names to save. `nothing` preserves the
  legacy full output. Supported nonlinear names are `:n`, `:a1`, `:b1`, `:jx`, `:jy`.

Returns:
- `filename::AbstractString`.
"""
function save_for_analysis(sol, semi, filename; nvisnodes=400, observables=nothing)
    final_time = sol.t[end]
    grids = compute_analysis_grids(sol.u[end], semi; nvisnodes=nvisnodes)
    mesh, _, _, _ = Trixi.mesh_equations_solver_cache(semi)
    @info "Analysis: writing HDF5" file=filename
    analysis_write_hdf5(filename, grids.density, grids.a1, grids.b1, grids.jx, grids.jy,
                        grids.x, grids.y, grids.mask, final_time, grids.equations;
                        band_grids=grids.bands,
                        observables=observables,
                        mesh_metadata=ElectronKinetics.mesh_provenance_attributes(mesh.current_filename))
    @info "Analysis: write complete" file=filename
    return filename
end

"""
    save_mesh_native_analysis(sol, semi, filename; refine=6, observables=nothing)

Save observables on a mesh-following unstructured visualization grid created by
subdividing each DG element in reference space.
"""
function save_mesh_native_analysis(sol, semi, filename; refine=6, observables=nothing)
    final_time = sol.t[end]
    mesh_data = compute_mesh_native_analysis(sol.u[end], semi; refine=refine)
    mesh, _, _, _ = Trixi.mesh_equations_solver_cache(semi)
    @info "Analysis: writing mesh-native HDF5" file=filename
    analysis_write_mesh_native_hdf5(filename, mesh_data, final_time, semi.equations;
                                    observables=observables,
                                    mesh_metadata=ElectronKinetics.mesh_provenance_attributes(mesh.current_filename))
    @info "Analysis: mesh-native write complete" file=filename
    return filename
end

Base.@deprecate save_observables_for_python save_for_analysis
@doc "Deprecated alias for [`save_for_analysis`](@ref)." save_observables_for_python

export save_solution_custom,
       save_for_analysis,
       save_mesh_native_analysis,
       save_observables_for_python,
       evaluate_solution,
       evaluate_observables

function normalize_analysis_observables(observables, equations)
    observables === nothing && return nothing

    normalized = Symbol[]
    seen = Set{Symbol}()
    density_name = analysis_density_name(equations)
    allowed = Set(analysis_default_observables(equations))
    if equations isa MultiBandFermiHarmonics2D
        for band in equations.bands
            push!(allowed, Symbol("$(band.name)_n"))
            push!(allowed, Symbol("$(band.name)_jx"))
            push!(allowed, Symbol("$(band.name)_jy"))
        end
    end

    for observable in observables
        name = Symbol(observable)
        if !transport_is_nonlinear(equations) && !(equations isa MultiBandFermiHarmonics2D) && name === :n
            name = :a0
        elseif equations isa MultiBandFermiHarmonics2D && name === :a0
            name = :n
        end
        name in allowed || throw(ArgumentError("unsupported analysis observable $(repr(name))"))
        if !(name in seen)
            push!(normalized, name)
            push!(seen, name)
        end
    end

    density_name in seen || throw(ArgumentError("analysis output must include $(density_name)"))
    return normalized
end

@inline analysis_density_name(equations) =
    transport_is_nonlinear(equations) || equations isa MultiBandFermiHarmonics2D ? :n : :a0

@inline function analysis_default_observables(equations)
    if transport_is_nonlinear(equations) || equations isa MultiBandFermiHarmonics2D
        return (:n, :a1, :b1, :jx, :jy)
    end
    return (:a0, :a1, :b1, :jx, :jy)
end

function analysis_grid_axes(solution_vector, semi, nvisnodes::Int)
    mesh, equations, solver, cache = Trixi.mesh_equations_solver_cache(semi)
    num_vars = Trixi.nvariables(equations)
    num_vars >= 3 || throw(ArgumentError("analysis output requires at least 3 variables (a0, a1, b1)"))
    num_nodes = Trixi.nnodes(solver)
    num_elements = Trixi.nelements(solver, cache)

    x_min = Inf
    x_max = -Inf
    y_min = Inf
    y_max = -Inf
    @inbounds for element_index in 1:num_elements, node_j in 1:num_nodes, node_i in 1:num_nodes
        coords = Trixi.get_node_coords(cache.elements.node_coordinates, equations, solver,
                                       node_i, node_j, element_index)
        x = coords[1]
        y = coords[2]
        x_min = min(x_min, x)
        x_max = max(x_max, x)
        y_min = min(y_min, y)
        y_max = max(y_max, y)
    end

    return (
        x = range(x_min, x_max, length=nvisnodes),
        y = range(y_min, y_max, length=nvisnodes),
        equations = equations,
    )
end

@inline reference_visualization_nodes(refine::Int) =
    collect(range(-1.0, 1.0, length=refine + 1))

function interpolate_element_state!(
    out::AbstractVector{Float64},
    basis_xi::AbstractVector{<:Real},
    basis_eta::AbstractVector{<:Real},
    element_u_values::AbstractArray{<:Real, 3},
)
    fill!(out, 0.0)
    nvars = size(element_u_values, 1)
    num_nodes = length(basis_xi)
    @inbounds for node_j in 1:num_nodes, node_i in 1:num_nodes
        weight = basis_xi[node_i] * basis_eta[node_j]
        for var_index in 1:nvars
            out[var_index] += weight * element_u_values[var_index, node_i, node_j]
        end
    end
    return out
end

function map_element_point(
    basis_xi::AbstractVector{<:Real},
    basis_eta::AbstractVector{<:Real},
    element_x_coords::AbstractMatrix{<:Real},
    element_y_coords::AbstractMatrix{<:Real},
)
    x_mapped = 0.0
    y_mapped = 0.0
    num_nodes = length(basis_xi)
    @inbounds for node_j in 1:num_nodes, node_i in 1:num_nodes
        weight = basis_xi[node_i] * basis_eta[node_j]
        x_mapped += weight * element_x_coords[node_i, node_j]
        y_mapped += weight * element_y_coords[node_i, node_j]
    end
    return x_mapped, y_mapped
end

function compute_mesh_native_analysis(solution_vector, semi; refine=6)
    refine_int = Int(refine)
    refine_int >= 1 || throw(ArgumentError("refine must be >= 1"))

    mesh, equations, solver, cache = Trixi.mesh_equations_solver_cache(semi)
    num_vars = Trixi.nvariables(equations)
    num_nodes = Trixi.nnodes(solver)
    num_elements = Trixi.nelements(solver, cache)
    solution_wrapped = Trixi.wrap_array(solution_vector, semi)
    basis_nodes = solver.basis.nodes
    visual_nodes = reference_visualization_nodes(refine_int)
    side_points = refine_int + 1
    points_per_element = side_points^2
    num_points = num_elements * points_per_element
    num_triangles = num_elements * 2 * refine_int^2

    x_points = Vector{Float64}(undef, num_points)
    y_points = Vector{Float64}(undef, num_points)
    density_points = Vector{Float64}(undef, num_points)
    a1_points = Vector{Float64}(undef, num_points)
    b1_points = Vector{Float64}(undef, num_points)
    jx_points = Vector{Float64}(undef, num_points)
    jy_points = Vector{Float64}(undef, num_points)
    band_points = if equations isa MultiBandFermiHarmonics2D
        Dict(
            band.name => (
                n = Vector{Float64}(undef, num_points),
                jx = Vector{Float64}(undef, num_points),
                jy = Vector{Float64}(undef, num_points),
            ) for band in equations.bands
        )
    else
        Dict{Symbol, Any}()
    end
    triangles = Matrix{Int32}(undef, num_triangles, 3)

    basis_cache = [lagrange_basis(basis_nodes, xi) for xi in visual_nodes]
    state_buffer = zeros(Float64, num_vars)
    triangle_index = 1

    @inbounds for element_index in 1:num_elements
        element_x_coords = zeros(num_nodes, num_nodes)
        element_y_coords = zeros(num_nodes, num_nodes)
        element_u_values = zeros(num_vars, num_nodes, num_nodes)
        for node_j in 1:num_nodes, node_i in 1:num_nodes
            coords = Trixi.get_node_coords(cache.elements.node_coordinates, equations, solver,
                                           node_i, node_j, element_index)
            element_x_coords[node_i, node_j] = coords[1]
            element_y_coords[node_i, node_j] = coords[2]
            vars = Trixi.get_node_vars(solution_wrapped, equations, solver,
                                       node_i, node_j, element_index)
            for var_index in 1:num_vars
                element_u_values[var_index, node_i, node_j] = vars[var_index]
            end
        end

        point_base = (element_index - 1) * points_per_element
        for eta_index in 1:side_points, xi_index in 1:side_points
            basis_xi = basis_cache[xi_index]
            basis_eta = basis_cache[eta_index]
            point_index = point_base + (eta_index - 1) * side_points + xi_index

            x_mapped, y_mapped = map_element_point(
                basis_xi, basis_eta, element_x_coords, element_y_coords,
            )
            interpolate_element_state!(state_buffer, basis_xi, basis_eta, element_u_values)

            if transport_is_nonlinear(equations)
                if equations isa FermiAngles2D
                    density_value, jx_value, jy_value, _ = anglegrid_moments(state_buffer, equations)
                    _, a1_value, b1_value = derived_harmonics(state_buffer, equations)
                else
                    density_value = nonlinear_density(state_buffer, equations)
                    _, a1_value, b1_value = derived_harmonics(state_buffer, equations)
                    jx_value, jy_value = nonlinear_current(state_buffer, equations)
                end
                a1_points[point_index] = a1_value
                b1_points[point_index] = b1_value
                jx_points[point_index] = jx_value
                jy_points[point_index] = jy_value
            elseif equations isa MultiBandFermiHarmonics2D
                obs = multiband_observables(state_buffer, equations)
                density_value = obs.n
                a1_points[point_index] = obs.a1
                b1_points[point_index] = obs.b1
                jx_points[point_index] = obs.jx
                jy_points[point_index] = obs.jy
                for band in equations.bands
                    band_obs = getproperty(obs.bands, band.name)
                    band_points[band.name].n[point_index] = band_obs.n
                    band_points[band.name].jx[point_index] = band_obs.jx
                    band_points[band.name].jy[point_index] = band_obs.jy
                end
            else
                density_value = state_buffer[1]
                a1_value = length(state_buffer) >= 2 ? state_buffer[2] : 0.0
                b1_value = length(state_buffer) >= 3 ? state_buffer[3] : 0.0
                a1_points[point_index] = a1_value
                b1_points[point_index] = b1_value
                jx_points[point_index] = a1_value
                jy_points[point_index] = b1_value
            end

            x_points[point_index] = x_mapped
            y_points[point_index] = y_mapped
            density_points[point_index] = density_value
        end

        for cell_eta in 1:refine_int, cell_xi in 1:refine_int
            lower_left = point_base + (cell_eta - 1) * side_points + cell_xi
            lower_right = lower_left + 1
            upper_left = lower_left + side_points
            upper_right = upper_left + 1

            triangles[triangle_index, 1] = Int32(lower_left - 1)
            triangles[triangle_index, 2] = Int32(lower_right - 1)
            triangles[triangle_index, 3] = Int32(upper_right - 1)
            triangle_index += 1
            triangles[triangle_index, 1] = Int32(lower_left - 1)
            triangles[triangle_index, 2] = Int32(upper_right - 1)
            triangles[triangle_index, 3] = Int32(upper_left - 1)
            triangle_index += 1
        end
    end

    return (
        x = x_points,
        y = y_points,
        n = density_points,
        a1 = a1_points,
        b1 = b1_points,
        jx = jx_points,
        jy = jy_points,
        bands = band_points,
        triangles = triangles,
        refine = refine_int,
    )
end

function compute_analysis_grids(solution_vector, semi; nvisnodes=400, log::Bool=true)
    nvisnodes = Int(nvisnodes)
    axes = analysis_grid_axes(solution_vector, semi, nvisnodes)
    equations = axes.equations
    x_uniform = axes.x
    y_uniform = axes.y

    log && @info "Analysis: direct grid evaluation" nvisnodes
    num_x = length(x_uniform)
    num_y = length(y_uniform)
    density_grid = fill(NaN, num_x, num_y)
    a1_grid = fill(NaN, num_x, num_y)
    b1_grid = fill(NaN, num_x, num_y)
    jx_grid = fill(NaN, num_x, num_y)
    jy_grid = fill(NaN, num_x, num_y)
    band_grids = if equations isa MultiBandFermiHarmonics2D
        Dict(
            band.name => (
                n = fill(NaN, num_x, num_y),
                jx = fill(NaN, num_x, num_y),
                jy = fill(NaN, num_x, num_y),
            ) for band in equations.bands
        )
    else
        Dict{Symbol, Any}()
    end
    in_domain_mask = fill(false, num_x, num_y)
    @inbounds for y_index in 1:num_y
        for x_index in 1:num_x
            x_target = x_uniform[x_index]
            y_target = y_uniform[y_index]
            density_value, a1_value, b1_value, jx_value, jy_value, band_values, in_domain =
                evaluate_analysis_observables(solution_vector, semi, x_target, y_target)
            density_grid[x_index, y_index] = density_value
            a1_grid[x_index, y_index] = a1_value
            b1_grid[x_index, y_index] = b1_value
            jx_grid[x_index, y_index] = jx_value
            jy_grid[x_index, y_index] = jy_value
            if equations isa MultiBandFermiHarmonics2D
                for band in equations.bands
                    band_obs = getproperty(band_values, band.name)
                    band_grids[band.name].n[x_index, y_index] = band_obs.n
                    band_grids[band.name].jx[x_index, y_index] = band_obs.jx
                    band_grids[band.name].jy[x_index, y_index] = band_obs.jy
                end
            end
            in_domain_mask[x_index, y_index] = in_domain
        end
    end

    return (
        density = density_grid,
        n = density_grid,
        a0 = density_grid,
        a1 = a1_grid,
        b1 = b1_grid,
        jx = jx_grid,
        jy = jy_grid,
        bands = band_grids,
        x = x_uniform,
        y = y_uniform,
        mask = in_domain_mask,
        equations = equations,
    )
end

"""
    load_restart_compatible(filename, semi)

Read a restart HDF5 file and return a solution vector compatible with the current
semidiscretization described by `semi`. If the saved file has fewer variables than
the current semidiscretization the remaining variables are zero-padded. If the
saved file has more variables they are truncated (highest indices dropped).
If the saved file stored `source_index` attributes (written by `save_solution_custom`
when saving a subset) those indices are honoured.
"""
function load_restart_compatible(filename::AbstractString, semi)
    mesh, equations, solver, cache = Trixi.mesh_equations_solver_cache(semi)
    target_nvars = Trixi.nvariables(equations)
    num_elements = Trixi.nelements(solver, cache)
    num_nodes = Trixi.nnodes(solver)
    block = num_elements * num_nodes^2

    # Target array in (nvars, block) layout.
    target_array = zeros(Float64, target_nvars, block)

    h5open(filename, "r") do f
        # Number of saved variable datasets.
        saved_nvars = try
            Int(attributes(f)["n_vars"])
        catch
            # Fallback: count variables_* datasets.
            count = 0
            for name in keys(f)
                startswith(String(name), "variables_") && (count += 1)
            end
            count
        end

        # Read saved variables into target layout.
        for i in 1:saved_nvars
            dname = "variables_$i"
            if haskey(f, dname)
                data = read(f[dname])
                vecdata = vec(data)
                length(vecdata) != block && warn("saved variable size does not match expected block size: $dname")

                # Explicit source index for subset saves.
                src_idx = try
                    Int(get(attributes(f[dname]), "source_index", nothing))
                catch
                    nothing
                end

                if src_idx !== nothing && 1 <= src_idx <= target_nvars
                    target_array[src_idx, 1:length(vecdata)] .= vecdata[1:min(end, block)]
                else
                    # Sequential placement when source index is not present.
                    if i <= target_nvars
                        target_array[i, 1:length(vecdata)] .= vecdata[1:min(end, block)]
                    end
                end
            end
        end
    end

    # Flatten in variables-major ordering expected by Trixi.
    return vec(target_array)
end




"""
    save_solution_custom(sol, semi, filename; variable_names=nothing)

Save solution to HDF5 file with custom filename. Restart-compatible with Trixi.
Always saves conservative variables (no transformations). Use `variable_names`
to save a subset for lightweight analysis (not restart-compatible).

Arguments:
- `sol`: time integration solution.
- `semi`: semidiscretization object.
- `filename`: output HDF5 path.
- `variable_names`: optional subset of variable names to save.

Returns:
- `filename::AbstractString`.
"""
function save_solution_custom(sol, semi, filename; variable_names=nothing)
    # Extract solution and metadata from Trixi containers.
    mesh, equations, solver, cache = Trixi.mesh_equations_solver_cache(semi)
    solution_vector = sol.u[end]
    final_time = sol.t[end]
    solution_array = Trixi.wrap_array_native(solution_vector, mesh, equations, solver, cache)
    num_vars = Trixi.nvariables(equations)
    all_variable_names = collect(Trixi.varnames(Trixi.cons2cons, equations))
    is_subset = variable_names !== nothing
    variable_indices = Int[]
    if is_subset
        for variable_name in variable_names
            variable_index = findfirst(==(variable_name), all_variable_names)
            variable_index === nothing && error("Variable name not found: $variable_name")
            push!(variable_indices, variable_index)
        end
    else
        variable_indices = collect(1:num_vars)
    end
    
    # Write restart-compatible HDF5.
    h5open(filename, "w") do file
        # Minimal attributes for Trixi restart.
        attributes(file)["ndims"] = Trixi.ndims(mesh)
        attributes(file)["equations"] = Trixi.get_name(equations)
        attributes(file)["polydeg"] = Trixi.polydeg(solver)
        attributes(file)["n_vars"] = length(variable_indices)
        attributes(file)["n_elements"] = Trixi.nelements(solver, cache)
        attributes(file)["mesh_type"] = Trixi.get_name(mesh)
        for (key, value) in ElectronKinetics.mesh_provenance_attributes(mesh.current_filename)
            attributes(file)[key] = value
        end
        attributes(file)["time"] = Float64(final_time)
        attributes(file)["dt"] = 0.0
        attributes(file)["timestep"] = 0
        if is_subset
            attributes(file)["subset_of_n_vars"] = num_vars
            attributes(file)["subset_names"] = join(variable_names, ",")
        end
        # Conservative variable datasets.
        for (output_index, variable_index) in enumerate(variable_indices)
            file["variables_$output_index"] = vec(solution_array[variable_index, .., :])
            attributes(file["variables_$output_index"])["name"] = all_variable_names[variable_index]
            if is_subset
                attributes(file["variables_$output_index"])["source_index"] = variable_index
            end
        end
    end
    return filename
end

# Evaluate Lagrange basis functions at reference coordinate `xi`.
@inline function lagrange_basis(nodes, xi)
    n = length(nodes)
    basis = zeros(n)
    @inbounds for i in 1:n
        li = 1.0
        for j in 1:n
            if j != i
                li *= (xi - nodes[j]) / (nodes[i] - nodes[j])
            end
        end
        basis[i] = li
    end
    return basis
end


function interpolate_state_at_point(
    solution_vector,
    semi, 
    x_target, 
    y_target;
    max_newton::Int=10, tol::Float64=1e-12
)
    mesh, equations, solver, cache = Trixi.mesh_equations_solver_cache(semi)
    num_vars = Trixi.nvariables(equations)
    num_nodes = Trixi.nnodes(solver)
    num_elements = Trixi.nelements(solver, cache)
    nodes = solver.basis.nodes
    solution_wrapped = Trixi.wrap_array(solution_vector, semi)
    
    # Search candidate element and invert mapping.
    @inbounds for element_index in 1:num_elements
        element_x_coords = zeros(num_nodes, num_nodes)
        element_y_coords = zeros(num_nodes, num_nodes)
        element_u_values = zeros(num_vars, num_nodes, num_nodes)
        
        # Gather element coordinates and values.
        for node_j in 1:num_nodes, node_i in 1:num_nodes
            coords = Trixi.get_node_coords(cache.elements.node_coordinates, equations, solver,
                                           node_i, node_j, element_index)
            element_x_coords[node_i, node_j] = coords[1]
            element_y_coords[node_i, node_j] = coords[2]
            vars = Trixi.get_node_vars(solution_wrapped, equations, solver, node_i, node_j, element_index)
            for var_index in 1:num_vars
                element_u_values[var_index, node_i, node_j] = vars[var_index]
            end
        end

        # Fast bounding-box reject.
        x_min_element, x_max_element = extrema(element_x_coords)
        y_min_element, y_max_element = extrema(element_y_coords)
        if x_target < x_min_element - tol || x_target > x_max_element + tol ||
           y_target < y_min_element - tol || y_target > y_max_element + tol
            continue
        end

        # Residual for reference-coordinate inversion.
        function reference_residual!(residual, reference_coords)
            xi = reference_coords[1]
            eta = reference_coords[2]
            basis_xi = lagrange_basis(nodes, xi)
            basis_eta = lagrange_basis(nodes, eta)
            x_mapped = 0.0
            y_mapped = 0.0
            for node_j in 1:num_nodes, node_i in 1:num_nodes
                weight = basis_xi[node_i] * basis_eta[node_j]
                x_mapped += weight * element_x_coords[node_i, node_j]
                y_mapped += weight * element_y_coords[node_i, node_j]
            end
            residual[1] = x_mapped - x_target
            residual[2] = y_mapped - y_target
            return residual
        end

        # Solve for reference coordinates.
        initial_guess = [0.0, 0.0]
        result = nlsolve(
            reference_residual!, 
            initial_guess;
            method=:newton, 
            ftol=tol, 
            xtol=tol, 
            iterations=max_newton)
        
        # Interpolate if inversion converged and point is inside reference element.
        if converged(result)
            xi = result.zero[1]
            eta = result.zero[2]
            if abs(xi) <= 1.0 + 1e-10 && abs(eta) <= 1.0 + 1e-10
                basis_xi = lagrange_basis(nodes, xi)
                basis_eta = lagrange_basis(nodes, eta)
                state_value = zeros(Float64, num_vars)
                for node_j in 1:num_nodes, node_i in 1:num_nodes
                    weight = basis_xi[node_i] * basis_eta[node_j]
                    for var_index in 1:num_vars
                        state_value[var_index] += weight * element_u_values[var_index, node_i, node_j]
                    end
                end

                return state_value, true
            end
        end
    end
    
    # Point not in domain.
    return Float64[], false
end

"""
    evaluate_solution(sol, semi, x_target, y_target; max_newton=10, tol=1e-12)

Evaluate `(a0, a1, b1, in_domain)` at one Cartesian point by element search and
reference-coordinate solve.
"""
function evaluate_solution(sol, semi, x_target, y_target; max_newton::Int=10, tol::Float64=1e-12)
    state_value, in_domain = interpolate_state_at_point(
        sol.u[end], semi, x_target, y_target; max_newton=max_newton, tol=tol,
    )
    if !in_domain
        return NaN, NaN, NaN, false
    end
    if transport_is_nonlinear(semi.equations) || semi.equations isa MultiBandFermiHarmonics2D
        a0_value, a1_value, b1_value = derived_harmonics(state_value, semi.equations)
        return a0_value, a1_value, b1_value, true
    end
    a0_value = state_value[1]
    a1_value = length(state_value) >= 2 ? state_value[2] : 0.0
    b1_value = length(state_value) >= 3 ? state_value[3] : 0.0
    return a0_value, a1_value, b1_value, true
end

"""
    evaluate_observables(sol, semi, x_target, y_target; max_newton=10, tol=1e-12)

Evaluate the physical observables at one Cartesian point.

Returns a named tuple with:
- `n`
- `a1`
- `b1`
- `jx`
- `jy`
- `in_domain`

For linear transport, `jx == a1` and `jy == b1`.
"""
function evaluate_observables(sol, semi, x_target, y_target; max_newton::Int=10, tol::Float64=1e-12)
    density_value, a1_value, b1_value, jx_value, jy_value, band_values, in_domain = evaluate_analysis_observables(
        sol.u[end], semi, x_target, y_target; max_newton=max_newton, tol=tol,
    )
    return (
        n = density_value,
        a0 = density_value,
        a1 = a1_value,
        b1 = b1_value,
        jx = jx_value,
        jy = jy_value,
        bands = band_values,
        in_domain = in_domain,
    )
end

function evaluate_analysis_observables(solution_vector, semi, x_target, y_target; max_newton::Int=10, tol::Float64=1e-12)
    state_value, in_domain = interpolate_state_at_point(
        solution_vector, semi, x_target, y_target; max_newton=max_newton, tol=tol,
    )
    if !in_domain
        return NaN, NaN, NaN, NaN, NaN, (;), false
    end

    equations = semi.equations
    if transport_is_nonlinear(equations)
        if equations isa FermiAngles2D
            density_value, jx_value, jy_value, _ = anglegrid_moments(state_value, equations)
            _, a1_value, b1_value = derived_harmonics(state_value, equations)
        else
            density_value = nonlinear_density(state_value, equations)
            _, a1_value, b1_value = derived_harmonics(state_value, equations)
            jx_value, jy_value = nonlinear_current(state_value, equations)
        end
        return density_value, a1_value, b1_value, jx_value, jy_value, (;), true
    elseif equations isa MultiBandFermiHarmonics2D
        obs = multiband_observables(state_value, equations)
        return obs.n, obs.a1, obs.b1, obs.jx, obs.jy, obs.bands, true
    end

    density_value = state_value[1]
    a1_value = length(state_value) >= 2 ? state_value[2] : 0.0
    b1_value = length(state_value) >= 3 ? state_value[3] : 0.0
    return density_value, a1_value, b1_value, a1_value, b1_value, (;), true
end

function analysis_write_hdf5(filename, density_grid, a1_grid, b1_grid, jx_grid, jy_grid, x_uniform, y_uniform,
                              in_domain_mask, t, equations;
                              observables=nothing,
                              band_grids=Dict{Symbol, Any}(),
                              mesh_metadata=Dict{String, Any}())
    requested = normalize_analysis_observables(observables, equations)
    density_name = String(analysis_density_name(equations))
    h5open(filename, "w") do file
        if isnothing(requested) || Symbol(density_name) in requested
            file[density_name] = density_grid
        end
        if isnothing(requested) || :a1 in requested
            file["a1"] = a1_grid
        end
        if isnothing(requested) || :b1 in requested
            file["b1"] = b1_grid
        end
        if !isnothing(jx_grid) && (isnothing(requested) || :jx in requested)
            file["jx"] = jx_grid
        end
        if !isnothing(jy_grid) && (isnothing(requested) || :jy in requested)
            file["jy"] = jy_grid
        end
        if equations isa MultiBandFermiHarmonics2D
            for band in equations.bands
                band_name = band.name
                band_n_name = Symbol("$(band_name)_n")
                band_jx_name = Symbol("$(band_name)_jx")
                band_jy_name = Symbol("$(band_name)_jy")
                band_data = band_grids[band_name]
                if isnothing(requested) || band_n_name in requested
                    file[string(band_n_name)] = band_data.n
                end
                if isnothing(requested) || band_jx_name in requested
                    file[string(band_jx_name)] = band_data.jx
                end
                if isnothing(requested) || band_jy_name in requested
                    file[string(band_jy_name)] = band_data.jy
                end
            end
        end
        file["x"] = collect(x_uniform)
        file["y"] = collect(y_uniform)
        file["mask"] = collect(in_domain_mask)

        for (key, value) in mesh_metadata
            attributes(file)[key] = value
        end
        attributes(file)["time"] = Float64(t)
        attributes(file)["nx"] = length(x_uniform)
        attributes(file)["ny"] = length(y_uniform)
        attributes(file)["grid_type"] = "uniform_cartesian"
        attributes(file)["mask_method"] = "direct"
        attributes(file)["saved_observables"] = isnothing(requested) ?
            join(string.(analysis_default_observables(equations)), ",") :
            join(string.(requested), ",")
        if transport_is_nonlinear(equations)
            attributes(file)["description"] = isnothing(requested) ?
                "Nonlinear observables: density n, currents jx and jy, plus harmonic reference fields a1 and b1" :
                "Selected nonlinear observables on a uniform Cartesian grid"
        elseif equations isa MultiBandFermiHarmonics2D
            attributes(file)["description"] = isnothing(requested) ?
                "Linear multiband observables: total density n, total currents jx and jy, plus total current-like reference fields a1 and b1" :
                "Selected linear multiband observables on a uniform Cartesian grid"
        else
            attributes(file)["description"] = isnothing(requested) ?
                "Observable harmonics: a0 (density), a1 (x-current), b1 (y-current)" :
                "Selected linear observables on a uniform Cartesian grid"
        end
    end
end

function analysis_write_mesh_native_hdf5(filename, mesh_data, t, equations;
                                         observables=nothing,
                                         mesh_metadata=Dict{String, Any}())
    requested = normalize_analysis_observables(observables, equations)
    density_name = String(analysis_density_name(equations))
    h5open(filename, "w") do file
        file["x"] = mesh_data.x
        file["y"] = mesh_data.y
        file["triangles"] = mesh_data.triangles
        if isnothing(requested) || Symbol(density_name) in requested
            file[density_name] = mesh_data.n
        end
        if isnothing(requested) || :a1 in requested
            file["a1"] = mesh_data.a1
        end
        if isnothing(requested) || :b1 in requested
            file["b1"] = mesh_data.b1
        end
        if isnothing(requested) || :jx in requested
            file["jx"] = mesh_data.jx
        end
        if isnothing(requested) || :jy in requested
            file["jy"] = mesh_data.jy
        end
        if equations isa MultiBandFermiHarmonics2D
            for band in equations.bands
                band_name = band.name
                band_n_name = Symbol("$(band_name)_n")
                band_jx_name = Symbol("$(band_name)_jx")
                band_jy_name = Symbol("$(band_name)_jy")
                band_data = mesh_data.bands[band_name]
                if isnothing(requested) || band_n_name in requested
                    file[string(band_n_name)] = band_data.n
                end
                if isnothing(requested) || band_jx_name in requested
                    file[string(band_jx_name)] = band_data.jx
                end
                if isnothing(requested) || band_jy_name in requested
                    file[string(band_jy_name)] = band_data.jy
                end
            end
        end

        for (key, value) in mesh_metadata
            attributes(file)[key] = value
        end
        attributes(file)["time"] = Float64(t)
        attributes(file)["grid_type"] = "mesh_native_triangles"
        attributes(file)["refine"] = mesh_data.refine
        attributes(file)["connectivity_index_base"] = 0
        attributes(file)["saved_observables"] = isnothing(requested) ?
            join(string.(analysis_default_observables(equations)), ",") :
            join(string.(requested), ",")
        attributes(file)["description"] = transport_is_nonlinear(equations) ?
            (isnothing(requested) ?
             "Mesh-native nonlinear observables on a refined unstructured visualization grid" :
             "Selected mesh-native nonlinear observables on a refined unstructured visualization grid") :
            equations isa MultiBandFermiHarmonics2D ?
            (isnothing(requested) ?
             "Mesh-native linear multiband observables on a refined unstructured visualization grid" :
             "Selected mesh-native linear multiband observables on a refined unstructured visualization grid") :
            (isnothing(requested) ?
             "Mesh-native linear observables on a refined unstructured visualization grid" :
             "Selected mesh-native linear observables on a refined unstructured visualization grid")
    end
end
