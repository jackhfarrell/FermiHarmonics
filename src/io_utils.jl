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
    save_for_analysis(sol, semi, filename; nvisnodes=400)

Save ``a0``, ``a1``, ``b1`` observables on a uniform Cartesian grid in a lightweight HDF5 format
for post-processing and analysis.  The grid is determined from the simulation domain bounds.

Arguments:
- `sol`: time integration solution.
- `semi`: semidiscretization object.
- `filename`: output HDF5 path.
- `nvisnodes`: number of nodes per axis of uniform sampling grid.

Returns:
- `filename::AbstractString`.
"""
function save_for_analysis(sol, semi, filename; nvisnodes=400)
    solution_vector = sol.u[end]
    final_time = sol.t[end]
    nvisnodes = Int(nvisnodes)
    mesh, equations, solver, cache = Trixi.mesh_equations_solver_cache(semi)
    num_vars = Trixi.nvariables(equations)
    num_vars >= 3 || throw(ArgumentError("save_for_analysis requires at least 3 variables (a0, a1, b1)"))
    num_nodes = Trixi.nnodes(solver)
    num_elements = Trixi.nelements(solver, cache)
    
    # Determine domain bounds from DG nodes.
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

    # Uniform Cartesian evaluation grid.
    x_uniform = range(x_min, x_max, length=nvisnodes)
    y_uniform = range(y_min, y_max, length=nvisnodes)

    @info "Analysis: direct grid evaluation" nvisnodes
    num_x = length(x_uniform)
    num_y = length(y_uniform)
    a0_grid = fill(NaN, num_x, num_y)
    a1_grid = fill(NaN, num_x, num_y)
    b1_grid = fill(NaN, num_x, num_y)
    in_domain_mask = fill(false, num_x, num_y)

    @inbounds for y_index in 1:num_y, x_index in 1:num_x
        x_target = x_uniform[x_index]
        y_target = y_uniform[y_index]
        a0_value, a1_value, b1_value, in_domain = evaluate_solution(sol, semi, x_target, y_target)
        a0_grid[x_index, y_index] = a0_value
        a1_grid[x_index, y_index] = a1_value
        b1_grid[x_index, y_index] = b1_value
        in_domain_mask[x_index, y_index] = in_domain
    end
    @info "Analysis: writing HDF5" file=filename
    analysis_write_hdf5(filename, a0_grid, a1_grid, b1_grid, x_uniform, y_uniform,
                        in_domain_mask, final_time)
    @info "Analysis: write complete" file=filename
    return filename
end

Base.@deprecate save_observables_for_python save_for_analysis
@doc "Deprecated alias for [`save_for_analysis`](@ref)." save_observables_for_python

export save_solution_custom, save_for_analysis, save_observables_for_python, evaluate_solution

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
        attributes(file)["mesh_file"] = splitdir(mesh.current_filename)[2]
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


"""
    evaluate_solution(sol, semi, x_target, y_target; max_newton=10, tol=1e-12)

Evaluate `(a0, a1, b1, in_domain)` at one Cartesian point by element search and
reference-coordinate solve.
"""
function evaluate_solution(
    sol, 
    semi, 
    x_target, 
    y_target;
    max_newton::Int=10, tol::Float64=1e-12
)
    # Solution and mesh data.
    solution_vector = sol.u[end]
    mesh, equations, solver, cache = Trixi.mesh_equations_solver_cache(semi)
    num_nodes = Trixi.nnodes(solver)
    num_elements = Trixi.nelements(solver, cache)
    nodes = solver.basis.nodes
    solution_wrapped = Trixi.wrap_array(solution_vector, semi)
    
    # Search candidate element and invert mapping.
    @inbounds for element_index in 1:num_elements
        element_x_coords = zeros(num_nodes, num_nodes)
        element_y_coords = zeros(num_nodes, num_nodes)
        element_u_values = zeros(3, num_nodes, num_nodes)
        
        # Gather element coordinates and values.
        for node_j in 1:num_nodes, node_i in 1:num_nodes
            coords = Trixi.get_node_coords(cache.elements.node_coordinates, equations, solver,
                                           node_i, node_j, element_index)
            element_x_coords[node_i, node_j] = coords[1]
            element_y_coords[node_i, node_j] = coords[2]
            vars = Trixi.get_node_vars(solution_wrapped, equations, solver, node_i, node_j, element_index)
            element_u_values[1, node_i, node_j] = vars[1]
            element_u_values[2, node_i, node_j] = vars[2]
            element_u_values[3, node_i, node_j] = vars[3]
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
                a0_value = 0.0
                a1_value = 0.0
                b1_value = 0.0
                for node_j in 1:num_nodes, node_i in 1:num_nodes
                    weight = basis_xi[node_i] * basis_eta[node_j]
                    a0_value += weight * element_u_values[1, node_i, node_j]
                    a1_value += weight * element_u_values[2, node_i, node_j]
                    b1_value += weight * element_u_values[3, node_i, node_j]
                end

                return a0_value, a1_value, b1_value, true
            end
        end
    end
    
    # Point not in domain.
    return NaN, NaN, NaN, false
end

function analysis_write_hdf5(filename, a0_grid, a1_grid, b1_grid, x_uniform, y_uniform,
                              in_domain_mask, t)
    h5open(filename, "w") do file
        # a0_grid[i, j] is at (x[i], y[j])
        file["a0"] = a0_grid
        file["a1"] = a1_grid
        file["b1"] = b1_grid
        file["x"] = collect(x_uniform)
        file["y"] = collect(y_uniform)
        file["mask"] = collect(in_domain_mask)

        attributes(file)["time"] = Float64(t)
        attributes(file)["nx"] = length(x_uniform)
        attributes(file)["ny"] = length(y_uniform)
        attributes(file)["grid_type"] = "uniform_cartesian"
        attributes(file)["mask_method"] = "direct"
        attributes(file)["description"] = "Observable harmonics: a0 (density), a1 (x-current), b1 (y-current)"
    end
end
