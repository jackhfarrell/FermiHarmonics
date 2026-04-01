function solve(args...; kwargs...)
    throw(ArgumentError("No backend solve method is available. Load a backend package such as Trixi to activate extensions."))
end

function solve_status(args...; kwargs...)
    throw(ArgumentError("solve_status is provided by backend extensions. Load Trixi to use the current backend implementation."))
end

function save_solution_custom(args...; kwargs...)
    throw(ArgumentError("save_solution_custom is provided by backend extensions. Load Trixi to use the current backend implementation."))
end

function save_for_analysis(args...; kwargs...)
    throw(ArgumentError("save_for_analysis is provided by backend extensions. Load Trixi to use the current backend implementation."))
end

function save_mesh_native_analysis(args...; kwargs...)
    throw(ArgumentError("save_mesh_native_analysis is provided by backend extensions. Load Trixi to use the current backend implementation."))
end

function compute_analysis_grids(args...; kwargs...)
    throw(ArgumentError("compute_analysis_grids is provided by backend extensions. Load Trixi to use the current backend implementation."))
end

function current_norm_variables(args...; kwargs...)
    throw(ArgumentError("current_norm_variables is provided by backend extensions. Load Trixi to use the current backend implementation."))
end

function visualization_callback(args...; kwargs...)
    throw(ArgumentError("visualization_callback is provided by backend extensions. Load Trixi to use the current backend implementation."))
end

function evaluate_solution(args...; kwargs...)
    throw(ArgumentError("evaluate_solution is provided by backend extensions. Load Trixi to use the current backend implementation."))
end

function evaluate_observables(args...; kwargs...)
    throw(ArgumentError("evaluate_observables is provided by backend extensions. Load Trixi to use the current backend implementation."))
end

function create_live_dashboard(args...; kwargs...)
    throw(ArgumentError("Live visualization requires loading GLMakie to activate the Makie extension."))
end

function update_live_dashboard!(args...; kwargs...)
    throw(ArgumentError("Live visualization requires loading GLMakie to activate the Makie extension."))
end

function finalize_live_dashboard!(args...; kwargs...)
    return nothing
end

function live_dashboard_is_open(args...; kwargs...)
    return true
end

function default_callbacks_builder(args...; kwargs...)
    throw(ArgumentError("default_callbacks_builder is provided by backend extensions. Load Trixi to use the current backend implementation."))
end
