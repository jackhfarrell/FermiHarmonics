Base.@kwdef struct LiveVisualizationConfig
    field::Union{Nothing, Symbol} = nothing
    colormap::Symbol = :magma
    geometry_mode::Symbol = :mesh_native
    nvisnodes::Int = 120
    refine::Int = 4
    accepted_step_interval::Int = 500
    min_update_seconds::Float64 = 0.5
    show_window::Bool = true
end

function validate(config::LiveVisualizationConfig)
    config.geometry_mode in (:cartesian, :mesh_native) ||
        throw(ArgumentError("live visualization geometry_mode must be :cartesian or :mesh_native"))
    config.nvisnodes >= 8 || throw(ArgumentError("live visualization nvisnodes must be >= 8"))
    config.refine >= 1 || throw(ArgumentError("live visualization refine must be >= 1"))
    config.accepted_step_interval >= 1 ||
        throw(ArgumentError("live visualization accepted_step_interval must be >= 1"))
    config.min_update_seconds >= 0.0 ||
        throw(ArgumentError("live visualization min_update_seconds must be >= 0"))
    return config
end

struct LiveProgressSnapshot
    accepted_steps::Int
    current_time::Float64
    final_time::Float64
    residual::Float64
    residual_tol::Float64
    residual_progress::Float64
    time_progress::Float64
    leading_stop_condition::Symbol
    stop_reason::Symbol
end

struct LiveFieldSnapshot
    field::Symbol
    geometry_mode::Symbol
    label::String
    x
    y
    values
    mask
    triangles
end

struct LiveVisualizationSnapshot
    progress::LiveProgressSnapshot
    field::LiveFieldSnapshot
end

function residual_progress_fraction(initial_residual::Real, residual::Real, residual_tol::Real)
    residual_tol > 0 || return 0.0
    residual <= residual_tol && return 1.0
    initial_residual <= residual_tol && return 0.0
    initial_ratio = max(Float64(initial_residual) / Float64(residual_tol), 1.0)
    current_ratio = max(Float64(residual) / Float64(residual_tol), 1.0)
    initial_ratio == 1.0 && return 0.0
    fraction = 1.0 - log10(current_ratio) / log10(initial_ratio)
    return clamp(fraction, 0.0, 1.0)
end
