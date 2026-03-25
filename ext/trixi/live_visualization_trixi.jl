mutable struct SolveMonitorState
    live_visualization::Union{Nothing, LiveVisualizationConfig}
    dashboard
    initial_residual::Float64
    last_update_time::Float64
    last_accepted_steps::Int
    window_closed::Bool
end

@inline function has_current_contact_boundary(semi)
    for bc in semi.boundary_conditions.boundary_condition_types
        if bc isa CurrentContactBC
            return true
        end
    end
    return false
end

function current_contact_target_summary(semi)
    targets = Float64[]
    for bc in semi.boundary_conditions.boundary_condition_types
        if bc isa CurrentContactBC
            push!(targets, bc.target_outward_flux)
        end
    end
    isempty(targets) && return nothing
    return (
        sum = sum(targets),
        abs_sum = sum(abs, targets),
        max_abs = maximum(abs, targets),
        count = length(targets),
    )
end

function trapz_line(x_values::AbstractVector{<:Real}, y_values::AbstractVector{<:Real})
    length(x_values) == length(y_values) || return 0.0
    length(x_values) >= 2 || return 0.0
    total = 0.0
    @inbounds for i in 1:(length(x_values) - 1)
        total += 0.5 * (Float64(x_values[i + 1]) - Float64(x_values[i])) *
                 (Float64(y_values[i + 1]) + Float64(y_values[i]))
    end
    return total
end

function measure_midline_jy(solution_vector, semi; target_y::Float64=0.0, nvisnodes::Int=181)
    grids = compute_analysis_grids(solution_vector, semi; nvisnodes=nvisnodes, log=false)
    y_index = argmin(abs.(grids.y .- target_y))
    line_mask = vec(grids.mask[:, y_index])
    x_line = collect(grids.x[line_mask])
    jy_line = collect(grids.jy[line_mask, y_index])
    length(x_line) >= 2 || return (average=0.0, integrated=0.0, y=Float64(grids.y[y_index]), points=length(x_line))
    order = sortperm(x_line)
    x_sorted = x_line[order]
    jy_sorted = jy_line[order]
    integrated = trapz_line(x_sorted, jy_sorted)
    span = Float64(maximum(x_sorted) - minimum(x_sorted))
    average = iszero(span) ? Float64(jy_sorted[1]) : integrated / span
    return (average=average, integrated=integrated, y=Float64(grids.y[y_index]), points=length(x_sorted))
end

function default_live_field(equations)
    if transport_is_nonlinear(equations)
        return :current_magnitude
    elseif equations isa MultiBandFermiHarmonics2D
        return :n
    end
    return :a0
end

function field_label(field::Symbol)
    field === :a0 && return "a0"
    field === :n && return "n"
    field === :a1 && return "a1"
    field === :b1 && return "b1"
    field === :jx && return "jx"
    field === :jy && return "jy"
    field === :current_magnitude && return "|j|"
    return replace(String(field), "_" => " ")
end

function select_live_field(snapshot_source, field::Symbol)
    if field === :current_magnitude
        return sqrt.(snapshot_source.jx .^ 2 .+ snapshot_source.jy .^ 2)
    elseif field === :a0 || field === :n
        return hasproperty(snapshot_source, :n) ? snapshot_source.n : snapshot_source.a0
    elseif field === :a1
        return snapshot_source.a1
    elseif field === :b1
        return snapshot_source.b1
    elseif field === :jx
        return snapshot_source.jx
    elseif field === :jy
        return snapshot_source.jy
    elseif hasproperty(snapshot_source, :bands)
        name = String(field)
        for suffix in ("_n", "_jx", "_jy")
            if endswith(name, suffix)
                band_name = Symbol(first(name, lastindex(name) - length(suffix)))
                band_data = get(snapshot_source.bands, band_name, nothing)
                isnothing(band_data) && break
                suffix === "_n" && return band_data.n
                suffix === "_jx" && return band_data.jx
                suffix === "_jy" && return band_data.jy
            end
        end
    end
    throw(ArgumentError("unsupported live visualization field $(repr(field))"))
end

function build_live_field_snapshot(solution_vector, semi, live_visualization::LiveVisualizationConfig)
    equations = semi.equations
    field = something(live_visualization.field, default_live_field(equations))
    label = field_label(field)

    if live_visualization.geometry_mode === :cartesian
        grids = compute_analysis_grids(solution_vector, semi; nvisnodes=live_visualization.nvisnodes, log=false)
        values = select_live_field(grids, field)
        return LiveFieldSnapshot(
            field,
            :cartesian,
            label,
            collect(grids.x),
            collect(grids.y),
            values,
            grids.mask,
            nothing,
        )
    end

    mesh_data = compute_mesh_native_analysis(solution_vector, semi; refine=live_visualization.refine)
    values = select_live_field(mesh_data, field)
    return LiveFieldSnapshot(
        field,
        :mesh_native,
        label,
        mesh_data.x,
        mesh_data.y,
        values,
        nothing,
        mesh_data.triangles,
    )
end

function build_progress_snapshot(
    accepted_steps::Integer,
    current_time::Real,
    final_time::Real,
    residual::Real,
    residual_tol::Real,
    initial_residual::Real;
    stop_reason::Symbol=:running,
)
    residual_fraction = residual_progress_fraction(initial_residual, residual, residual_tol)
    time_fraction = clamp(Float64(current_time) / Float64(final_time), 0.0, 1.0)
    leading = if residual <= residual_tol
        :steady_state
    elseif current_time >= final_time
        :final_time
    else
        residual_fraction >= time_fraction ? :steady_state : :final_time
    end

    return LiveProgressSnapshot(
        Int(accepted_steps),
        Float64(current_time),
        Float64(final_time),
        Float64(residual),
        Float64(residual_tol),
        residual_fraction,
        time_fraction,
        leading,
        stop_reason,
    )
end

function compute_residual(integrator, semi)
    du_ode = Trixi.get_du(integrator)
    integrator.f(du_ode, integrator.u, integrator.p, integrator.t)
    du = Trixi.wrap_array(du_ode, semi)
    return Trixi.residual_steady_state(du, semi.equations)
end

function compute_residual(prob, state, t, semi)
    du_ode = similar(state)
    prob.f(du_ode, state, prob.p, t)
    du = Trixi.wrap_array(du_ode, semi)
    return Trixi.residual_steady_state(du, semi.equations)
end

function create_monitor_state(ode, semi, live_visualization::Union{Nothing, LiveVisualizationConfig})
    initial_residual = compute_residual(ode, ode.u0, 0.0, semi)
    return SolveMonitorState(live_visualization, nothing, initial_residual, -Inf, 0, false)
end

function initialize_live_dashboard!(state::SolveMonitorState, solution_vector, semi, config::SolverConfig, name::AbstractString)
    isnothing(state.live_visualization) && return nothing
    progress = build_progress_snapshot(
        0,
        0.0,
        config.tspan_end,
        state.initial_residual,
        config.residual_tol,
        state.initial_residual;
        stop_reason=:running,
    )
    snapshot = LiveVisualizationSnapshot(progress, build_live_field_snapshot(solution_vector, semi, state.live_visualization))
    state.dashboard = create_live_dashboard(state.live_visualization, snapshot; name=name)
    state.last_update_time = time()
    return nothing
end

function should_update_live_visualization(state::SolveMonitorState, accepted_steps::Integer)
    isnothing(state.live_visualization) && return false
    accepted_steps % state.live_visualization.accepted_step_interval == 0 || return false
    return (time() - state.last_update_time) >= state.live_visualization.min_update_seconds
end

function solve_monitor_callback(config::SolverConfig, semi, state::SolveMonitorState)
    return SciMLBase.DiscreteCallback(
        (u, t, integrator) -> true,
        integrator -> begin
            accepted_steps = integrator.stats.naccept
            state.last_accepted_steps = accepted_steps

            if !isnothing(state.dashboard) && !live_dashboard_is_open(state.dashboard)
                state.window_closed = true
                SciMLBase.terminate!(integrator)
                return nothing
            end

            log_due = accepted_steps > 0 && accepted_steps % config.log_every == 0
            update_due = should_update_live_visualization(state, accepted_steps)
            if !(log_due || update_due)
                return nothing
            end

            residual = compute_residual(integrator, semi)
            du_ode = Trixi.get_du(integrator)
            integrator.f(du_ode, integrator.u, integrator.p, integrator.t)
            du = Trixi.wrap_array(du_ode, semi)
            u = Trixi.wrap_array(integrator.u, semi)
            u_norm = Trixi.residual_steady_state(u, semi.equations)
            rel_residual = u_norm > 0 ? residual / u_norm : Inf
            progress = build_progress_snapshot(
                accepted_steps,
                integrator.t,
                config.tspan_end,
                residual,
                config.residual_tol,
                state.initial_residual;
                stop_reason=:running,
            )

            if log_due
                @info "Progress" iter=accepted_steps t=round(integrator.t, digits=4) dt=round(integrator.dt, digits=6) residual=round(residual, sigdigits=3) rel_residual=round(rel_residual, sigdigits=3) tolerance=config.residual_tol residual_progress=round(progress.residual_progress, digits=3) time_progress=round(progress.time_progress, digits=3) leading_stop_condition=progress.leading_stop_condition
                if has_current_contact_boundary(semi)
                    targets = current_contact_target_summary(semi)
                    measured = measure_midline_jy(integrator.u, semi)
                    @info "Current-drive progress" iter=accepted_steps target_outward_flux_sum=targets.sum target_outward_flux_abs_sum=targets.abs_sum target_outward_flux_max_abs=targets.max_abs target_contact_count=targets.count measured_midline_avg_jy=measured.average measured_midline_integrated_jy=measured.integrated measured_line_y=measured.y measured_points=measured.points
                end
                flush(stdout)
                flush(stderr)
            end

            if update_due
                snapshot = LiveVisualizationSnapshot(
                    progress,
                    build_live_field_snapshot(integrator.u, semi, state.live_visualization),
                )
                update_live_dashboard!(state.dashboard, snapshot)
                state.last_update_time = time()
            end
            return nothing
        end;
        save_positions=(false, false),
    )
end

function finalize_live_dashboard!(state::SolveMonitorState, solution_vector, semi, status, config::SolverConfig)
    isnothing(state.live_visualization) && return nothing
    isnothing(state.dashboard) && return nothing
    progress = build_progress_snapshot(
        state.last_accepted_steps,
        status.final_time,
        config.tspan_end,
        status.final_residual,
        config.residual_tol,
        state.initial_residual;
        stop_reason=status.stop_reason,
    )
    snapshot = LiveVisualizationSnapshot(progress, build_live_field_snapshot(solution_vector, semi, state.live_visualization))
    live_dashboard_is_open(state.dashboard) && finalize_live_dashboard!(state.dashboard, snapshot)
    return nothing
end
