module FermiFlowsMakieExt

using FermiFlows
using GLMakie

import FermiFlows: LiveVisualizationConfig,
                   LiveVisualizationSnapshot,
                   create_live_dashboard,
                   finalize_live_dashboard!,
                   update_live_dashboard!

mutable struct MakieLiveDashboard
    figure
    field_plot
    field_values
    colorrange
    title_text
    panel_text
    name::String
    screen
end

function ascii_progress_bar(fraction::Real; width::Int=18)
    frac = clamp(Float64(fraction), 0.0, 1.0)
    filled = clamp(round(Int, frac * width), 0, width)
    return "[" * repeat("#", filled) * repeat("-", width - filled) * "]"
end

function finite_colorrange(values)
    finite_values = filter(isfinite, vec(values))
    isempty(finite_values) && return (0.0, 1.0)
    lo = minimum(finite_values)
    hi = maximum(finite_values)
    lo == hi && return (lo - 1.0e-12, hi + 1.0e-12)
    return (lo, hi)
end

function render_values(snapshot)
    if snapshot.field.geometry_mode === :cartesian && snapshot.field.mask !== nothing
        return ifelse.(snapshot.field.mask, snapshot.field.values, NaN)
    end
    return snapshot.field.values
end

function field_title(snapshot::LiveVisualizationSnapshot, name::AbstractString)
    return "$(name): $(snapshot.field.label) at t=$(round(snapshot.progress.current_time, digits=3))"
end

function panel_summary(snapshot::LiveVisualizationSnapshot)
    progress = snapshot.progress
    return join([
        "status: $(progress.stop_reason)",
        "leading stop: $(progress.leading_stop_condition)",
        "steps: $(progress.accepted_steps)",
        "t = $(round(progress.current_time, digits=4)) / $(round(progress.final_time, digits=4))",
        "residual = $(round(progress.residual, sigdigits=4))",
        "tol = $(round(progress.residual_tol, sigdigits=4))",
        "convergence " * ascii_progress_bar(progress.residual_progress) * " $(round(progress.residual_progress * 100, digits=1))%",
        "time        " * ascii_progress_bar(progress.time_progress) * " $(round(progress.time_progress * 100, digits=1))%",
    ], "\n")
end

function create_live_dashboard(config::LiveVisualizationConfig, snapshot::LiveVisualizationSnapshot; name::AbstractString="run")
    field_values = Observable(render_values(snapshot))
    colorrange = Observable(finite_colorrange(field_values[]))
    title_text = Observable(field_title(snapshot, name))
    panel_text = Observable(panel_summary(snapshot))

    fig = Figure(; size=(1100, 720))
    axis = Axis(fig[1, 1]; title=title_text, xlabel="x", ylabel="y", aspect=DataAspect())

    field_plot = if snapshot.field.geometry_mode === :cartesian
        heatmap!(
            axis,
            snapshot.field.x,
            snapshot.field.y,
            field_values;
            colormap=:magma,
            colorrange=colorrange,
        )
    else
        scatter!(
            axis,
            snapshot.field.x,
            snapshot.field.y;
            color=field_values,
            colormap=:magma,
            colorrange=colorrange,
            markersize=6,
        )
    end

    Colorbar(fig[1, 2], field_plot; label=snapshot.field.label)
    Label(
        fig[1, 3],
        panel_text;
        tellwidth=false,
        halign=:left,
        valign=:top,
        justification=:left,
    )

    colsize!(fig.layout, 1, Relative(0.72))
    colsize!(fig.layout, 2, Relative(0.05))
    colsize!(fig.layout, 3, Relative(0.23))

    screen = config.show_window ? display(fig) : nothing
    return MakieLiveDashboard(fig, field_plot, field_values, colorrange, title_text, panel_text, String(name), screen)
end

function update_live_dashboard!(dashboard::MakieLiveDashboard, snapshot::LiveVisualizationSnapshot)
    dashboard.field_values[] = render_values(snapshot)
    dashboard.colorrange[] = finite_colorrange(dashboard.field_values[])
    dashboard.title_text[] = field_title(snapshot, dashboard.name)
    dashboard.panel_text[] = panel_summary(snapshot)
    return dashboard
end

function finalize_live_dashboard!(dashboard::MakieLiveDashboard, snapshot::LiveVisualizationSnapshot)
    update_live_dashboard!(dashboard, snapshot)
    return dashboard
end

end
