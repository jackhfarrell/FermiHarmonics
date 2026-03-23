module ElectronKineticsMakieExt

using ElectronKinetics
using GLMakie

import ElectronKinetics: LiveVisualizationConfig,
                         LiveVisualizationSnapshot,
                         create_live_dashboard,
                         live_dashboard_is_open,
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
        "status  $(progress.stop_reason)",
        "target  $(progress.leading_stop_condition)",
        "steps   $(progress.accepted_steps)",
        "time    $(round(progress.current_time, digits=4)) / $(round(progress.final_time, digits=4))",
        "resid   $(round(progress.residual, sigdigits=4))",
        "tol     $(round(progress.residual_tol, sigdigits=4))",
    ], "\n")
end

function create_live_dashboard(config::LiveVisualizationConfig, snapshot::LiveVisualizationSnapshot; name::AbstractString="run")
    field_values = Observable(render_values(snapshot))
    colorrange = Observable(finite_colorrange(field_values[]))
    title_text = Observable(field_title(snapshot, name))
    panel_text = Observable(panel_summary(snapshot))

    fig = Figure(; size=(1020, 680), figure_padding=(8, 8, 8, 8))
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
        fig[2, 1:2],
        panel_text;
        tellwidth=false,
        tellheight=false,
        halign=:left,
        valign=:top,
        justification=:left,
    )

    colsize!(fig.layout, 1, Relative(0.94))
    colsize!(fig.layout, 2, Relative(0.06))
    rowsize!(fig.layout, 1, Relative(0.9))
    rowsize!(fig.layout, 2, Relative(0.1))
    colgap!(fig.layout, 8)
    rowgap!(fig.layout, 6)

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

function live_dashboard_is_open(dashboard::MakieLiveDashboard)
    isnothing(dashboard.screen) && return true
    return isopen(dashboard.screen)
end

end
