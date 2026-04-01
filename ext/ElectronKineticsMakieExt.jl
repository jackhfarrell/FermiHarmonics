module ElectronKineticsMakieExt

using ElectronKinetics
using GLMakie

import ElectronKinetics: LiveVisualizationConfig,
                         LiveVisualizationSnapshot,
                         create_mesh_preview,
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

function mesh_native_geometry(snapshot::LiveVisualizationSnapshot)
    points = GLMakie.Point3f[
        GLMakie.Point3f(Float32(x), Float32(y), 0.0f0) for (x, y) in zip(snapshot.field.x, snapshot.field.y)
    ]
    faces = GLMakie.GLTriangleFace[
        GLMakie.GLTriangleFace((
            UInt32(snapshot.field.triangles[index, 1] + 1),
            UInt32(snapshot.field.triangles[index, 2] + 1),
            UInt32(snapshot.field.triangles[index, 3] + 1),
        )) for index in axes(snapshot.field.triangles, 1)
    ]
    return GLMakie.GeometryBasics.Mesh(points, faces)
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
            colormap=config.colormap,
            colorrange=colorrange,
        )
    elseif config.mesh_outline_only
        nothing
    else
        mesh!(
            axis,
            mesh_native_geometry(snapshot);
            color=field_values,
            colormap=config.colormap,
            colorrange=colorrange,
            shading=GLMakie.NoShading,
        )
    end

    if config.mesh_outline && snapshot.field.geometry_mode === :mesh_native
        wireframe!(
            axis,
            mesh_native_geometry(snapshot);
            color=:black,
            linewidth=0.5,
        )
    end

    if field_plot !== nothing
        Colorbar(fig[1, 2], field_plot; label=snapshot.field.label)
    end
    Label(
        fig[2, 1:2],
        panel_text;
        tellwidth=false,
        tellheight=false,
        halign=:left,
        valign=:top,
        justification=:left,
    )

    colsize!(fig.layout, 1, Relative(field_plot === nothing ? 1.0 : 0.94))
    colsize!(fig.layout, 2, Relative(field_plot === nothing ? 0.0 : 0.06))
    rowsize!(fig.layout, 1, Relative(0.9))
    rowsize!(fig.layout, 2, Relative(0.1))
    colgap!(fig.layout, 8)
    rowgap!(fig.layout, 6)

    screen = config.show_window ? display(fig) : nothing
    return MakieLiveDashboard(fig, field_plot, field_values, colorrange, title_text, panel_text, String(name), screen)
end

function create_mesh_preview(
    nodes::AbstractMatrix{<:Real},
    quads::AbstractMatrix{<:Integer};
    scale::Real=1.0,
    name::AbstractString="mesh_preview",
)
    point_count = size(nodes, 1)
    quad_count = size(quads, 1)
    point_count > 0 || throw(ArgumentError("mesh preview requires non-empty node list"))
    quad_count > 0 || throw(ArgumentError("mesh preview requires quad elements"))

    scale_value = Float32(scale)
    points = Vector{GLMakie.Point2f}(undef, point_count)
    @inbounds for i in 1:point_count
        points[i] = GLMakie.Point2f(scale_value * Float32(nodes[i, 1]), scale_value * Float32(nodes[i, 2]))
    end

    segments = GLMakie.Point2f[]
    sizehint!(segments, quad_count * 8)
    @inbounds for row in 1:quad_count
        n1 = quads[row, 1]
        n2 = quads[row, 2]
        n3 = quads[row, 3]
        n4 = quads[row, 4]
        push!(segments, points[n1]); push!(segments, points[n2])
        push!(segments, points[n2]); push!(segments, points[n3])
        push!(segments, points[n3]); push!(segments, points[n4])
        push!(segments, points[n4]); push!(segments, points[n1])
    end

    fig = Figure(; size=(920, 680), figure_padding=(8, 8, 8, 8))
    axis = Axis(fig[1, 1]; title="$(name): mesh outline", xlabel="x", ylabel="y", aspect=DataAspect())
    linesegments!(axis, segments; color=:black, linewidth=0.5)
    screen = display(fig)
    return (figure=fig, axis=axis, screen=screen)
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
