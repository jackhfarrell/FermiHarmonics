# single run script for the Ludwig compass device with live visualization.
# usage:
#   JULIA_NUM_THREADS=4 julia --project=. projects/ludwig/scripts/run_compass_live_viz.jl

using Plots
using Trixi
using FermiHarmonics
using DrWatson

gr()
default(reuse = true)

# ======================================================================================================================
# Configuration and Setup
# ======================================================================================================================

name = get(ENV, "FH_COMPASS_NAME", "ludwig_compass_live_viz")
project_root = normpath(joinpath(@__DIR__, ".."))
mesh_path = joinpath(project_root, "mesh", "compass.inp")

visualize = true
save_analysis = true

bias = parse(Float64, get(ENV, "FH_COMPASS_BIAS", "1.0"))
p_scatter = parse(Float64, get(ENV, "FH_COMPASS_P_SCATTER", "1.0"))
gamma_mr = parse(Float64, get(ENV, "FH_COMPASS_GAMMA_MR", "0.0"))
gamma_ee = parse(Float64, get(ENV, "FH_COMPASS_GAMMA_EE", "100.0"))

boundary_conditions = Dict(
    :walls => MaxwellWallBC(p_scatter),
    :S => MaxwellWallBC(p_scatter),
    :SE => MaxwellWallBC(p_scatter),
    :E => MaxwellWallBC(p_scatter),
    :NE => OhmicContactBC(-bias / 2),
    :N => MaxwellWallBC(p_scatter),
    :NW => OhmicContactBC(bias / 2),
    :W => MaxwellWallBC(p_scatter),
    :SW => MaxwellWallBC(p_scatter),
)

params = SolveParams(;
    min_harmonic = 4,
    max_harmonic_auto = 100,
    polydeg = 3,
    tspan_end = 100.0,
    residual_tol = 1e-5,
    residual_mode = :absolute_all,
    cfl = 0.8,
    log_every = 500,
)

output_dir = joinpath(project_root, "data", "compass_live_viz")
mkpath(output_dir)

function compass_window_plot(plot_data, variable_names;
                             show_mesh = false, plot_arguments = Dict{Symbol, Any}(),
                             time = nothing, timestep = nothing)
    plots = Any[]
    title_suffix = isnothing(timestep) ? "" : " (step $(timestep))"

    for v in variable_names
        push!(plots, Plots.plot(plot_data[v]; plot_arguments..., title = "$(v)$(title_suffix)"))
    end
    if show_mesh
        push!(plots, Plots.plot(Trixi.getmesh(plot_data); plot_arguments..., title = "mesh"))
    end

    cols = ceil(Int, sqrt(length(plots)))
    rows = div(length(plots), cols, RoundUp)
    fig = Plots.plot(plots..., layout = (rows, cols), size = (900, 700))
    display(fig)
    gui(fig)
    return nothing
end

function FermiHarmonics.visualization_callback(params::FermiHarmonics.SolveParams, semi, name::AbstractString)
    return Trixi.VisualizationCallback(
        semi;
        interval = params.log_every,
        solution_variables = FermiHarmonics.current_norm_variables,
        variable_names = ["j_norm"],
        plot_creator = compass_window_plot,
        show_mesh = true,
        seriescolor = :magma,
    )
end

isfile(mesh_path) || error("Compass .inp not found: $mesh_path")

# ======================================================================================================================
# Solve
# ======================================================================================================================

sol, semi = FermiHarmonics.solve(
    mesh_path,
    boundary_conditions,
    params,
    gamma_mr,
    gamma_ee;
    max_harmonic = :auto,
    visualize = visualize,
    name = name,
)

# ======================================================================================================================
# Save Analysis Data
# ======================================================================================================================

if save_analysis
    file_params = (
        device = "compass",
        bias = bias,
        p_scatter = p_scatter,
        gamma_mr = gamma_mr,
        gamma_ee = gamma_ee,
    )
    small_filename = joinpath(output_dir, "observables_" * DrWatson.savename(file_params, "h5"))
    FermiHarmonics.save_for_analysis(sol, semi, small_filename)
    @info "Saved analysis: $(basename(small_filename))"
end
