# single run script for the simple_geometries dogleg device with live visualization.
# usage:
#   JULIA_NUM_THREADS=4 julia --project=. projects/simple_geometries/scripts/run_dogleg_live_viz.jl

using Plots
using Trixi
using FermiHarmonics
using DrWatson

# ======================================================================================================================
# Configuration and Setup
# ======================================================================================================================

name = "simple_geometries_dogleg_live_viz"
project_root = normpath(joinpath(@__DIR__, ".."))
mesh_path = joinpath(project_root, "meshes", "dogleg", "dogleg.inp")

visualize = true
save_analysis = true

bias = 1.0
p_scatter = 1.0
gamma_mr = 1e-2
gamma_mc = 1.0

boundary_conditions = Dict(
    :walls => MaxwellWallBC(p_scatter),
    :source => OhmicContactBC(bias / 2),
    :drain => OhmicContactBC(-bias / 2),
)

params = SolveParams(;
    min_harmonic = 4,
    max_harmonic_auto = 100,
    polydeg = 3,
    tspan_end = 100.0,
    residual_tol = 1e-5,
    residual_mode = :absolute_all,
    cfl = 0.5,
    log_every = 500,
)

output_dir = joinpath(project_root, "data", "dogleg_live_viz")
mkpath(output_dir)

# Write live visualization files into the dogleg_live_viz folder.
function FermiHarmonics.visualization_callback(params::FermiHarmonics.SolveParams, semi, name::AbstractString)
    return Trixi.VisualizationCallback(
        semi;
        interval = params.log_every,
        solution_variables = FermiHarmonics.current_norm_variables,
        variable_names = ["j_norm"],
        filename = joinpath(output_dir, "live_viz_$(name)"),
        overwrite = true,
        seriescolor = :magma,
    )
end

isfile(mesh_path) || error("Dogleg .inp not found: $mesh_path")

# ======================================================================================================================
# Solve
# ======================================================================================================================

sol, semi = FermiHarmonics.solve(
    mesh_path,
    boundary_conditions,
    params,
    gamma_mr,
    gamma_mc;
    max_harmonic = :auto,
    visualize = visualize,
    name = name,
)

# ======================================================================================================================
# Save Analysis Data
# ======================================================================================================================

if save_analysis
    file_params = (
        device = "dogleg",
        bias = bias,
        p_scatter = p_scatter,
        gamma_mr = gamma_mr,
        gamma_mc = gamma_mc,
    )
    small_filename = joinpath(output_dir, "observables_" * DrWatson.savename(file_params, "h5"))
    FermiHarmonics.save_for_analysis(sol, semi, small_filename)
    @info "Saved analysis: $(basename(small_filename))"
end
