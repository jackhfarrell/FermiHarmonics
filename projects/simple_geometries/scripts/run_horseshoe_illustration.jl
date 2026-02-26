# single run script for the simple_geometries horseshoe device with live visualization.
# usage:
#   JULIA_NUM_THREADS=4 julia --project=. projects/simple_geometries/scripts/run_horseshoe_illustration.jl

using Plots
using Trixi
using FermiHarmonics
using DrWatson

# ======================================================================================================================
# Configuration and Setup
# ======================================================================================================================

name = "horseshoe_ullustration"
project_root = normpath(joinpath(@__DIR__, ".."))
mesh_path = joinpath(project_root, "meshes", "snake", "horseshoe.inp")

visualize = true
save_analysis = true

bias = 1.0
p_scatter = 1.0
gamma_mr = 0.1
gamma_mc = 10.0

boundary_conditions = Dict(
    :walls => MaxwellWallBC(p_scatter),
    :contact_in => OhmicContactBC(bias / 2),
    :contact_out => OhmicContactBC(-bias / 2),
)

params = SolveParams(;
    polydeg = 3,
    tspan_end = 100.0,
    residual_tol = 1e-5,
    cfl = 0.8,
    log_every = 500,
)

output_dir = joinpath(project_root, "data", "horseshoe_ullustration")
mkpath(output_dir)

# Write live visualization files into the horseshoe_ullustration folder.
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

isfile(mesh_path) || error("Horseshoe .inp not found: $mesh_path")

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
        device = "horseshoe",
        bias = bias,
        p_scatter = p_scatter,
        gamma_mr = gamma_mr,
        gamma_mc = gamma_mc,
    )
    small_filename = joinpath(output_dir, "observables_" * DrWatson.savename(file_params, "h5"))
    FermiHarmonics.save_for_analysis(sol, semi, small_filename)
    @info "Saved analysis: $(basename(small_filename))"
end
