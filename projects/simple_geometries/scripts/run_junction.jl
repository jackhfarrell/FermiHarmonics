# single run script for the simple_geometries junction device, as simulated by FermiHarmonics.jl. optional live
# visualization through visualize = true.
# usage: to run using 4 threads for example,
#   JULIA_NUM_THREADS=4 julia --project=. projects/simple_geometries/scripts/run_junction.jl

using Plots
using FermiHarmonics
using DrWatson


# ======================================================================================================================
# Configuration and Setup
# ======================================================================================================================

name = "simple_geometries_junction"
project_root = normpath(joinpath(@__DIR__, ".."))
mesh_path = joinpath(project_root, "meshes", "junction", "junction_coarse.inp")

visualize = true
save_analysis = true

bias = 1.0
p_scatter = 1.0
gamma_mr = 1e-3
gamma_mc = 1e-3

boundary_conditions = Dict(
    :walls => MaxwellWallBC(p_scatter),
    :middle => OhmicContactBC(-bias / 2),
    :bottom => OhmicContactBC(bias / 2),
    :left => OhmicContactBC(-bias / 2),
    :right => OhmicContactBC(-bias / 2),
)

params = SolveParams(;
    polydeg = 3,
    tspan_end = 100.0,
    residual_tol = 1e-5,
    cfl = 0.8,
    log_every = 1000,
)

output_dir = joinpath(project_root, "data", "sims")
mkpath(output_dir)

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
    file_params = (bias = bias, p_scatter = p_scatter, gamma_mr = gamma_mr, gamma_mc = gamma_mc)
    small_filename = joinpath(output_dir, "observables_" * DrWatson.savename(file_params, "h5"))
    FermiHarmonics.save_for_analysis(sol, semi, small_filename)
    @info "Saved analysis: $(basename(small_filename))"
end
