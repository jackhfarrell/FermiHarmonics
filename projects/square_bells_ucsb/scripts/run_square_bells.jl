# single run script for the square bells device, as simulated by FermiHarmonics.jl. optional live visualization through
# visualization = true.
# usage: to run using 4 threads for example,
#   JULIA_NUM_THREADS=4 julia --project=. projects/square_bells_ucsb/scripts/run_square_bells.jl

using FermiHarmonics
using DrWatson

# ======================================================================================================================
# Configuration and Setup
# ======================================================================================================================

name = "square_bells_ucsb"
project_root = normpath(joinpath(@__DIR__, ".."))
mesh_path = joinpath(project_root, "mesh", "square_bells.inp")

visualize = true
save_analysis = true

bias = 1.0
p_scatter = 1.0
gamma_mr = 0.0
gamma_mc = 0.0

boundary_conditions = Dict(
    :walls => MaxwellWallBC(p_scatter),
    :contact_top => OhmicContactBC(-bias / 2),
    :contact_bottom => OhmicContactBC(bias / 2),
)

params = SolveParams(;
    max_harmonic = 60,
    polydeg = 3,
    tspan_end = 100.0,
    residual_tol = 1e-5,
    cfl = 0.8,
    log_every = 500,
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
