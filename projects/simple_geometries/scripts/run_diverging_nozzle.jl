# simple single-run script for the simple_geometries diverging nozzle.
# usage:
#   JULIA_NUM_THREADS=9 julia --project=. projects/simple_geometries/scripts/run_diverging_nozzle.jl

using Plots
using FermiHarmonics

name = "simple_geometries_diverging_nozzle"
project_root = normpath(joinpath(@__DIR__, ".."))
mesh_path = joinpath(project_root, "meshes", "diverging", "diverging.inp")

visualize = true

bias = 1.0
p_scatter = 1.0
gamma_mr = 0.01
gamma_mc = 10.0

boundary_conditions = Dict(
    :walls => MaxwellWallBC(p_scatter),
    :contact_top => OhmicContactBC(-bias / 2),
    :contact_bottom => OhmicContactBC(bias / 2),
)

params = SolveParams(;
    polydeg = 4,
    tspan_end = 100.0,
    residual_tol = 1e-9,
    cfl = 0.8,
    log_every = 500,
)

max_harmonic = 150

isfile(mesh_path) || error("Diverging nozzle .inp not found: $mesh_path")

FermiHarmonics.solve(
    mesh_path,
    boundary_conditions,
    params,
    gamma_mr,
    gamma_mc;
    max_harmonic=:auto,
    visualize=visualize,
    name=name,
)
