# single run script for the simple_geometries intersection device, as simulated by FermiHarmonics.jl. optional live
# visualization through visualize = true.
# usage: to run using 4 threads for example,
#   JULIA_NUM_THREADS=4 julia --project=. projects/simple_geometries/scripts/run_intersection.jl

using Plots
using Trixi
using StaticArrays
using FermiHarmonics
using DrWatson


# ======================================================================================================================
# Live Visualization Override (a1 instead of j_norm)
# ======================================================================================================================

@inline function a1_only_variables(u, equations::FermiHarmonics.FermiHarmonics2D{NVARS}) where {NVARS}
    a1 = length(u) >= 2 ? u[2] : 0.0
    return SVector{NVARS, Float64}(ntuple(i -> i == 1 ? a1 : 0.0, NVARS))
end

function Trixi.varnames(::typeof(a1_only_variables), equations::FermiHarmonics.FermiHarmonics2D{NVARS}) where {NVARS}
    return ntuple(i -> i == 1 ? "a1" : "_viz_pad_$(i)", NVARS)
end

function FermiHarmonics.visualization_callback(params::FermiHarmonics.SolveParams, semi, name::AbstractString)
    return Trixi.VisualizationCallback(
        semi;
        interval=params.log_every,
        solution_variables=a1_only_variables,
        variable_names=["a1"],
        filename="live_viz_$(name)_a1",
        overwrite=true,
        seriescolor=:magma,
    )
end


# ======================================================================================================================
# Configuration and Setup
# ======================================================================================================================

name = "simple_geometries_intersection"
project_root = normpath(joinpath(@__DIR__, ".."))
mesh_path = joinpath(project_root, "meshes", "intersection", "intersection.inp")

visualize = true
save_analysis = true

bias = 1.0
p_scatter = 1.0
gamma_mr = 0.0
gamma_mc = 1e3

boundary_conditions = Dict(
    :walls => MaxwellWallBC(p_scatter),
    :contact_top => OhmicContactBC(-bias / 2),
    :contact_bottom => OhmicContactBC(bias / 2),
    :contact_right => OhmicContactBC(-bias / 2),
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
    max_harmonic = 5,
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
