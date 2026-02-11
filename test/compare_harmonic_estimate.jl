# compare harmonic estimate against a 100-harmonic reference at fixed scattering.
# this runs both solves serially to convergence, with live visualization enabled, then logs physical-mode differences.
#
# usage:
#   julia --project=. test/compare_harmonic_estimate.jl

using FermiHarmonics
using Trixi
using Plots
using Statistics

# ======================================================================================================================
# Configuration and Setup
# ======================================================================================================================

name = "compare_harmonic_estimate"
project_root = normpath(joinpath(@__DIR__, ".."))
mesh_path = joinpath(project_root, "projects", "square_bells_ucsb", "mesh", "square_bells.inp")

visualize = true

bias = 1.0
p_scatter = 1.0
gamma_mr = 0.0
gamma_mc = 50.0

reference_max_harmonic = 100
estimated_max_harmonic = FermiHarmonics.estimate_max_harmonic(gamma_mr, gamma_mc)

boundary_conditions = Dict(
    :walls => MaxwellWallBC(p_scatter),
    :contact_top => OhmicContactBC(-bias / 2),
    :contact_bottom => OhmicContactBC(bias / 2),
)

params = SolveParams(;
    max_harmonic = reference_max_harmonic,
    min_harmonic = 4,
    max_harmonic_auto = 150,
    polydeg = 3,
    tspan_end = 150.0,
    residual_tol = 1e-5,
    cfl = 0.8,
    log_every = 500,
)

@info "Harmonic estimate comparison setup" gamma_mr gamma_mc reference_max_harmonic estimated_max_harmonic visualize

# ======================================================================================================================
# Solve: 100-harmonic reference
# ======================================================================================================================

@info "Running reference solve" max_harmonic=reference_max_harmonic
sol_ref, semi_ref = FermiHarmonics.solve(
    mesh_path,
    boundary_conditions,
    params,
    gamma_mr,
    gamma_mc;
    max_harmonic = reference_max_harmonic,
    visualize = visualize,
    name = "$(name)_reference_m$(reference_max_harmonic)",
)

# ======================================================================================================================
# Solve: estimated harmonic count
# ======================================================================================================================

@info "Running estimated solve" max_harmonic=estimated_max_harmonic
sol_est, semi_est = FermiHarmonics.solve(
    mesh_path,
    boundary_conditions,
    params,
    gamma_mr,
    gamma_mc;
    max_harmonic = estimated_max_harmonic,
    visualize = visualize,
    name = "$(name)_estimated_m$(estimated_max_harmonic)",
)

# ======================================================================================================================
# Compare physical modes (a0, a1, b1)
# ======================================================================================================================

u_ref = Trixi.wrap_array(sol_ref.u[end], semi_ref)
u_est = Trixi.wrap_array(sol_est.u[end], semi_est)

u_ref_phys = Array(u_ref[1:3, .., :])
u_est_phys = Array(u_est[1:3, .., :])

size(u_ref_phys) == size(u_est_phys) ||
    throw(ArgumentError("Grid mismatch between runs: ref=$(size(u_ref_phys)), est=$(size(u_est_phys))"))

@info "Converged run summary" reference_max_harmonic estimated_max_harmonic reference_iters=sol_ref.destats.naccept estimated_iters=sol_est.destats.naccept

physical_mode_names = ("a0", "a1", "b1")
aggregate_linf = 0.0
aggregate_rel_linf = 0.0

for (mode_index, mode_name) in enumerate(physical_mode_names)
    ref_mode = u_ref_phys[mode_index, .., :]
    est_mode = u_est_phys[mode_index, .., :]
    delta = est_mode .- ref_mode

    linf = maximum(abs.(delta))
    mean_abs = mean(abs.(delta))
    rms = sqrt(mean(delta .^ 2))
    ref_linf = max(maximum(abs.(ref_mode)), eps(Float64))
    rel_linf = linf / ref_linf

    aggregate_linf = max(aggregate_linf, linf)
    aggregate_rel_linf = max(aggregate_rel_linf, rel_linf)

    @info "Physical mode difference stats" mode=mode_name linf mean_abs rms rel_linf
end

@info "Physical-mode aggregate stats" linf=aggregate_linf rel_linf=aggregate_rel_linf
