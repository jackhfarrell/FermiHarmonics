using Test
using FermiHarmonics
using Trixi
using StaticArrays
using LinearAlgebra

@testset "BLG reference convention" begin
    reference = blg_reference_setup()
    @test reference.convention_name == "blg_reference_dimensionless"
    @test reference.channel_length == 1.0
    @test reference.mu0 == 1.0
    @test reference.mass == 2.0
    @test reference.gamma_mr == 0.0
    @test reference.gamma_mc == 0.0
    @test reference.vF ≈ 1.0 atol=1e-12 rtol=1e-12
    @test reference.left_probe_x == -0.3
    @test reference.right_probe_x == 0.3

    geo_path = normpath(joinpath(@__DIR__, "..", "demo", "mesh", "straight_channel.geo"))
    geo_contents = read(geo_path, String)
    @test occursin("length_x = 1.0;", geo_contents)
end

@testset "FermiHarmonics smoke tests" begin
    params = SolveParams()
    @test params.max_harmonic >= params.min_harmonic
    @test params.max_harmonic_auto >= params.min_harmonic

    @test estimate_max_harmonic(0.0, 0.0; min_harmonic=4, max_harmonic=100) == 100
    @test estimate_max_harmonic(0.0, 500.0; min_harmonic=4, max_harmonic=100) == 4

    eq = FermiHarmonics2D(9; gamma_mr=0.1, gamma_mc=1.0, max_harmonic=4)
    @test typeof(eq) <: Trixi.AbstractEquations{2, 9}
end

@testset "Nonlinear transport utilities" begin
    eq = FermiHarmonics2D(
        9;
        gamma_mr=0.1,
        gamma_mc=1.0,
        max_harmonic=4,
        transport=:parabolic_nonlinear,
        mu0=2.0,
        mass=8.0,
        theta_oversample=2,
    )
    @test eq.collision_model === :exact_bgk

    state = [0.3, -0.08, 0.05, 0.03, -0.02, 0.015, -0.01, 0.005, -0.004]
    cache = FermiHarmonics.get_nonlinear_cache(eq)
    FermiHarmonics.harmonic_state_to_samples!(cache.samples, state, eq)
    recovered = zeros(Float64, length(state))
    FermiHarmonics.samples_to_harmonics!(recovered, cache.samples, eq)
    @test recovered ≈ state atol=1e-10 rtol=1e-10

    flux_zero_x = Trixi.flux(zeros(9), 1, eq)
    flux_zero_y = Trixi.flux(zeros(9), 2, eq)
    @test flux_zero_x ≈ zeros(9)
    @test flux_zero_y ≈ zeros(9)

    v0 = sqrt(2.0 * eq.mu0 / eq.mass)
    Ax, Ay = FermiHarmonics.streaming_matrices(4, v0)
    epsilon = 1.0e-7
    nonlinear_x = collect(Trixi.flux(epsilon .* state, 1, eq)) ./ epsilon
    nonlinear_y = collect(Trixi.flux(epsilon .* state, 2, eq)) ./ epsilon
    @test nonlinear_x ≈ Ax * state atol=1e-6 rtol=1e-6
    @test nonlinear_y ≈ Ay * state atol=1e-6 rtol=1e-6

    expected_density = eq.mass * (eq.mu0 + 0.5 * state[1]) / (2.0 * pi)
    @test FermiHarmonics.nonlinear_density(state, eq) ≈ expected_density atol=1e-12 rtol=1e-12
    @test collect(FermiHarmonics.nonlinear_current([state[1]; zeros(8)], eq)) ≈ [0.0, 0.0] atol=1e-12 rtol=1e-12

    mu_target = 2.15
    velocity_target = SVector(0.12, -0.05)
    equilibrium_state = zeros(Float64, 9)
    FermiHarmonics.local_equilibrium_state!(equilibrium_state, mu_target, velocity_target, eq)
    recovered_mu, recovered_velocity = FermiHarmonics.recover_mu_u(equilibrium_state, eq)
    expected_density_eq = eq.mass * mu_target / (2.0 * pi)
    @test FermiHarmonics.nonlinear_density(equilibrium_state, eq) ≈ expected_density_eq atol=2e-10 rtol=2e-10
    expected_current_eq = 0.25 * eq.mass * mu_target .* velocity_target
    @test collect(FermiHarmonics.nonlinear_current(equilibrium_state, eq)) ≈ collect(expected_current_eq) atol=5e-9 rtol=5e-9
    @test recovered_mu ≈ mu_target atol=2e-10 rtol=2e-10
    @test recovered_velocity ≈ velocity_target atol=5e-9 rtol=5e-9

    reconstructed_state = zeros(Float64, 9)
    FermiHarmonics.local_equilibrium_state!(reconstructed_state, equilibrium_state, eq)
    @test reconstructed_state ≈ equilibrium_state atol=5e-9 rtol=5e-9

    small_velocity = SVector(0.03, 0.0)
    small_drift_state = zeros(Float64, 9)
    FermiHarmonics.local_equilibrium_state!(small_drift_state, eq.mu0, small_velocity, eq)
    @test small_drift_state[FermiHarmonics.cosine_index(2)] ≈
          small_drift_state[FermiHarmonics.cosine_index(1)]^2 / (4.0 * eq.mu0) atol=2e-4 rtol=2e-3

    @test_throws ArgumentError FermiHarmonics2D(
        3;
        gamma_mr=0.1,
        gamma_mc=0.2,
        max_harmonic=1,
        transport=:parabolic_nonlinear,
        mass=1.0,
    )
    @test_throws DomainError Trixi.flux([-2.5], 1, FermiHarmonics2D(
        1;
        gamma_mr=0.1,
        gamma_mc=0.2,
        transport=:parabolic_nonlinear,
        mu0=1.0,
        mass=1.0,
    ))
end

@testset "Nonlinear BGK source terms" begin
    eq_bgk = FermiHarmonics2D(
        9;
        gamma_mr=0.0,
        gamma_mc=0.8,
        max_harmonic=4,
        transport=:parabolic_nonlinear,
        mu0=2.0,
        mass=8.0,
        theta_oversample=2,
    )
    equilibrium_state = zeros(Float64, 9)
    FermiHarmonics.local_equilibrium_state!(equilibrium_state, 2.1, SVector(0.08, -0.03), eq_bgk)
    source_eq = FermiHarmonics.physical_sources(equilibrium_state, nothing, 0.0, eq_bgk)
    @test source_eq ≈ zeros(9) atol=5e-9 rtol=5e-9

    perturbed = copy(equilibrium_state)
    perturbed[FermiHarmonics.cosine_index(3)] += 0.04
    source_perturbed = FermiHarmonics.physical_sources(perturbed, nothing, 0.0, eq_bgk)
    @test source_perturbed[FermiHarmonics.cosine_index(3)] < 0.0

    eq_two_rate = FermiHarmonics2D(
        9;
        gamma_mr=0.2,
        gamma_mc=0.5,
        max_harmonic=4,
        transport=:parabolic_nonlinear,
        mu0=2.0,
        mass=8.0,
        theta_oversample=2,
    )
    two_rate_state = zeros(Float64, 9)
    FermiHarmonics.local_equilibrium_state!(two_rate_state, 2.05, SVector(0.07, 0.01), eq_two_rate)
    two_rate_source = FermiHarmonics.physical_sources(two_rate_state, nothing, 0.0, eq_two_rate)
    @test norm(two_rate_source[2:3]) > 0.0
    @test norm(two_rate_source[4:end]) > 0.0
end

@testset "Nonlinear boundary conditions" begin
    eq = FermiHarmonics2D(
        9;
        gamma_mr=0.0,
        gamma_mc=0.0,
        max_harmonic=4,
        transport=:parabolic_nonlinear,
        mu0=2.0,
        mass=8.0,
    )
    unit_normal = SVector(1.0, 0.0)
    state = zeros(Float64, 9)
    out = similar(state)
    scratch = similar(state)

    FermiHarmonics.nonlinear_maxwell_wall!(out, state, unit_normal, 1.0, scratch, eq, 1.0e-12)
    @test out ≈ zeros(9)

    FermiHarmonics.nonlinear_ohmic_contact!(out, state, unit_normal, 1.0, 1.0, scratch, eq, 1.0e-12)
    cache = FermiHarmonics.get_nonlinear_cache(eq)
    FermiHarmonics.harmonic_state_to_samples!(cache.samples, out, eq)
    data = FermiHarmonics.nonlinear_data(eq)
    incoming = [real(cache.samples[j]) for j in eachindex(cache.samples)
                if data.cos_theta[j] < -1.0e-12]
    outgoing = [real(cache.samples[j]) for j in eachindex(cache.samples)
                if data.cos_theta[j] > 1.0e-12]
    @test !isempty(incoming)
    @test !isempty(outgoing)
    @test 0.9 < out[1] < 1.0
    @test sum(incoming) / length(incoming) > sum(outgoing) / length(outgoing)
    @test all(isfinite, incoming)
    @test all(isfinite, outgoing)

    out_top = similar(state)
    out_bottom = similar(state)
    FermiHarmonics.nonlinear_ohmic_contact!(out_top, state, SVector(0.0, -1.0), 1.0, 1.0, scratch, eq, 1.0e-12)
    FermiHarmonics.nonlinear_ohmic_contact!(out_bottom, state, SVector(0.0, 1.0), 1.0, 1.0, scratch, eq, 1.0e-12)
    @test out_top[1] ≈ out_bottom[1] atol=1e-12 rtol=1e-12
    @test out_top[2] ≈ out_bottom[2] atol=1e-12 rtol=1e-12
    @test out_top[3] ≈ -out_bottom[3] atol=1e-12 rtol=1e-12
end

@testset "Solve smoke tests" begin
    mesh_path = normpath(joinpath(@__DIR__, "..", "projects", "square_bells_ucsb", "mesh", "square_bells.inp"))
    boundary_conditions = Dict(
        :walls => MaxwellWallBC(1.0),
        :contact_top => OhmicContactBC(-0.1),
        :contact_bottom => OhmicContactBC(0.1),
    )
    params = SolveParams(;
        polydeg=1,
        tspan_end=0.02,
        residual_tol=1e-3,
        cfl=0.4,
        log_every=10_000,
        min_harmonic=2,
        max_harmonic_auto=4,
    )

    sol_linear, semi_linear = solve(
        mesh_path,
        boundary_conditions,
        params,
        0.0,
        0.5;
        max_harmonic=2,
        name="test_linear",
    )
    @test length(sol_linear.u[end]) == length(Trixi.wrap_array(sol_linear.u[end], semi_linear))

    sol_nonlinear, semi_nonlinear = solve(
        mesh_path,
        boundary_conditions,
        params,
        0.0,
        0.5;
        max_harmonic=2,
        transport=:parabolic_nonlinear,
        mu0=2.0,
        mass=8.0,
        theta_oversample=2,
        name="test_nonlinear",
    )
    @test semi_nonlinear.equations.transport === :parabolic_nonlinear
    @test semi_nonlinear.equations.collision_model === :exact_bgk
    @test length(sol_nonlinear.u[end]) == length(Trixi.wrap_array(sol_nonlinear.u[end], semi_nonlinear))

    linear_probe = evaluate_observables(sol_linear, semi_linear, 0.0, 0.0)
    @test linear_probe.in_domain
    @test linear_probe.jx ≈ linear_probe.a1 atol=1e-10 rtol=1e-10
    @test linear_probe.jy ≈ linear_probe.b1 atol=1e-10 rtol=1e-10

    nonlinear_probe = evaluate_observables(sol_nonlinear, semi_nonlinear, 0.0, 0.0)
    @test nonlinear_probe.in_domain
    @test isfinite(nonlinear_probe.n)
    @test isfinite(nonlinear_probe.jx)
    @test isfinite(nonlinear_probe.jy)
    @test Trixi.varnames(FermiHarmonics.analysis_variables, semi_nonlinear.equations) == ("n", "jx", "jy")
end
