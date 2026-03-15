using Test
using FermiHarmonics
using Trixi
using StaticArrays

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

@testset "Nonlinear boundary conditions" begin
    eq = FermiHarmonics2D(
        5;
        gamma_mr=0.0,
        gamma_mc=0.0,
        max_harmonic=2,
        transport=:parabolic_nonlinear,
        mu0=2.0,
        mass=8.0,
    )
    unit_normal = SVector(1.0, 0.0)
    state = zeros(Float64, 5)
    out = similar(state)
    scratch = similar(state)
    P_in = FermiHarmonics.incoming_projector(eq, unit_normal)

    FermiHarmonics.nonlinear_maxwell_wall!(out, state, unit_normal, P_in, 1.0, scratch, eq)
    @test out ≈ zeros(5)

    FermiHarmonics.nonlinear_ohmic_contact!(out, state, unit_normal, P_in, 1.0, 1.0, scratch)
    cache = FermiHarmonics.get_nonlinear_cache(eq)
    FermiHarmonics.harmonic_state_to_samples!(cache.samples, out, eq)
    data = FermiHarmonics.nonlinear_data(eq)
    incoming = [real(cache.samples[j]) for j in eachindex(cache.samples)
                if data.cos_theta[j] < -1.0e-8]
    outgoing = [real(cache.samples[j]) for j in eachindex(cache.samples)
                if data.cos_theta[j] > 1.0e-8]
    @test !isempty(incoming)
    @test !isempty(outgoing)
    @test out[1] ≈ 0.5 atol=1e-10 rtol=1e-10
    @test sum(incoming) / length(incoming) > sum(outgoing) / length(outgoing)
    @test all(isfinite, incoming)
    @test all(isfinite, outgoing)
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
    @test length(sol_nonlinear.u[end]) == length(Trixi.wrap_array(sol_nonlinear.u[end], semi_nonlinear))
end
