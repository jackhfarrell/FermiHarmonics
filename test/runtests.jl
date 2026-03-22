using Test
using FermiHarmonics
using Trixi

@testset "FermiHarmonics smoke tests" begin
    params = SolveParams()
    @test params.max_harmonic >= params.min_harmonic
    @test params.max_harmonic_auto >= params.min_harmonic

    @test estimate_max_harmonic(0.0, 0.0; min_harmonic=4, max_harmonic=100) == 100
    @test estimate_max_harmonic(0.0, 500.0; min_harmonic=4, max_harmonic=100) == 4

    eq = FermiHarmonics2D(9; gamma_mr=0.1, gamma_mc=1.0, max_harmonic=4)
    eq_mag = FermiHarmonics2D(9; gamma_mr=0.1, gamma_mc=1.0, omega_c=0.25, max_harmonic=4)
    @test typeof(eq) <: Trixi.AbstractEquations{2, 9}
    @test eq.omega_c == 0.0
    @test eq_mag.omega_c ≈ 0.25 atol=1e-12 rtol=1e-12
    @test occursin("ω_c=0.25", sprint(show, eq_mag))

    u = [1.0, 0.5, -0.2, 0.3, -0.4, 0.1, 0.2, -0.05, 0.07]
    source_zero = FermiHarmonics.physical_sources(u, nothing, 0.0, eq)
    source_mag = FermiHarmonics.physical_sources(u, nothing, 0.0, eq_mag)
    @test source_zero[1] == 0.0
    @test source_mag[1] == 0.0
    @test source_mag[2] ≈ source_zero[2] - 0.25 * u[3] atol=1e-12 rtol=1e-12
    @test source_mag[3] ≈ source_zero[3] + 0.25 * u[2] atol=1e-12 rtol=1e-12
    @test source_mag[4] ≈ source_zero[4] - 2 * 0.25 * u[5] atol=1e-12 rtol=1e-12
    @test source_mag[5] ≈ source_zero[5] + 2 * 0.25 * u[4] atol=1e-12 rtol=1e-12

    pure_mode = zeros(9)
    pure_mode[4] = 1.0
    pure_mode_source = FermiHarmonics.physical_sources(pure_mode, nothing, 0.0, FermiHarmonics2D(
        9; gamma_mr=0.0, gamma_mc=0.0, omega_c=0.5, max_harmonic=4,
    ))
    @test pure_mode_source[1] == 0.0
    @test pure_mode_source[4] == 0.0
    @test pure_mode_source[5] ≈ 1.0 atol=1e-12 rtol=1e-12

    mesh_path = normpath(joinpath(@__DIR__, "..", "projects", "square_bells_ucsb", "mesh", "square_bells.inp"))
    boundary_conditions = Dict(
        :walls => MaxwellWallBC(1.0),
        :contact_top => OhmicContactBC(-0.025),
        :contact_bottom => OhmicContactBC(0.025),
    )
    params = SolveParams(
        polydeg=1,
        tspan_end=0.01,
        residual_tol=1e-3,
        cfl=0.2,
        log_every=10,
        min_harmonic=2,
        max_harmonic_auto=4,
    )
    sol_mag, semi_mag = solve(
        mesh_path,
        boundary_conditions,
        params,
        0.0,
        0.5;
        max_harmonic=2,
        omega_c=0.3,
        name="test_magnetic_field",
    )
    @test semi_mag.equations.omega_c ≈ 0.3 atol=1e-12 rtol=1e-12
    @test sol_mag.t[end] ≈ params.tspan_end atol=1e-12 rtol=1e-12
    @test all(isfinite, sol_mag.u[end])

    sol_mag_neg, semi_mag_neg = solve(
        mesh_path,
        boundary_conditions,
        params,
        0.0,
        0.0;
        max_harmonic=2,
        omega_c=-0.3,
        name="test_magnetic_field_negative",
    )
    @test semi_mag_neg.equations.omega_c ≈ -0.3 atol=1e-12 rtol=1e-12
    @test sol_mag_neg.t[end] ≈ params.tspan_end atol=1e-12 rtol=1e-12
    @test all(isfinite, sol_mag_neg.u[end])
end
