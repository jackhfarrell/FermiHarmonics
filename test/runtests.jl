using Test
using FermiFlows

@testset "Core Package Loads Without Trixi" begin
    @test Base.get_extension(FermiFlows, :FermiFlowsTrixiExt) === nothing

    surface = Isotropic2DFermiSurface(; vF=1.2, nu=1.5, mass=2.0, charge=-1.0)
    basis = HarmonicBasis(:auto)
    streaming = IsotropicHarmonicStreaming()
    collision = LinearBGKCollision(0.1, TwoRateProfile(0.4))
    model = KineticModel2D(surface, basis, streaming, collision)

    @test FermiFlows.transport_symbol(model) === :linear
    @test FermiFlows.collision_symbol(model.collision) === :linear_mrt
    @test harmonic_state_nvars(4) == 9
    @test estimate_max_harmonic(0.0, 0.0; min_harmonic=4, max_harmonic=100) == 100
    @test mode_rate(TwoRateProfile(0.7), 2) ≈ 0.7
    @test mode_rate(OddQuarticRateProfile(1.2, 0.05), 3) ≈ min(1.2, 0.05 * 3^4)
    @test mode_rate(ConstantModeRateProfile(0.9), 8) ≈ 0.9
    @test mode_rate(CustomModeRateProfile(m -> 0.2m), 5) ≈ 1.0
    @test model.discretization.max_harmonic === :auto
end

@testset "Core Multiband Bookkeeping" begin
    bands = [
        BandSpec(
            name=:light,
            vF=1.0,
            nu=1.5,
            mass=1.0,
            charge=-1.0,
            gamma_mr=0.0,
            gamma_mc=0.2,
        ),
        BandSpec(
            name=:heavy,
            vF=0.8,
            nu=2.0,
            mass=3.0,
            charge=1.0,
            gamma_mr=0.1,
            gamma_mc=0.4,
        ),
    ]
    model = KineticModel2D(
        Isotropic2DFermiSurface(),
        HarmonicBasis(2),
        IsotropicHarmonicStreaming(),
        LinearBGKCollision(0.0, TwoRateProfile(0.0));
        bands=bands,
        gamma_drag=0.7,
    )

    @test length(model.bands) == 2
    @test band_momentum_weight(model.bands[1]) ≈ 1.5
    @test band_momentum_weight(model.bands[2]) ≈ 4.8
    @test model.gamma_drag ≈ 0.7
end

using Trixi

@testset "Trixi Extension Loads" begin
    @test Base.get_extension(FermiFlows, :FermiFlowsTrixiExt) !== nothing
end

const TESLA_MESH = normpath(joinpath(@__DIR__, "..", "projects", "nonlinearities", "mesh", "tesla_valve.inp"))
const TESLA_BCS = Dict(
    :walls => MaxwellWallBC(1.0),
    :inlet => OhmicContactBC(0.05),
    :outlet => OhmicContactBC(-0.05),
)

@testset "Linear Harmonic Solve" begin
    config = SolverConfig(;
        polydeg=1,
        tspan_end=0.01,
        residual_tol=1e-3,
        cfl=0.2,
        log_every=10_000,
        min_harmonic=2,
        max_harmonic_auto=4,
    )
    model = KineticModel2D(
        Isotropic2DFermiSurface(; vF=1.0, nu=1.0, mass=1.0, charge=-1.0),
        HarmonicBasis(2),
        IsotropicHarmonicStreaming(),
        LinearBGKCollision(0.0, ConstantModeRateProfile(0.5)),
    )
    sol, semi = solve(TrixiProblem(; mesh_path=TESLA_MESH, boundary_conditions=TESLA_BCS), model, config; name="test_linear")
    status = solve_status(sol, semi, config)

    @test length(sol.u[end]) == length(Trixi.wrap_array(sol.u[end], semi))
    @test status.successful
end

@testset "Nonlinear Harmonic Solve And Export" begin
    config = SolverConfig(;
        polydeg=1,
        tspan_end=0.01,
        residual_tol=1e-3,
        cfl=0.2,
        log_every=10_000,
        min_harmonic=2,
        max_harmonic_auto=4,
    )
    model = KineticModel2D(
        Isotropic2DFermiSurface(; vF=1.0, nu=1.0, mass=2.0, charge=-1.0),
        HarmonicBasis(2),
        IsotropicHarmonicStreaming(),
        QuadraticBGKCollision(0.0, OddQuarticRateProfile(0.5); mu0=1.0, mass=2.0, electrostatic_coupling=0.0),
    )
    sol, semi = solve(TrixiProblem(; mesh_path=TESLA_MESH, boundary_conditions=TESLA_BCS), model, config; name="test_quadratic")
    obs = evaluate_observables(sol, semi, 0.0, 0.0)

    @test hasproperty(obs, :in_domain)

    mktempdir() do dir
        mesh_native_path = joinpath(dir, "tesla_mesh_native.h5")
        save_mesh_native_analysis(sol, semi, mesh_native_path; refine=2, observables=[:n, :jx, :jy])
        @test isfile(mesh_native_path)
    end
end

@testset "Angle BGK Solve" begin
    config = SolverConfig(;
        polydeg=1,
        tspan_end=0.005,
        residual_tol=1e-3,
        cfl=0.2,
        log_every=10_000,
        min_harmonic=2,
        max_harmonic_auto=4,
    )
    model = KineticModel2D(
        Isotropic2DFermiSurface(; vF=1.0, nu=1.0, mass=2.0, charge=-1.0),
        AngleGrid(16),
        IsotropicAngleStreaming(),
        ExactAngleBGKCollision(; gamma_mr=0.0, gamma_mc=0.5, mu0=1.0, mass=2.0, electrostatic_coupling=0.0),
    )
    sol, semi = solve(TrixiProblem(; mesh_path=TESLA_MESH, boundary_conditions=TESLA_BCS), model, config; name="test_exact_angle")

    @test length(sol.u[end]) == length(Trixi.wrap_array(sol.u[end], semi))
end
