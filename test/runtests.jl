using Test
using ElectronKinetics
using StaticArrays

@testset "Core Package Loads Without Trixi" begin
    @test Base.get_extension(ElectronKinetics, :ElectronKineticsTrixiExt) === nothing

    surface = Isotropic2DFermiSurface(; vF=1.2, nu=1.5, mass=2.0, charge=-1.0)
    basis = HarmonicBasis(:auto)
    streaming = IsotropicHarmonicStreaming()
    collision = LinearBGKCollision(0.1, TwoRateProfile(0.4))
    model = KineticModel2D(surface, basis, streaming, collision)

    @test ElectronKinetics.transport_symbol(model) === :linear
    @test ElectronKinetics.collision_symbol(model.collision) === :linear_mrt
    @test harmonic_state_nvars(4) == 9
    @test estimate_max_harmonic(0.0, 0.0; min_harmonic=4, max_harmonic=100) == 100
    @test mode_rate(TwoRateProfile(0.7), 2) ≈ 0.7
    @test mode_rate(OddQuarticRateProfile(1.2, 0.05), 3) ≈ min(1.2, 0.05 * 3^4)
    @test mode_rate(ConstantModeRateProfile(0.9), 8) ≈ 0.9
    @test mode_rate(CustomModeRateProfile(m -> 0.2m), 5) ≈ 1.0
    @test model.discretization.max_harmonic === :auto
    @test residual_progress_fraction(10.0, 1.0e-3, 1.0e-2) == 1.0
    @test 0.0 <= residual_progress_fraction(10.0, 1.0, 1.0e-2) < 1.0

    live_config = LiveVisualizationConfig(; accepted_step_interval=25, min_update_seconds=0.0, show_window=false)
    @test live_config.geometry_mode === :mesh_native
    @test live_config.accepted_step_interval == 25
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
    @test Base.get_extension(ElectronKinetics, :ElectronKineticsTrixiExt) !== nothing
end

using GLMakie

@testset "Makie Extension Loads" begin
    @test Base.get_extension(ElectronKinetics, :ElectronKineticsMakieExt) !== nothing
end

const TrixiExt = Base.get_extension(ElectronKinetics, :ElectronKineticsTrixiExt)
const MakieExt = Base.get_extension(ElectronKinetics, :ElectronKineticsMakieExt)

const TESLA_MESH = normpath(joinpath(@__DIR__, "..", "projects", "nonlinearities", "mesh", "tesla_valve.inp"))
const TESLA_BCS = Dict(
    :walls => MaxwellWallBC(1.0),
    :inlet => OhmicContactBC(0.05),
    :outlet => OhmicContactBC(-0.05),
)

@testset "Progress And Cadence Helpers" begin
    steady_progress = TrixiExt.build_progress_snapshot(20, 0.1, 1.0, 1.0e-4, 1.0e-3, 1.0)
    time_progress = TrixiExt.build_progress_snapshot(20, 0.95, 1.0, 1.0, 1.0e-3, 10.0)

    @test steady_progress.leading_stop_condition === :steady_state
    @test time_progress.leading_stop_condition === :final_time

    state = TrixiExt.SolveMonitorState(
        LiveVisualizationConfig(; accepted_step_interval=5, min_update_seconds=0.2, show_window=false),
        nothing,
        1.0,
        time() - 1.0,
        0,
        false,
    )
    @test TrixiExt.should_update_live_visualization(state, 5)
    state.last_update_time = time()
    @test !TrixiExt.should_update_live_visualization(state, 5)
    @test !TrixiExt.should_update_live_visualization(state, 4)
end

@testset "Makie Dashboard Updates In Place" begin
    config = LiveVisualizationConfig(; geometry_mode=:cartesian, show_window=false)
    snapshot0 = LiveVisualizationSnapshot(
        LiveProgressSnapshot(0, 0.0, 1.0, 1.0, 1.0e-3, 0.0, 0.0, :steady_state, :running),
        LiveFieldSnapshot(
            :a0,
            :cartesian,
            "a0",
            collect(range(-1.0, 1.0; length=8)),
            collect(range(-1.0, 1.0; length=8)),
            rand(8, 8),
            trues(8, 8),
            nothing,
        ),
    )
    dashboard = ElectronKinetics.create_live_dashboard(config, snapshot0; name="dashboard_test")
    original_plot = dashboard.field_plot

    snapshot1 = LiveVisualizationSnapshot(
        LiveProgressSnapshot(10, 0.5, 1.0, 1.0e-2, 1.0e-3, 0.5, 0.5, :steady_state, :running),
        LiveFieldSnapshot(
            :a0,
            :cartesian,
            "a0",
            snapshot0.field.x,
            snapshot0.field.y,
            fill(2.0, 8, 8),
            trues(8, 8),
            nothing,
        ),
    )
    ElectronKinetics.update_live_dashboard!(dashboard, snapshot1)

    @test dashboard.field_plot === original_plot
    @test all(dashboard.field_values[] .== 2.0)
end

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
    sol, semi = solve(
        TrixiProblem(; mesh_path=TESLA_MESH, boundary_conditions=TESLA_BCS),
        model,
        config;
        name="test_linear",
        live_visualization=LiveVisualizationConfig(; geometry_mode=:cartesian, accepted_step_interval=5, min_update_seconds=0.0, show_window=false, nvisnodes=24),
    )
    status = solve_status(sol, semi, config)
    overridden_status = solve_status(sol, semi, config; stop_reason_override=:window_closed)

    @test length(sol.u[end]) == length(Trixi.wrap_array(sol.u[end], semi))
    @test status.successful
    @test overridden_status.stop_reason === :window_closed
end

@testset "Legacy SolveParams Wrapper" begin
    params = SolveParams(;
        polydeg=1,
        tspan_end=0.005,
        residual_tol=1e-3,
        cfl=0.2,
        log_every=10_000,
        min_harmonic=2,
        max_harmonic_auto=4,
    )
    sol, semi = ElectronKinetics.solve(
        TESLA_MESH,
        TESLA_BCS,
        params,
        0.0,
        0.5;
        max_harmonic=2,
        visualize=false,
        name="test_legacy_wrapper",
    )
    status = ElectronKinetics.solve_status(sol, semi, params)

    @test ElectronKinetics.SolveParams === ElectronKinetics.SolverConfig
    @test status.successful
    @test length(sol.u[end]) == length(Trixi.wrap_array(sol.u[end], semi))
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
    sol, semi = solve(
        TrixiProblem(; mesh_path=TESLA_MESH, boundary_conditions=TESLA_BCS),
        model,
        config;
        name="test_quadratic",
        live_visualization=LiveVisualizationConfig(; geometry_mode=:mesh_native, accepted_step_interval=5, min_update_seconds=0.0, show_window=false, refine=2),
    )
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
    snapshot = TrixiExt.build_live_field_snapshot(sol.u[end], semi, LiveVisualizationConfig(; geometry_mode=:mesh_native, show_window=false, refine=2))

    @test length(sol.u[end]) == length(Trixi.wrap_array(sol.u[end], semi))
    @test snapshot.field === :current_magnitude
end

@testset "Multiband Snapshot Extraction" begin
    config = SolverConfig(;
        polydeg=1,
        tspan_end=0.005,
        residual_tol=1e-3,
        cfl=0.2,
        log_every=10_000,
        min_harmonic=2,
        max_harmonic_auto=4,
    )
    bands = [
        BandSpec(name=:light, vF=1.0, nu=1.5, mass=1.0, charge=-1.0, gamma_mr=0.0, gamma_mc=0.2),
        BandSpec(name=:heavy, vF=0.8, nu=2.0, mass=3.0, charge=1.0, gamma_mr=0.1, gamma_mc=0.4),
    ]
    model = KineticModel2D(
        Isotropic2DFermiSurface(),
        HarmonicBasis(2),
        IsotropicHarmonicStreaming(),
        LinearBGKCollision(0.0, TwoRateProfile(0.0));
        bands=bands,
        gamma_drag=0.7,
    )
    equations, _ = TrixiExt.build_equations(model, config)
    boundary_symbols = sort(collect(keys(TESLA_BCS)))
    boundary_conditions = (; TESLA_BCS...)
    solver = Trixi.DGSEM(polydeg=config.polydeg, surface_flux=Trixi.flux_lax_friedrichs)
    mesh = Trixi.P4estMesh{2}(TESLA_MESH; boundary_symbols=boundary_symbols)
    semi = Trixi.SemidiscretizationHyperbolic(
        mesh, equations, (x, t, eq) -> zeros(StaticArrays.SVector{Trixi.nvariables(equations), Float64}), solver;
        boundary_conditions=boundary_conditions,
        source_terms=TrixiExt.source_terms,
    )
    ode = Trixi.semidiscretize(semi, (0.0, config.tspan_end))
    snapshot = TrixiExt.build_live_field_snapshot(ode.u0, semi, LiveVisualizationConfig(; geometry_mode=:mesh_native, show_window=false, refine=2))

    @test snapshot.field === :n
end
