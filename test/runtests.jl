using Test
using ElectronKinetics
using StaticArrays
using HDF5

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
        Band(Isotropic2DFermiSurface(; vF=1.0, nu=1.5, mass=1.0, charge=-1.0);
             name=:light, gamma_mr=0.0, gamma_mc=0.2),
        Band(Isotropic2DFermiSurface(; vF=0.8, nu=2.0, mass=3.0, charge=1.0);
             name=:heavy, gamma_mr=0.1, gamma_mc=0.4),
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

@testset "Surface Interface and New Surface Types" begin
    # Isotropic surface interface
    s_iso = Isotropic2DFermiSurface(; vF=1.2, nu=1.5, mass=2.0, charge=-1.0)
    @test surface_vF(s_iso) ≈ 1.2
    @test surface_max_speed(s_iso) ≈ 1.2
    @test surface_vF_angle(s_iso, 0.0) ≈ 1.2
    @test surface_density_of_states(s_iso) ≈ 1.5
    @test surface_mass(s_iso) ≈ 2.0
    @test surface_charge(s_iso) ≈ -1.0

    # Elliptic surface — isotropic limit matches scalar streaming_matrices
    s_ell = EllipticFermiSurface2D(; vF0=1.0, aspect=1.0, nu=1.0, mass=1.0, charge=-1.0)
    s_ref = Isotropic2DFermiSurface(; vF=1.0, nu=1.0, mass=1.0, charge=-1.0)
    for M in 1:3
        Ax_e, Ay_e = streaming_matrices(M, s_ell)
        Ax_r, Ay_r = streaming_matrices(M, s_ref)
        @test Ax_e ≈ Ax_r atol=1e-10
        @test Ay_e ≈ Ay_r atol=1e-10
    end
    # vF(θ) = vF0/hypot(cos θ, sin θ/aspect)
    # aspect>1: max vF = vF0*aspect (at θ=π/2, hypot=1/aspect); aspect<1: max = vF0 (at θ=0)
    s_aniso = EllipticFermiSurface2D(; vF0=1.0, aspect=2.0, nu=1.0, mass=1.0, charge=-1.0)
    @test surface_max_speed(s_aniso) ≈ 2.0        # vF0*aspect = 2.0
    @test surface_vF_angle(s_aniso, 0.0) ≈ 1.0   # θ=0: hypot(1,0)=1
    @test surface_vF_angle(s_aniso, π/2) ≈ 2.0   # θ=π/2: hypot(0,0.5)=0.5 → vF=1/0.5=2

    # GeneralFermiSurface2D with constant vF matches scalar formula
    s_gen = GeneralFermiSurface2D(θ -> 1.3; max_vF=1.3, nu=1.0, mass=1.0, charge=-1.0)
    for M in 1:2
        Ax_g, Ay_g = streaming_matrices(M, s_gen)
        Ax_r, Ay_r = streaming_matrices(M, 1.3)
        @test Ax_g ≈ Ax_r atol=1e-8
        @test Ay_g ≈ Ay_r atol=1e-8
    end
end

@testset "Band with EllipticFermiSurface2D" begin
    s = EllipticFermiSurface2D(; vF0=1.0, aspect=2.0, nu=1.5, mass=1.0, charge=-1.0)
    b = Band(s; name=:electron, gamma_mr=0.1, gamma_mc=0.2)
    @test b.surface === s
    @test surface_vF(b.surface) ≈ 1.0
    @test surface_max_speed(b.surface) ≈ 2.0     # vF0 * aspect = 2.0
    @test surface_density_of_states(b.surface) ≈ 1.5
    @test band_momentum_weight(b) ≈ 1.5 * 1.0 * 1.0
end

@testset "3-band linear model" begin
    bands = [
        Band(Isotropic2DFermiSurface(; vF=1.0, nu=1.0, mass=1.0, charge=-1.0); name=:a, gamma_mr=0.1, gamma_mc=0.2),
        Band(Isotropic2DFermiSurface(; vF=0.8, nu=1.2, mass=1.5, charge=-1.0); name=:b, gamma_mr=0.1, gamma_mc=0.3),
        Band(Isotropic2DFermiSurface(; vF=0.6, nu=0.9, mass=2.0, charge=-1.0); name=:c, gamma_mr=0.05, gamma_mc=0.1),
    ]
    model = KineticModel2D(
        Isotropic2DFermiSurface(),
        HarmonicBasis(2),
        IsotropicHarmonicStreaming(),
        LinearBGKCollision(0.0, TwoRateProfile(0.0));
        bands=bands,
        gamma_drag=0.0,
    )
    @test length(model.bands) == 3
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
const STRAIGHT_CHANNEL_GEO = normpath(joinpath(@__DIR__, "..", "demo", "mesh", "straight_channel.geo"))
const TESLA_BCS = Dict(
    :walls => MaxwellWallBC(1.0),
    :inlet => OhmicContactBC(0.05),
    :outlet => OhmicContactBC(-0.05),
)
const STRAIGHT_CHANNEL_BCS = Dict(
    :walls => MaxwellWallBC(1.0),
    :inlet => OhmicContactBC(0.05),
    :outlet => OhmicContactBC(-0.05),
)

@testset "Gmsh Geo Mesh Generation" begin
    mesh_path_1 = ElectronKinetics.generate_mesh_from_geo(STRAIGHT_CHANNEL_GEO)
    mesh_path_2 = ElectronKinetics.generate_mesh_from_geo(STRAIGHT_CHANNEL_GEO)
    mesh_contents = read(mesh_path_1, String)

    @test isfile(mesh_path_1)
    @test isfile(mesh_path_2)
    @test mesh_path_1 != mesh_path_2
    @test endswith(lowercase(mesh_path_1), ".inp")
    @test occursin("type=CPS4", mesh_contents) || occursin("type=CPE4", mesh_contents)
    @test occursin("*NSET,NSET=inlet", mesh_contents)
    @test occursin("*NSET,NSET=outlet", mesh_contents)
    @test occursin("*NSET,NSET=walls", mesh_contents)

    @test_throws ArgumentError ElectronKinetics.generate_mesh_from_geo(joinpath(tempdir(), "missing.geo"))
    @test_throws ArgumentError ElectronKinetics.resolve_mesh_path(
        STRAIGHT_CHANNEL_GEO,
        Dict(:missing_contact => MaxwellWallBC(1.0)),
    )

    mktempdir() do dir
        tri_geo = joinpath(dir, "triangle_only.geo")
        open(tri_geo, "w") do io
            write(io, """
SetFactory("OpenCASCADE");
Point(1) = {0, 0, 0, 0.2};
Point(2) = {1, 0, 0, 0.2};
Point(3) = {1, 1, 0, 0.2};
Point(4) = {0, 1, 0, 0.2};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Physical Surface("domain") = {1};
Physical Curve("walls") = {1, 2, 3, 4};
""")
        end

        @test_throws ArgumentError ElectronKinetics.generate_mesh_from_geo(
            tri_geo;
            config=MeshBuildConfig(recombine_all=false, algorithm=6, output_mode=:persistent, output_dir=dir),
        )
    end
end

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

@testset "Geometry-Backed TrixiProblem Solve" begin
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
        Isotropic2DFermiSurface(; vF=1.0, nu=1.0, mass=1.0, charge=-1.0),
        HarmonicBasis(2),
        IsotropicHarmonicStreaming(),
        LinearBGKCollision(0.0, ConstantModeRateProfile(0.5)),
    )
    problem = TrixiProblem(; geometry_path=STRAIGHT_CHANNEL_GEO, boundary_conditions=STRAIGHT_CHANNEL_BCS)
    sol, semi = solve(
        problem,
        model,
        config;
        name="test_linear_geometry_problem",
        live_visualization=LiveVisualizationConfig(; geometry_mode=:cartesian, accepted_step_interval=5, min_update_seconds=0.0, show_window=false, nvisnodes=24),
    )
    status = solve_status(sol, semi, config)
    mesh, _, _, _ = Trixi.mesh_equations_solver_cache(semi)
    provenance_attrs = ElectronKinetics.mesh_provenance_attributes(mesh.current_filename)

    @test status.successful
    @test provenance_attrs["mesh_source_geometry"] == basename(STRAIGHT_CHANNEL_GEO)
    @test provenance_attrs["mesh_build_output_mode"] == "temporary"

    mktempdir() do dir
        analysis_path = joinpath(dir, "straight_channel_analysis.h5")
        save_for_analysis(sol, semi, analysis_path; nvisnodes=32)
        h5open(analysis_path, "r") do file
            @test read(HDF5.attributes(file)["mesh_source_geometry"]) == basename(STRAIGHT_CHANNEL_GEO)
            @test read(HDF5.attributes(file)["mesh_build_output_mode"]) == "temporary"
        end
    end
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

@testset "Legacy Geo Solve Wrapper" begin
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
        STRAIGHT_CHANNEL_GEO,
        STRAIGHT_CHANNEL_BCS,
        params,
        0.0,
        0.5;
        max_harmonic=2,
        visualize=false,
        name="test_legacy_geo_wrapper",
    )
    status = ElectronKinetics.solve_status(sol, semi, params)
    mesh, _, _, _ = Trixi.mesh_equations_solver_cache(semi)

    @test status.successful
    @test ElectronKinetics.mesh_provenance_attributes(mesh.current_filename)["mesh_source_geometry"] ==
        basename(STRAIGHT_CHANNEL_GEO)
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
        Band(Isotropic2DFermiSurface(; vF=1.0, nu=1.5, mass=1.0, charge=-1.0); name=:light, gamma_mr=0.0, gamma_mc=0.2),
        Band(Isotropic2DFermiSurface(; vF=0.8, nu=2.0, mass=3.0, charge=1.0); name=:heavy, gamma_mr=0.1, gamma_mc=0.4),
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
