using Plots
using ElectronKinetics
using Trixi

function main()
    project_root = normpath(joinpath(@__DIR__, ".."))
    mesh_path = joinpath(project_root, "projects", "square_bells_ucsb", "mesh", "square_bells.inp")
    output_dir = joinpath(@__DIR__, "data_nonlinear")
    mkpath(output_dir)

    reference = blg_reference_setup()
    mu0 = reference.mu0
    mass = reference.mass
    bias = 0.2
    p_scatter = reference.p_scatter

    boundary_conditions = Dict(
        :walls => MaxwellWallBC(p_scatter),
        :contact_top => OhmicContactBC(-bias / 2),
        :contact_bottom => OhmicContactBC(bias / 2),
    )

    config = SolverConfig(;
        polydeg = 1,
        tspan_end = 5.0,
        residual_tol = 1e-4,
        cfl = 0.4,
        log_every = 10,
        min_harmonic = 4,
        max_harmonic_auto = 8,
    )

    gamma_mr = reference.gamma_mr
    gamma_ee = reference.gamma_ee
    run_name = "reference_nonlinear_live_demo"

    @info "Running nonlinear live demo" mu0 mass bias gamma_mr gamma_ee

    model = KineticModel2D(
        Isotropic2DFermiSurface(; vF=1.0, nu=1.0, mass=mass, charge=-1.0),
        HarmonicBasis(:auto),
        IsotropicHarmonicStreaming(),
        QuadraticBGKCollision(gamma_mr, OddQuarticRateProfile(gamma_ee); mu0=mu0, mass=mass),
        reference = reference,
    )

    sol, semi = solve(
        TrixiProblem(; mesh_path=mesh_path, boundary_conditions=boundary_conditions),
        model,
        config;
        visualize = true,
        name = run_name,
    )

    save_path = joinpath(output_dir, "$(run_name).h5")
    save_for_analysis(sol, semi, save_path)
    @info "Saved nonlinear analysis output" path = save_path

    return nothing
end

main()
