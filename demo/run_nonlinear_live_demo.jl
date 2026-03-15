using Plots
using FermiHarmonics

function main()
    project_root = normpath(joinpath(@__DIR__, ".."))
    mesh_path = joinpath(project_root, "projects", "square_bells_ucsb", "mesh", "square_bells.inp")
    output_dir = joinpath(@__DIR__, "data_nonlinear")
    mkpath(output_dir)

    reference = FermiHarmonics.blg_reference_setup()
    mu0 = reference.mu0
    mass = reference.mass
    bias = 0.2
    p_scatter = reference.p_scatter

    boundary_conditions = Dict(
        :walls => MaxwellWallBC(p_scatter),
        :contact_top => OhmicContactBC(-bias / 2),
        :contact_bottom => OhmicContactBC(bias / 2),
    )

    params = SolveParams(;
        polydeg = 1,
        tspan_end = 5.0,
        residual_tol = 1e-4,
        cfl = 0.4,
        log_every = 10,
        min_harmonic = 4,
        max_harmonic_auto = 8,
    )

    gamma_mr = reference.gamma_mr
    gamma_mc = reference.gamma_mc
    run_name = "reference_nonlinear_live_demo"

    @info "Running nonlinear live demo" mu0 mass bias gamma_mr gamma_mc

    sol, semi = FermiHarmonics.solve(
        mesh_path,
        boundary_conditions,
        params,
        gamma_mr,
        gamma_mc;
        max_harmonic = 4,
        transport = :parabolic_nonlinear,
        mu0 = mu0,
        mass = mass,
        theta_oversample = 2,
        visualize = true,
        name = run_name,
    )

    save_path = joinpath(output_dir, "$(run_name).h5")
    FermiHarmonics.save_for_analysis(sol, semi, save_path)
    @info "Saved nonlinear analysis output" path = save_path

    return nothing
end

main()
