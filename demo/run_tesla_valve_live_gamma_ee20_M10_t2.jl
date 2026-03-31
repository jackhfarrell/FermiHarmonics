using Plots
using ElectronKinetics, Trixi

function main()
    project_root = normpath(joinpath(@__DIR__, ".."))
    geo_path = joinpath(project_root, "projects", "nonlinearities", "mesh", "tesla_valve.geo")
    mesh_path = ElectronKinetics.generate_mesh_from_geo(
        geo_path;
        config=MeshBuildConfig(output_mode=:persistent),
    )
    output_dir = joinpath(@__DIR__, "data_nonlinear")
    mkpath(output_dir)

    reference = ElectronKinetics.blg_reference_setup()
    mu0 = reference.mu0
    mass = reference.mass
    gamma_mr = reference.gamma_mr
    gamma_ee = 20.0
    bias = 1.0
    chi = 1.0
    max_harmonic = 10
    p_scatter = reference.p_scatter

    boundary_conditions = Dict(
        :walls => MaxwellWallBC(p_scatter),
        :inlet => OhmicContactBC(bias / 2),
        :outlet => OhmicContactBC(-bias / 2),
    )

    params = SolveParams(;
        polydeg = 3,
        tspan_end = 2.0,
        residual_tol = 1e-3,
        cfl = 0.4,
        log_every = 500,
        min_harmonic = 4,
        max_harmonic_auto = 20,
    )

    run_name = "tesla_valve_live_gamma_ee20_M10_bias10_chi10_poly3_t2"

    @info "Running Tesla valve nonlinear live visualization case" mu0 mass gamma_mr gamma_ee bias chi max_harmonic polydeg=params.polydeg tspan_end=params.tspan_end residual_tol=params.residual_tol

    sol, semi = ElectronKinetics.solve(
        mesh_path,
        boundary_conditions,
        params,
        gamma_mr,
        gamma_ee;
        transport = :parabolic_nonlinear,
        max_harmonic = max_harmonic,
        mu0 = mu0,
        mass = mass,
        chi = chi,
        visualize = true,
        name = run_name,
    )

    save_path = joinpath(output_dir, "$(run_name).h5")
    ElectronKinetics.save_for_analysis(sol, semi, save_path)
    @info "Saved Tesla valve nonlinear analysis output" path = save_path final_time = sol.t[end]

    return nothing
end

main()
