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
    bias = 0.5
    chi = 1.0
    max_harmonic = 20
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
        log_every = 1000,
        min_harmonic = 4,
        max_harmonic_auto = 20,
    )

    run_name = "tesla_valve_live_gamma_ee20_M20_bias05_chi10_poly3_t2"

    @info "Running Tesla valve nonlinear batch case" mu0 mass gamma_mr gamma_ee bias chi max_harmonic polydeg=params.polydeg tspan_end=params.tspan_end residual_tol=params.residual_tol

    sol, semi = ElectronKinetics.solve(
        mesh_path,
        boundary_conditions,
        params,
        gamma_mr,
        gamma_ee;
    callbacks=default_callbacks_builder(),
        transport = :parabolic_nonlinear,
        max_harmonic = max_harmonic,
        mu0 = mu0,
        mass = mass,
        chi = chi,
        visualize = false,
        name = run_name,
    )

    cartesian_path = joinpath(output_dir, "$(run_name).h5")
    mesh_native_path = joinpath(output_dir, "$(run_name)_mesh_native.h5")
    ElectronKinetics.save_for_analysis(sol, semi, cartesian_path)
    ElectronKinetics.save_mesh_native_analysis(sol, semi, mesh_native_path; refine = 6)
    @info "Saved Tesla valve nonlinear analysis output" cartesian_path mesh_native_path final_time = sol.t[end]

    return nothing
end

main()
