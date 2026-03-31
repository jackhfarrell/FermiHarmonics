using ElectronKinetics, Trixi

function main()
    project_root = normpath(joinpath(@__DIR__, "..", "..", ".."))
    mesh_path = joinpath(project_root, "projects", "nonlinearities", "mesh", "tesla_valve.inp")
    output_dir = joinpath(project_root, "projects", "nonlinear", "data")
    mkpath(output_dir)

    reference = ElectronKinetics.blg_reference_setup()
    mu0 = reference.mu0
    mass = reference.mass
    gamma_mr = 0.01
    gamma_ee = 0.01
    bias = 0.2
    chi = 0.0
    max_harmonic = 50
    p_scatter = reference.p_scatter

    boundary_conditions = Dict(
        :walls => MaxwellWallBC(p_scatter),
        :inlet => OhmicContactBC(bias / 2),
        :outlet => OhmicContactBC(-bias / 2),
    )

    params = SolveParams(;
        polydeg=3,
        tspan_end=2.0,
        residual_tol=1e-3,
        cfl=0.8,
        log_every=200,
        min_harmonic=8,
        max_harmonic_auto=50,
    )

    run_name = "tesla_valve_live_mesh_native_gamma_ee001_M50_bias02_chi000_poly3_t2"
    live_h5_path = joinpath(project_root, "live_viz_$(run_name)_mesh_native.h5")
    watcher_script = joinpath(project_root, "projects", "nonlinear", "scripts", "watch_mesh_native_live.py")

    @info "Running Tesla valve nonlinear live mesh-native case" mu0 mass gamma_mr gamma_ee bias chi max_harmonic polydeg=params.polydeg tspan_end=params.tspan_end residual_tol=params.residual_tol

    if isfile(watcher_script)
        run(Cmd(["python3", watcher_script, live_h5_path]); wait=false)
        @info "Started live viewer window" watcher_script live_h5_path
    end

    sol, semi = ElectronKinetics.solve(
        mesh_path,
        boundary_conditions,
        params,
        gamma_mr,
        gamma_ee;
        transport=:parabolic_nonlinear,
        max_harmonic=max_harmonic,
        mu0=mu0,
        mass=mass,
        chi=chi,
        visualize=true,
        visualize_every=2000,
        visualization_mode=:mesh_native,
        name=run_name,
    )

    cartesian_path = joinpath(output_dir, "$(run_name).h5")
    mesh_native_path = joinpath(output_dir, "$(run_name)_mesh_native.h5")
    ElectronKinetics.save_for_analysis(sol, semi, cartesian_path)
    ElectronKinetics.save_mesh_native_analysis(sol, semi, mesh_native_path; refine=6)
    @info "Saved Tesla valve nonlinear analysis output" cartesian_path mesh_native_path final_time=sol.t[end]

    return nothing
end

main()
