using ElectronKinetics
using Trixi
using GLMakie

function main()
    project_root = normpath(joinpath(@__DIR__, ".."))
    mesh_path = joinpath(project_root, "projects", "square_bells_ucsb", "mesh", "square_bells.inp")
    output_dir = joinpath(@__DIR__, "data_square_bells_angle_rate")
    mkpath(output_dir)

    reference = blg_reference_setup()
    mu0 = reference.mu0
    mass = reference.mass
    gamma_mr = 0.01
    gamma_mc = 0.1
    bias = 0.2
    chi = 10.0
    n_angles = 80
    p_scatter = reference.p_scatter

    boundary_conditions = Dict(
        :walls => MaxwellWallBC(p_scatter),
        :contact_top => OhmicContactBC(-bias / 2),
        :contact_bottom => OhmicContactBC(bias / 2),
    )

    params = SolveParams(;
        polydeg=3,
        tspan_end=100.0,
        residual_tol=1e-3,
        cfl=0.8,
        log_every=200,
        min_harmonic=8,
        max_harmonic_auto=50,
    )

    run_name = "square_bells_live_mesh_native_angle_rate_n80_bias02_chi10_gamma_mc01"
    live_h5_path = joinpath(project_root, "live_viz_$(run_name)_mesh_native.h5")

    @info "Running square bells angle-rate live mesh-native case" mu0 mass gamma_mr gamma_mc bias chi n_angles polydeg=params.polydeg tspan_end=params.tspan_end residual_tol=params.residual_tol

    sol, semi = ElectronKinetics.solve(
        mesh_path,
        boundary_conditions,
        params,
        gamma_mr,
        gamma_mc;
        transport=:parabolic_nonlinear,
        collision_model=:angle_rate_bgk,
        n_angles=n_angles,
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
    @info "Saved square bells angle-rate analysis output" cartesian_path mesh_native_path final_time=sol.t[end]

    return nothing
end

main()
