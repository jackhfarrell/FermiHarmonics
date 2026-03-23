ENV["GKSwstype"] = "100"
ENV["GKS_NO_GUI"] = "1"

using FermiHarmonics

function main()
    project_root = normpath(joinpath(@__DIR__, ".."))
    mesh_path = joinpath(project_root, "projects", "square_bells_ucsb", "mesh", "square_bells.inp")
    output_dir = joinpath(@__DIR__, "data_nonlinear")
    mkpath(output_dir)

    reference = FermiHarmonics.blg_reference_setup()
    mu0 = reference.mu0
    mass = reference.mass
    gamma_mr = 0.0
    gamma_mc = 200.0
    bias = 0.1
    chi = 10.0
    n_angles = 50
    p_scatter = reference.p_scatter

    boundary_conditions = Dict(
        :walls => MaxwellWallBC(p_scatter),
        :contact_top => OhmicContactBC(-bias / 2),
        :contact_bottom => OhmicContactBC(bias / 2),
    )

    params = SolveParams(;
        polydeg=3,
        tspan_end=2.0,
        residual_tol=1e-3,
        cfl=0.8,
        log_every=1000,
        min_harmonic=4,
        max_harmonic_auto=20,
    )

    run_name = "square_bells_exact_live_streamlines_gamma_mc200_A50_bias10_chi100_poly3_t2"

    @info "Running square-bells exact-angle nonlinear live streamline case" mu0 mass gamma_mr gamma_mc bias chi n_angles polydeg=params.polydeg tspan_end=params.tspan_end residual_tol=params.residual_tol

    sol, semi = FermiHarmonics.solve(
        mesh_path,
        boundary_conditions,
        params,
        gamma_mr,
        gamma_mc;
        transport=:parabolic_nonlinear,
        collision_model=:exact_bgk,
        n_angles=n_angles,
        mu0=mu0,
        mass=mass,
        chi=chi,
        visualize=true,
        name=run_name,
    )

    save_path = joinpath(output_dir, "$(run_name).h5")
    FermiHarmonics.save_for_analysis(sol, semi, save_path)
    @info "Saved exact-angle nonlinear analysis output" path=save_path final_time=sol.t[end]

    return nothing
end

main()
