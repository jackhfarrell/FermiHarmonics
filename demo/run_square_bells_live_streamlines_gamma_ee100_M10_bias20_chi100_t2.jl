ENV["GKSwstype"] = "100"
ENV["GKS_NO_GUI"] = "1"

using ElectronKinetics, Trixi

function main()
    project_root = normpath(joinpath(@__DIR__, ".."))
    mesh_path = joinpath(project_root, "projects", "square_bells_ucsb", "mesh", "square_bells.inp")
    output_dir = joinpath(@__DIR__, "data_nonlinear")
    mkpath(output_dir)

    reference = ElectronKinetics.blg_reference_setup()
    mu0 = reference.mu0
    mass = reference.mass
    gamma_mr = reference.gamma_mr
    gamma_ee = 100.0
    bias = 0.2
    chi = 10.0
    max_harmonic = 10
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

    run_name = "square_bells_live_streamlines_gamma_ee100_M10_bias20_chi100_poly3_t2"

    @info "Running square-bells nonlinear live streamline case" mu0 mass gamma_mr gamma_ee bias chi max_harmonic polydeg=params.polydeg tspan_end=params.tspan_end residual_tol=params.residual_tol

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
        name=run_name,
    )

    save_path = joinpath(output_dir, "$(run_name).h5")
    ElectronKinetics.save_for_analysis(sol, semi, save_path)
    @info "Saved nonlinear analysis output" path=save_path final_time=sol.t[end]

    return nothing
end

main()
