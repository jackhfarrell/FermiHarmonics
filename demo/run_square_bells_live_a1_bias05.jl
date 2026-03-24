using ElectronKinetics, Trixi, GLMakie

function main()
    project_root = normpath(joinpath(@__DIR__, ".."))
    mesh_path = joinpath(project_root, "projects", "square_bells_ucsb", "mesh", "square_bells.inp")

    reference = ElectronKinetics.blg_reference_setup()
    mu0       = reference.mu0
    mass      = reference.mass
    gamma_mr  = reference.gamma_mr
    gamma_mc  = 100.0
    bias      = 0.5
    chi       = 0.0
    p_scatter = reference.p_scatter

    boundary_conditions = Dict(
        :walls          => MaxwellWallBC(p_scatter),
        :contact_top    => OhmicContactBC(-bias / 2),
        :contact_bottom => OhmicContactBC( bias / 2),
    )

    params = SolveParams(;
        polydeg           = 3,
        tspan_end         = 1e8,
        residual_tol      = 1e-3,
        cfl               = 0.8,
        log_every         = 500,
        min_harmonic      = 4,
        max_harmonic_auto = 20,
    )

    sol, semi = ElectronKinetics.solve(
        mesh_path,
        boundary_conditions,
        params,
        gamma_mr,
        gamma_mc;
        transport          = :parabolic_nonlinear,
        max_harmonic       = :auto,
        mu0                = mu0,
        mass               = mass,
        chi                = chi,
        visualize     = true,
        viz_field     = :a1,
        viz_colormap  = :RdBu,
        name               = "square_bells_a1_bias05",
    )

    return nothing
end

main()
