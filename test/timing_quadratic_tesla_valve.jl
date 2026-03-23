using ElectronKinetics, Trixi

function env_float(name::AbstractString, default::Float64)
    return haskey(ENV, name) ? parse(Float64, ENV[name]) : default
end

function env_int(name::AbstractString, default::Int)
    return haskey(ENV, name) ? parse(Int, ENV[name]) : default
end

function main()
    project_root = normpath(joinpath(@__DIR__, ".."))
    mesh_path = joinpath(project_root, "projects", "nonlinearities", "mesh", "tesla_valve.inp")

    reference = ElectronKinetics.blg_reference_setup()
    mu0 = reference.mu0
    mass = reference.mass
    gamma_mr = env_float("FH_TIMING_GAMMA_MR", 0.01)
    gamma_mc = env_float("FH_TIMING_GAMMA_MC", 0.01)
    bias = env_float("FH_TIMING_BIAS", 0.2)
    chi = env_float("FH_TIMING_CHI", 0.0)
    max_harmonic = env_int("FH_TIMING_MAX_HARMONIC", 50)
    tspan_end = env_float("FH_TIMING_TSPAN_END", 0.1)
    residual_tol = env_float("FH_TIMING_RESIDUAL_TOL", 1.0e-3)
    cfl = env_float("FH_TIMING_CFL", 0.8)
    polydeg = env_int("FH_TIMING_POLYDEG", 3)

    boundary_conditions = Dict(
        :walls => MaxwellWallBC(reference.p_scatter),
        :inlet => OhmicContactBC(bias / 2),
        :outlet => OhmicContactBC(-bias / 2),
    )

    params = SolveParams(;
        polydeg=polydeg,
        tspan_end=tspan_end,
        residual_tol=residual_tol,
        cfl=cfl,
        log_every=10_000,
        min_harmonic=min(8, max_harmonic),
        max_harmonic_auto=max_harmonic,
    )

    println("Quadratic Tesla timing run")
    println("mesh_path         = ", mesh_path)
    println("gamma_mr         = ", gamma_mr)
    println("gamma_mc         = ", gamma_mc)
    println("bias             = ", bias)
    println("chi              = ", chi)
    println("max_harmonic     = ", max_harmonic)
    println("polydeg          = ", polydeg)
    println("tspan_end        = ", tspan_end)
    println("cfl              = ", cfl)
    println("residual_tol     = ", residual_tol)
    println("threads          = ", Threads.nthreads())

    reset_nonlinear_timing!()
    enable_nonlinear_timing!()
    wall_t0 = time_ns()
    sol, semi = solve(
        mesh_path,
        boundary_conditions,
        params,
        gamma_mr,
        gamma_mc;
        transport=:parabolic_nonlinear,
        max_harmonic=max_harmonic,
        mu0=mu0,
        mass=mass,
        chi=chi,
        visualize=false,
        name="timing_quadratic_tesla_valve",
    )
    wall_ns = time_ns() - wall_t0
    disable_nonlinear_timing!()

    println()
    println("Solve complete")
    println("final_time        = ", sol.t[end])
    println("wall_time_s       = ", round(Float64(wall_ns) * 1.0e-9; digits=4))
    print_nonlinear_timing_summary()

    snapshot = nonlinear_timing_snapshot()
    instrumented_s = Float64(snapshot.total_ns) * 1.0e-9
    wall_s = Float64(wall_ns) * 1.0e-9
    uncovered = max(wall_s - instrumented_s, 0.0)
    println("uncovered_wall_s  = ", round(uncovered; digits=4))

    return nothing
end

main()
