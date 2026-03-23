function solve(
    problem::TrixiProblem,
    model::KineticModel2D,
    config::SolverConfig;
    u0_override::Union{Nothing, AbstractVector}=nothing,
    live_visualization::Union{Nothing, LiveVisualizationConfig}=nothing,
    visualize::Bool=false,
    visualize_every::Union{Nothing, Integer}=nothing,
    visualization_mode::Symbol=:mesh_native,
    name::AbstractString="run",
)
    isfile(problem.mesh_path) || error("Mesh file not found: $(problem.mesh_path)")
    validate(config)
    if isnothing(live_visualization) && visualize
        interval = isnothing(visualize_every) ? config.log_every : Int(visualize_every)
        live_visualization = LiveVisualizationConfig(;
            accepted_step_interval=interval,
            geometry_mode=visualization_mode,
        )
    end
    !isnothing(live_visualization) && validate(live_visualization)

    equations, harmonic_mode = build_equations(model, config)
    nvars = Trixi.nvariables(equations)
    boundary_symbols = sort(collect(keys(problem.boundary_conditions)))
    boundary_conditions = (; problem.boundary_conditions...)
    solver = Trixi.DGSEM(polydeg=config.polydeg, surface_flux=Trixi.flux_lax_friedrichs)
    mesh = Trixi.P4estMesh{2}(problem.mesh_path; boundary_symbols=boundary_symbols)

    if equations isa FermiHarmonics2D && transport_is_nonlinear(equations)
        semi = Trixi.SemidiscretizationHyperbolic(
            mesh, equations, (x, t, eq) -> zeros(SVector{nvars, Float64}), solver;
            boundary_conditions=boundary_conditions,
            source_terms=source_terms,
        )
    elseif nonlinear_has_electrostatic_force(equations)
        equations_parabolic = ElectrostaticGradientEquation2D(equations)
        semi = Trixi.SemidiscretizationHyperbolicParabolic(
            mesh,
            (equations, equations_parabolic),
            (x, t, eq) -> zeros(SVector{nvars, Float64}),
            solver;
            solver_parabolic=Trixi.ViscousFormulationLocalDG(),
            source_terms=source_terms,
            source_terms_parabolic=source_terms,
            boundary_conditions=(boundary_conditions, boundary_conditions),
        )
    else
        semi = Trixi.SemidiscretizationHyperbolic(
            mesh, equations, (x, t, eq) -> zeros(SVector{nvars, Float64}), solver;
            boundary_conditions=boundary_conditions,
            source_terms=source_terms,
        )
    end

    boundary_types = Dict(key => boundary_condition_name(value) for (key, value) in problem.boundary_conditions)
    @info "Starting solve" name=name harmonic_mode=harmonic_mode polydeg=config.polydeg cfl=config.cfl residual_tol=config.residual_tol boundaries=boundary_types transport=transport_symbol(model) collision_model=collision_symbol(model.collision)
    flush(stdout)
    flush(stderr)

    tspan = (0.0, config.tspan_end)
    ode = Trixi.semidiscretize(semi, tspan)

    if !isnothing(u0_override)
        warm = if equations isa FermiAngles2D
            validate_nonlinear_warm_start(u0_override, ode.u0, nvars)
        elseif equations isa MultiBandFermiHarmonics2D
            resize_multiband_warm_start(u0_override, ode.u0, band_count(equations), band_nvars(equations))
        else
            resize_warm_start(u0_override, ode.u0, nvars)
        end
        if warm.mode != :same
            @info "Adjusted warm start for state mismatch" mode=warm.mode source_nvars=warm.source_nvars target_nvars=warm.target_nvars source_length=length(u0_override) target_length=length(warm.u0)
            flush(stdout)
            flush(stderr)
        end
        ode = SciMLBase.remake(ode; u0=warm.u0)
    end

    monitor_state = create_monitor_state(ode, semi, live_visualization)
    initialize_live_dashboard!(monitor_state, ode.u0, semi, config, name)

    stepsize_callback = Trixi.StepsizeCallback(cfl=config.cfl)
    steady_state_callback = Trixi.SteadyStateCallback(abstol=config.residual_tol, reltol=0.0)
    monitor = solve_monitor_callback(config, semi, monitor_state)

    callbacks = Any[stepsize_callback, steady_state_callback, monitor]

    sol = Trixi.solve(
        ode,
        Trixi.CarpenterKennedy2N54();
        dt=stepsize_callback(ode),
        callback=Trixi.CallbackSet(callbacks...),
        adaptive=false,
        save_everystep=false,
        save_start=false,
        save_end=true,
    )

    status = solve_status(sol, semi, config; stop_reason_override=monitor_state.window_closed ? :window_closed : nothing)
    finalize_live_dashboard!(monitor_state, sol.u[end], semi, status, config)
    @info "Solve complete" name=name stop_reason=status.stop_reason final_time=status.final_time target_final_time=status.target_final_time final_residual=status.final_residual tolerance=config.residual_tol converged=status.converged retcode=status.retcode successful=status.successful
    flush(stdout)
    flush(stderr)

    return sol, semi
end

function solve_status(sol, semi, config::SolverConfig; time_atol::Real=1e-10, stop_reason_override::Union{Nothing, Symbol}=nothing)
    final_u_ode = similar(sol.u[end])
    sol.prob.f(final_u_ode, sol.u[end], sol.prob.p, sol.t[end])
    final_du = Trixi.wrap_array(final_u_ode, semi)
    final_residual = Trixi.residual_steady_state(final_du, semi.equations)
    converged = final_residual <= config.residual_tol
    hit_final_time = isapprox(sol.t[end], config.tspan_end; atol=time_atol, rtol=0.0)
    stop_reason = isnothing(stop_reason_override) ? (converged ? :steady_state : (hit_final_time ? :final_time : :other)) : stop_reason_override
    retcode = hasproperty(sol, :retcode) ? getproperty(sol, :retcode) : nothing
    successful = isnothing(retcode) ? true : SciMLBase.successful_retcode(retcode)
    return (
        stop_reason=stop_reason,
        final_residual=final_residual,
        converged=converged,
        hit_final_time=hit_final_time,
        final_time=sol.t[end],
        target_final_time=config.tspan_end,
        successful=successful,
        retcode=retcode,
    )
end
