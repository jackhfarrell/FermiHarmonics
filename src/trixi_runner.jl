function solve(
    problem::TrixiProblem,
    model::KineticModel2D,
    config::SolverConfig;
    u0_override::Union{Nothing, AbstractVector}=nothing,
    visualize::Bool=false,
    visualize_every::Union{Nothing, Integer}=nothing,
    visualization_mode::Symbol=:cartesian,
    name::AbstractString="run",
)
    isfile(problem.mesh_path) || error("Mesh file not found: $(problem.mesh_path)")
    validate(config)

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

    stepsize_callback = Trixi.StepsizeCallback(cfl=config.cfl)
    steady_state_callback = Trixi.SteadyStateCallback(abstol=config.residual_tol, reltol=0.0)
    monitor = monitor_callback(config, semi)

    callbacks = Any[stepsize_callback, steady_state_callback, monitor]
    if visualize
        interval = isnothing(visualize_every) ? config.log_every : Int(visualize_every)
        interval > 0 || throw(ArgumentError("visualize_every must be positive"))
        push!(callbacks, visualization_callback(config, semi, name; interval=interval, mode=visualization_mode))
    end

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

    status = solve_status(sol, semi, config)
    @info "Solve complete" name=name stop_reason=status.stop_reason final_time=status.final_time target_final_time=status.target_final_time final_residual=status.final_residual tolerance=config.residual_tol converged=status.converged retcode=status.retcode successful=status.successful
    flush(stdout)
    flush(stderr)

    return sol, semi
end

function monitor_callback(config, semi)
    return SciMLBase.DiscreteCallback(
        (u, t, integrator) -> integrator.stats.naccept % config.log_every == 0,
        integrator -> begin
            du_ode = Trixi.get_du(integrator)
            integrator.f(du_ode, integrator.u, integrator.p, integrator.t)
            du = Trixi.wrap_array(du_ode, semi)
            residual = Trixi.residual_steady_state(du, semi.equations)
            @info "Progress" iter=integrator.stats.naccept t=round(integrator.t, digits=4) dt=round(integrator.dt, digits=6) residual=round(residual, sigdigits=3) tolerance=config.residual_tol
            flush(stdout)
            flush(stderr)
            nothing
        end;
        save_positions=(false, false),
    )
end

function solve_status(sol, semi, config::SolverConfig; time_atol::Real=1e-10)
    final_u_ode = similar(sol.u[end])
    sol.prob.f(final_u_ode, sol.u[end], sol.prob.p, sol.t[end])
    final_du = Trixi.wrap_array(final_u_ode, semi)
    final_residual = Trixi.residual_steady_state(final_du, semi.equations)
    converged = final_residual <= config.residual_tol
    hit_final_time = isapprox(sol.t[end], config.tspan_end; atol=time_atol, rtol=0.0)
    stop_reason = converged ? :steady_state : (hit_final_time ? :final_time : :other)
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

function visualization_callback(config, semi, name::AbstractString; interval::Int=config.log_every, mode::Symbol=:cartesian)
    if transport_is_nonlinear(semi.equations)
        return nonlinear_visualization_callback(config, semi, name; interval=interval, mode=mode)
    end

    variable_names = semi.equations isa MultiBandFermiHarmonics2D ?
        collect(Trixi.varnames(Trixi.cons2cons, semi.equations)) :
        ["a0", "a1", "b1"]
    return Trixi.VisualizationCallback(
        semi;
        interval=interval,
        variable_names=variable_names,
        filename="live_viz_$(name)",
        overwrite=true,
        seriescolor=:magma,
    )
end

function nonlinear_visualization_callback(config, semi, name::AbstractString; interval::Int=config.log_every, mode::Symbol=:cartesian)
    output_path = "live_viz_$(name).png"
    analysis_path = mode === :mesh_native ? "live_viz_$(name)_mesh_native.h5" : "live_viz_$(name).h5"
    project_root = normpath(joinpath(@__DIR__, ".."))
    plot_script = mode === :mesh_native ?
        joinpath(project_root, "demo", "plot_mesh_native_streamlines.py") :
        joinpath(project_root, "demo", "plot_nonlinear_streamlines.py")
    nvisnodes = 120

    return SciMLBase.DiscreteCallback(
        (u, t, integrator) -> integrator.stats.naccept % interval == 0,
        integrator -> begin
            if mode === :mesh_native
                mesh_data = compute_mesh_native_analysis(integrator.u, semi; refine=6)
                analysis_write_mesh_native_hdf5(analysis_path, mesh_data, integrator.t, semi.equations)
            elseif mode === :cartesian
                grids = compute_analysis_grids(integrator.u, semi; nvisnodes=nvisnodes)
                analysis_write_hdf5(
                    analysis_path,
                    grids.density,
                    grids.a1,
                    grids.b1,
                    grids.jx,
                    grids.jy,
                    grids.x,
                    grids.y,
                    grids.mask,
                    integrator.t,
                    grids.equations,
                    band_grids=grids.bands,
                )
            else
                throw(ArgumentError("unsupported visualization_mode=$(mode); use :cartesian or :mesh_native"))
            end
            run(`python3 $plot_script $analysis_path --output $output_path`)
            @info "Updated nonlinear live visualization" path=output_path t=round(integrator.t, digits=4)
            nothing
        end;
        save_positions=(false, false),
    )
end
