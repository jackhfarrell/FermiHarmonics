# This file contains the main solver function and callbacks for the Fermi Harmonics problem.  The `solve` function is 
# the main entry point, which takes in a mesh file, boundary conditions, parameters, and options for visualization and 
# callbacks.  The callbacks include a monitor to track progress and a visualization callback to visualize the solution 
# during the solve.


# ======================================================================================================================
# Solver
# ======================================================================================================================

const AUTO_HARMONIC_GAMMA_HIGH = 300.0

# ======================================================================================================================
# Solver Configuration
# ======================================================================================================================

"""
    SolveParams

Typed configuration for [`solve`](@ref).

Fields:
- `max_harmonic::Int=60`: default fixed harmonic cutoff used when `solve(...; max_harmonic=...)` is not provided.
- `min_harmonic::Int=4`: lower bound for auto harmonic estimation.
- `max_harmonic_auto::Int=100`: upper bound for auto harmonic estimation.
- `polydeg::Int=3`: DGSEM polynomial degree.
- `tspan_end::Float64=100.0`: final integration time.
- `residual_tol::Float64=1e-5`: steady-state absolute tolerance.
- `cfl::Float64=0.8`: CFL number for timestep control.
- `log_every::Int=500`: monitor/visualization logging interval (accepted steps).
"""
Base.@kwdef struct SolveParams
    max_harmonic::Int = 60
    min_harmonic::Int = 4
    max_harmonic_auto::Int = 100
    polydeg::Int = 3
    tspan_end::Float64 = 100.0
    residual_tol::Float64 = 1e-5
    cfl::Float64 = 0.8
    log_every::Int = 500
end

function validate(params::SolveParams)
    params.max_harmonic >= 1 || throw(ArgumentError("params.max_harmonic must be >= 1"))
    params.min_harmonic >= 1 || throw(ArgumentError("params.min_harmonic must be >= 1"))
    params.max_harmonic_auto >= params.min_harmonic ||
        throw(ArgumentError("params.max_harmonic_auto must be >= params.min_harmonic"))
    params.polydeg >= 1 || throw(ArgumentError("params.polydeg must be >= 1"))
    params.tspan_end > 0 || throw(ArgumentError("params.tspan_end must be > 0"))
    params.residual_tol > 0 || throw(ArgumentError("params.residual_tol must be > 0"))
    params.cfl > 0 || throw(ArgumentError("params.cfl must be > 0"))
    params.log_every >= 1 || throw(ArgumentError("params.log_every must be >= 1"))
    return params
end

# ======================================================================================================================
# Harmonic Cutoff Selection
# ======================================================================================================================

"""
    estimate_max_harmonic(gamma_mr, gamma_mc; min_harmonic=4, max_harmonic=100)

Estimate an efficient harmonic cutoff from physical scattering rates.

The estimate uses the total scattering rate
`gamma_total = gamma_mr + gamma_mc` and logarithmically interpolates:
- `gamma_total = 0` -> `max_harmonic`,
- `gamma_total >= 300` -> `min_harmonic`,
- intermediate values map via `log1p(gamma_total)` between those endpoints.
"""
function estimate_max_harmonic(
    gamma_mr::Real,
    gamma_mc::Real;
    min_harmonic::Integer = 4,
    max_harmonic::Integer = 100,
)::Int
    gamma_mr < 0 && throw(ArgumentError("gamma_mr must be >= 0"))
    gamma_mc < 0 && throw(ArgumentError("gamma_mc must be >= 0"))

    min_h = Int(min_harmonic)
    max_h = Int(max_harmonic)
    min_h >= 1 || throw(ArgumentError("min_harmonic must be >= 1"))
    max_h >= min_h || throw(ArgumentError("max_harmonic must be >= min_harmonic"))

    gamma_total = Float64(gamma_mr) + Float64(gamma_mc)
    gamma_total <= 0 && return max_h
    gamma_total >= AUTO_HARMONIC_GAMMA_HIGH && return min_h

    frac = log1p(gamma_total) / log1p(AUTO_HARMONIC_GAMMA_HIGH)
    estimate = max_h - (max_h - min_h) * frac
    return clamp(ceil(Int, estimate), min_h, max_h)
end

"""
    resolve_max_harmonic(max_harmonic_kw, params, gamma_mr, gamma_mc)

Resolve the harmonic cutoff for one solve.

- `Integer`: use fixed cutoff.
- `:auto` or `nothing`: estimate from `gamma_mr`, `gamma_mc`.

Auto mode can be tuned through optional `params` fields:
- `min_harmonic` (default `4`)
- `max_harmonic_auto` (default `100`)
"""
function resolve_max_harmonic(max_harmonic_kw, params::SolveParams, gamma_mr::Real, gamma_mc::Real)
    if max_harmonic_kw isa Integer
        M = Int(max_harmonic_kw)
        M >= 1 || throw(ArgumentError("max_harmonic must be >= 1"))
        return M, :manual
    end

    if max_harmonic_kw === :auto || isnothing(max_harmonic_kw)
        min_h = params.min_harmonic
        max_h = params.max_harmonic_auto
        M = estimate_max_harmonic(
            gamma_mr,
            gamma_mc;
            min_harmonic=min_h,
            max_harmonic=max_h,
        )
        return M, :auto
    end

    throw(ArgumentError("max_harmonic must be an Integer, :auto, or nothing"))
end

# ======================================================================================================================
# Warm Start Utilities
# ======================================================================================================================

"""
    resize_warm_start(u0_override, target_u0, target_nvars)

Resize a warm-start vector so it matches the current solve size.

If harmonic count changes but mesh/discretization stay the same, this preserves low-order
harmonics and either pads missing higher modes with zeros or truncates extra modes.
"""
function resize_warm_start(
    u0_override::AbstractVector,
    target_u0::AbstractVector,
    target_nvars::Integer,
)
    target_nvars_int = Int(target_nvars)
    target_len = length(target_u0)
    target_len % target_nvars_int == 0 ||
        throw(ArgumentError("Target state length $target_len is incompatible with nvars=$target_nvars_int"))
    target_block = target_len ÷ target_nvars_int
    source_len = length(u0_override)

    if source_len == target_len
        return (
            u0 = collect(Float64, u0_override),
            mode = :same,
            source_nvars = target_nvars_int,
            target_nvars = target_nvars_int,
        )
    end

    # Preferred path: same mesh/polydeg, only number of harmonics changed.
    if source_len % target_block == 0
        source_nvars = source_len ÷ target_block
        source_state = reshape(collect(Float64, u0_override), source_nvars, target_block)
        target_state = zeros(Float64, target_nvars_int, target_block)
        ncopy = min(source_nvars, target_nvars_int)
        @views target_state[1:ncopy, :] .= source_state[1:ncopy, :]
        mode = source_nvars < target_nvars_int ? :padded : :truncated
        return (
            u0 = vec(target_state),
            mode = mode,
            source_nvars = source_nvars,
            target_nvars = target_nvars_int,
        )
    end

    # Fallback path when source/target discretizations are incompatible.
    resized = zeros(Float64, target_len)
    ncopy = min(source_len, target_len)
    @inbounds resized[1:ncopy] .= u0_override[1:ncopy]
    mode = source_len < target_len ? :padded_flat : :truncated_flat
    return (
        u0 = resized,
        mode = mode,
        source_nvars = nothing,
        target_nvars = target_nvars_int,
    )
end

# ======================================================================================================================
# Solve Entry Point
# ======================================================================================================================

"""
    solve(mesh_path::AbstractString, boundary_conditions::Dict{Symbol, Any},
          params, gamma_mr::Real, gamma_mc::Real; kwargs...)

Solve one FermiHarmonics case.

Arguments:
- `mesh_path`: path to mesh input file.
- `boundary_conditions`: boundary-condition map used by Trixi.
- `params`: solver configuration, as [`SolveParams`](@ref).
- `gamma_mr`, `gamma_mc`: physical scattering rates.

Keywords:
- `max_harmonic`: harmonic cutoff. Use `:auto` (default) to estimate from `gamma_mr`, `gamma_mc`,
  or pass an integer for a fixed cutoff (`nvars = 1 + 2*max_harmonic`).
- `u0_override`: optional warm-start state vector.
- `visualize`: enable live visualization callback.
- `name`: run name used in logs/visualization filenames.

Returns:
- `(sol, semi)`: Trixi time-integration solution and semidiscretization.
"""
function solve(mesh_path::AbstractString, boundary_conditions::Dict{Symbol, Any},
               params::SolveParams, gamma_mr::Real, gamma_mc::Real;
               max_harmonic::Union{Integer, Symbol, Nothing}=:auto,
               u0_override::Union{Nothing, AbstractVector}=nothing,
               visualize::Bool=false,
               name::AbstractString="run")

    # ------------------------------------------------------------------------------------------------------------------
    # Input validation and equations setup
    # ------------------------------------------------------------------------------------------------------------------
    isfile(mesh_path) || error("Mesh file not found: $mesh_path")
    validate(params)

    max_harmonic_resolved, harmonic_mode = resolve_max_harmonic(max_harmonic, params, gamma_mr, gamma_mc)
    nvars = 1 + 2 * max_harmonic_resolved
    equations = FermiHarmonics2D(
        nvars;
        gamma_mr=gamma_mr,
        gamma_mc=gamma_mc,
        max_harmonic=max_harmonic_resolved,
    )

    boundary_symbols = sort(collect(keys(boundary_conditions)))
    solver = Trixi.DGSEM(polydeg=params.polydeg, surface_flux=Trixi.flux_lax_friedrichs)
    mesh = Trixi.P4estMesh{2}(mesh_path; boundary_symbols=boundary_symbols)
    semi = Trixi.SemidiscretizationHyperbolic(
        mesh, equations, (x, t, eq) -> zeros(SVector{nvars, Float64}), solver;
        boundary_conditions=boundary_conditions,
        source_terms=FermiHarmonics.source_terms,
    )

    boundary_types = Dict(key => boundary_condition_name(value) for (key, value) in boundary_conditions)
    @info "Starting solve" name=name max_harmonic=max_harmonic_resolved harmonic_mode=harmonic_mode gamma_mr=gamma_mr gamma_mc=gamma_mc polydeg=params.polydeg cfl=params.cfl residual_tol=params.residual_tol boundaries=boundary_types
    flush(stdout)
    flush(stderr)

    # ------------------------------------------------------------------------------------------------------------------
    # ODE construction and optional warm start
    # ------------------------------------------------------------------------------------------------------------------
    tspan = (0.0, params.tspan_end)
    ode = Trixi.semidiscretize(semi, tspan)

    # handle warm start from previous solution in memory if provided.
    if !isnothing(u0_override)
        warm = resize_warm_start(u0_override, ode.u0, nvars)
        if warm.mode != :same
            @info "Adjusted warm start for harmonic mismatch" mode=warm.mode source_nvars=warm.source_nvars target_nvars=warm.target_nvars source_length=length(u0_override) target_length=length(warm.u0)
            flush(stdout)
            flush(stderr)
        end
        ode = SciMLBase.remake(ode; u0=warm.u0)
    end

    # ------------------------------------------------------------------------------------------------------------------
    # Callback assembly
    # ------------------------------------------------------------------------------------------------------------------
    stepsize_callback = Trixi.StepsizeCallback(cfl=params.cfl)
    steady_state_callback = Trixi.SteadyStateCallback(abstol=params.residual_tol, reltol=0.0)
    monitor = monitor_callback(params, semi)

    callbacks = Any[stepsize_callback, steady_state_callback, monitor]

    if visualize
        push!(callbacks, visualization_callback(params, semi, name))
    end

    # ------------------------------------------------------------------------------------------------------------------
    # Time integration
    # ------------------------------------------------------------------------------------------------------------------
    # Solve the semidiscretized system with the configured explicit RK method.
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

    return sol, semi
end


# ======================================================================================================================
# Callbacks
# ======================================================================================================================

# ======================================================================================================================
# Progress Monitoring
# ======================================================================================================================

"""
    monitor_callback(params, semi)

Callback to monitor progress of the solve, printing out the current time, timestep, and residual every 
`params.log_every` steps.  Should be used in conjunction with `Trixi.SteadyStateCallback` to track convergence to steady 
state.
"""
function monitor_callback(params, semi)
    return SciMLBase.DiscreteCallback(
        (u, t, integrator) -> integrator.stats.naccept % params.log_every == 0,
        integrator -> begin
            du_ode = Trixi.get_du(integrator)
            integrator.f(du_ode, integrator.u, integrator.p, integrator.t)
            du = Trixi.wrap_array(du_ode, semi)
            residual = Trixi.residual_steady_state(du, semi.equations)
            @info "Progress" iter=integrator.stats.naccept t=round(integrator.t, digits=4) dt=round(integrator.dt, digits=6) residual=round(residual, sigdigits=3) tolerance=params.residual_tol
            flush(stdout)
            flush(stderr)
            nothing
        end;
        save_positions=(false, false),
    )
end

# ======================================================================================================================
# Live Visualization
# ======================================================================================================================

"""
    visualization_callback(params, semi, name::AbstractString)

Create a live visualization callback for current magnitude
`j_norm = sqrt(a1^2 + b1^2)` every `params.log_every` accepted steps.
"""
function visualization_callback(params, semi, name::AbstractString)
    return Trixi.VisualizationCallback(
        semi;
        interval=params.log_every,
        solution_variables=current_norm_variables,
        variable_names=["j_norm"],
        filename="live_viz_$(name)",
        overwrite=true,
        seriescolor=:magma,
    )
end
