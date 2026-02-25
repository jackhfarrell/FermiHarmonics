# run a T-only parameter sweep for the simple_geometries rectangle problem using
# SLURM array jobs. this follows the same submit/worker structure as the existing
# gamma sweeps, with warm-starting across monotonically increasing T.

using Dates
using DrWatson
using FermiHarmonics
using HDF5
using SHA
using Sockets: gethostname
using Trixi

# ======================================================================================================================
# Configuration and Setup
# ======================================================================================================================

# sweep configuration
name = get(ENV, "FH_RECT_T_NAME", "simple_geometries_rectangle_T")
project_root = normpath(joinpath(@__DIR__, "..", ".."))
main_project = normpath(joinpath(project_root, "..", ".."))
mesh_path = get(ENV, "FH_RECT_T_MESH_PATH", joinpath(project_root, "meshes", "rectangle", "rectangle.inp"))
results_root = get(ENV, "FH_RECT_T_RESULTS_ROOT", joinpath(project_root, "results"))
mkpath(results_root)

# local mode: run all T points in-process without submitting slurm (default false)
run_local = lowercase(get(ENV, "FH_RECT_T_RUN_LOCAL", "false")) in ("1", "true", "yes", "on")

# physical and sweep parameters
bias = parse(Float64, get(ENV, "FH_RECT_T_BIAS", "1.0"))
p_scatter = parse(Float64, get(ENV, "FH_RECT_T_P_SCATTER", "1.0"))
a_mc = parse(Float64, get(ENV, "FH_RECT_T_A_MC", "1.0"))
a_mr = parse(Float64, get(ENV, "FH_RECT_T_A_MR", "1.0"))

Tmin = parse(Float64, get(ENV, "FH_RECT_T_TMIN", "1e-3"))
Tmax = parse(Float64, get(ENV, "FH_RECT_T_TMAX", "1e-1"))
nT = parse(Int, get(ENV, "FH_RECT_T_NT", "2500"))
chunk_size = parse(Int, get(ENV, "FH_RECT_T_CHUNK_SIZE", "20"))
nvisnodes = parse(Int, get(ENV, "FH_RECT_T_NVISNODES", "600"))

warm_start = lowercase(get(ENV, "FH_RECT_T_WARM_START", "true")) in ("1", "true", "yes", "on")
residual_mode = Symbol(lowercase(strip(get(ENV, "FH_RECT_T_RESIDUAL_MODE", "relative_low_modes"))))

# solver parameters
solve_params = SolveParams(;
    min_harmonic = parse(Int, get(ENV, "FH_RECT_T_MIN_HARMONIC", "4")),
    max_harmonic_auto = parse(Int, get(ENV, "FH_RECT_T_MAX_HARMONIC_AUTO", "100")),
    polydeg = parse(Int, get(ENV, "FH_RECT_T_POLYDEG", "3")),
    tspan_end = parse(Float64, get(ENV, "FH_RECT_T_TSPAN_END", "100.0")),
    residual_tol = parse(Float64, get(ENV, "FH_RECT_T_RESIDUAL_TOL", "1e-4")),
    residual_mode = residual_mode,
    residual_nvars = parse(Int, get(ENV, "FH_RECT_T_RESIDUAL_NVARS", "3")),
    residual_scale_floor = parse(Float64, get(ENV, "FH_RECT_T_RESIDUAL_SCALE_FLOOR", "1e-10")),
    cfl = parse(Float64, get(ENV, "FH_RECT_T_CFL", "0.5")),
    log_every = parse(Int, get(ENV, "FH_RECT_T_LOG_EVERY", "500")),
)

max_harmonic_env = lowercase(strip(get(ENV, "FH_RECT_T_MAX_HARMONIC", "auto")))
max_harmonic = max_harmonic_env == "auto" ? :auto : parse(Int, max_harmonic_env)

n_threads = parse(Int, get(ENV, "FH_RECT_T_N_THREADS", "1"))

# slurm parameters
sbatch = Dict(
    :job_name => "$(name)_sweep",
    :time => get(ENV, "FH_RECT_T_SBATCH_TIME", "24:00:00"),
    :cpus_per_task => parse(Int, get(ENV, "FH_RECT_T_SBATCH_CPUS_PER_TASK", "1")),
    :mem_per_cpu => get(ENV, "FH_RECT_T_SBATCH_MEM_PER_CPU", "4G"),
    :account => get(ENV, "FH_RECT_T_SBATCH_ACCOUNT", "ucb485_asc2"),
    :partition => get(ENV, "FH_RECT_T_SBATCH_PARTITION", "amilan"),
    :qos => get(ENV, "FH_RECT_T_SBATCH_QOS", "normal"),
    :export => "ALL",
)

function validate_T(T::Real)
    (0 < T < 1) || error("Invalid T=$(T). Expected strict 0 < T < 1.")
    return nothing
end

function gamma_from_T(T::Float64, a_mc::Float64, a_mr::Float64)
    validate_T(T)
    gamma_mc = a_mc * T^2 * log(1 / T)
    gamma_mr = a_mr * T^3
    return gamma_mr, gamma_mc
end

function build_T_grid(Tmin::Float64, Tmax::Float64, nT::Int)
    nT >= 1 || error("nT must be >= 1")
    (0 < Tmin < 1) || error("Tmin must satisfy 0 < Tmin < 1, got $(Tmin)")
    (0 < Tmax < 1) || error("Tmax must satisfy 0 < Tmax < 1, got $(Tmax)")
    Tmin < Tmax || error("Require Tmin < Tmax, got Tmin=$(Tmin), Tmax=$(Tmax)")
    T_vals = collect(10.0 .^ range(log10(Tmin), log10(Tmax), length=nT))
    foreach(validate_T, T_vals)
    return T_vals
end

function param_hash(; a_mc, a_mr, Tmin, Tmax, nT, chunk_size, bias, p_scatter)
    payload = "a_mc=$(a_mc);a_mr=$(a_mr);Tmin=$(Tmin);Tmax=$(Tmax);nT=$(nT);chunk=$(chunk_size);bias=$(bias);p=$(p_scatter)"
    return bytes2hex(sha1(payload))[1:10]
end

function ensure_summary_header(path::AbstractString)
    if !isfile(path)
        open(path, "w") do io
            write(io, "T,gamma_mc,gamma_mr,status,iters,residual,runtime,output_file\n")
        end
    end
end

function append_summary_row(path::AbstractString; T, gamma_mc, gamma_mr, status, iters, residual, runtime, output_file)
    open(path, "a") do io
        write(io, "$(T),$(gamma_mc),$(gamma_mr),$(status),$(iters),$(residual),$(runtime),$(output_file)\n")
    end
end

function save_warm_state(path::AbstractString, u_final::AbstractVector, T::Float64, gamma_mr::Float64, gamma_mc::Float64)
    h5open(path, "w") do f
        f["u_final"] = collect(u_final)
        attributes(f)["T"] = T
        attributes(f)["gamma_mr"] = gamma_mr
        attributes(f)["gamma_mc"] = gamma_mc
        attributes(f)["saved_at"] = string(now())
    end
end

function load_warm_state(path::AbstractString)
    h5open(path, "r") do f
        return read(f["u_final"])
    end
end

function nearest_smaller_state(state_dir::AbstractString, T_target::Float64)
    isdir(state_dir) || return nothing, NaN
    best_path = nothing
    best_T = -Inf
    for fp in readdir(state_dir; join=true)
        endswith(fp, ".h5") || continue
        m = match(r"T=([0-9eE+\-.]+)", basename(fp))
        m === nothing && continue
        T = tryparse(Float64, m.captures[1])
        T === nothing && continue
        if T < T_target && T > best_T
            best_T = T
            best_path = fp
        end
    end
    return best_path, best_T
end

function final_residual(sol, semi, params::SolveParams)
    tspan_end = params.tspan_end
    ode = Trixi.semidiscretize(semi, (0.0, tspan_end))
    du_ode = similar(sol.u[end])
    ode.f(du_ode, sol.u[end], ode.p, sol.t[end])
    du = Trixi.wrap_array(du_ode, semi)
    du_norm = Trixi.residual_steady_state(du, semi.equations)

    if params.residual_mode == :relative_low_modes
        u_local = Trixi.wrap_array(sol.u[end], semi)
        return FermiHarmonics.relative_residual_steady_state(
            du,
            u_local,
            semi.equations;
            scale_floor=params.residual_scale_floor,
        )
    end

    return du_norm
end

function accepted_steps(sol)
    if hasproperty(sol, :destats)
        stats = getproperty(sol, :destats)
        return hasproperty(stats, :naccept) ? getproperty(stats, :naccept) : missing
    elseif hasproperty(sol, :stats)
        stats = getproperty(sol, :stats)
        return hasproperty(stats, :naccept) ? getproperty(stats, :naccept) : missing
    end
    return missing
end

T_vals = build_T_grid(Tmin, Tmax, nT)
gamma_mr_vals = map(T -> gamma_from_T(T, a_mc, a_mr)[1], T_vals)
gamma_mc_vals = map(T -> gamma_from_T(T, a_mc, a_mr)[2], T_vals)

total_cases = length(T_vals)
n_jobs = cld(total_cases, chunk_size)
if total_cases % chunk_size != 0
    @warn "total_cases is not divisible by chunk_size; last task will run fewer cases" total_cases chunk_size n_jobs
end

boundary_conditions = Dict(
    :walls => MaxwellWallBC(p_scatter),
    :contact_top => OhmicContactBC(-bias / 2),
    :contact_bottom => OhmicContactBC(bias / 2),
)

slurm_metadata = Dict(
    "n_jobs" => n_jobs,
    "chunk_size" => chunk_size,
    "job_name" => sbatch[:job_name],
    "time" => sbatch[:time],
    "cpus_per_task" => sbatch[:cpus_per_task],
    "mem_per_cpu" => sbatch[:mem_per_cpu],
    "account" => sbatch[:account],
    "partition" => sbatch[:partition],
    "qos" => sbatch[:qos],
)

sweep_metadata = Dict(
    "sweep_variable" => "T",
    "geometry" => "rectangle",
    "mesh_file" => basename(mesh_path),
    "T_values" => T_vals,
    "Tmin" => Tmin,
    "Tmax" => Tmax,
    "nT" => nT,
    "a_mc" => a_mc,
    "a_mr" => a_mr,
    "gamma_mc_formula" => "a_mc * T^2 * log(1/T)",
    "gamma_mr_formula" => "a_mr * T^3",
    "bias" => bias,
    "p_scatter" => p_scatter,
    "warm_start" => warm_start,
    "max_harmonic" => string(max_harmonic),
    "nvisnodes_analysis" => nvisnodes,
)

# ======================================================================================================================
# Submit mode
# ======================================================================================================================

if !haskey(ENV, "SLURM_ARRAY_TASK_ID") && !run_local
    @info "Rectangle T sweep submission" name mesh=basename(mesh_path) Tmin Tmax nT chunk_size n_jobs a_mc a_mr

    timestamp = Dates.format(now(), "yyyy-mm-dd_HHMMSS")
    tag = param_hash(; a_mc, a_mr, Tmin, Tmax, nT, chunk_size, bias, p_scatter)
    sweep_id = "$(name)_sweep_rectangle_$(tag)_$(timestamp)"
    sweep_dir = joinpath(results_root, sweep_id)
    data_dir = joinpath(sweep_dir, "data")
    state_dir = joinpath(sweep_dir, "state")
    log_dir = joinpath(sweep_dir, "logs")
    summary_dir = joinpath(sweep_dir, "summary")
    mkpath(data_dir)
    mkpath(state_dir)
    mkpath(log_dir)
    mkpath(summary_dir)

    sbatch_submit = copy(sbatch)
    sbatch_submit[:output] = joinpath(log_dir, "slurm-%A_%a.out")
    sbatch_submit[:error] = joinpath(log_dir, "slurm-%A_%a.err")

    submit_sweep!(;
        script=abspath(@__FILE__),
        n_jobs=n_jobs,
        project_dir=main_project,
        env=Dict(
            "FERMI_SWEEP_DIR" => sweep_dir,
            "JULIA_NUM_THREADS" => string(n_threads),
        ),
        sbatch=sbatch_submit,
    )

    write_sweep_metadata!(sweep_dir, solve_params, sweep_metadata, slurm_metadata, gamma_mr_vals, gamma_mc_vals)
    archive_mesh!(mesh_path, sweep_dir)

    @info "Submitted" output_dir=sweep_dir
    flush(stdout)
    flush(stderr)
    exit(0)
end

# ======================================================================================================================
# Worker mode (or local mode)
# ======================================================================================================================

if haskey(ENV, "SLURM_ARRAY_TASK_ID")
    job_id = ENV["SLURM_ARRAY_JOB_ID"]
    task_id = parse(Int, ENV["SLURM_ARRAY_TASK_ID"])
    sweep_dir = get(ENV, "FERMI_SWEEP_DIR", joinpath(results_root, "$(name)_sweep_$(job_id)"))
    case_start_idx = (task_id - 1) * chunk_size + 1
    case_end_idx = min(task_id * chunk_size, total_cases)
    assigned_indices = case_start_idx:case_end_idx
else
    task_id = 1
    timestamp = Dates.format(now(), "yyyy-mm-dd_HHMMSS")
    tag = param_hash(; a_mc, a_mr, Tmin, Tmax, nT, chunk_size, bias, p_scatter)
    sweep_id = "$(name)_sweep_rectangle_$(tag)_$(timestamp)"
    sweep_dir = joinpath(results_root, sweep_id)
    assigned_indices = 1:total_cases
    mkpath(sweep_dir)
    write_sweep_metadata!(sweep_dir, solve_params, sweep_metadata, slurm_metadata, gamma_mr_vals, gamma_mc_vals)
    archive_mesh!(mesh_path, sweep_dir)
    @info "Running local sweep" output_dir=sweep_dir cases=length(assigned_indices)
end

data_dir = joinpath(sweep_dir, "data")
state_dir = joinpath(sweep_dir, "state")
summary_dir = joinpath(sweep_dir, "summary")
mkpath(data_dir)
mkpath(state_dir)
mkpath(summary_dir)

summary_path = haskey(ENV, "SLURM_ARRAY_TASK_ID") ?
    joinpath(summary_dir, "summary_task_$(lpad(string(task_id), 4, '0')).csv") :
    joinpath(summary_dir, "summary_all.csv")
ensure_summary_header(summary_path)

assigned = collect(assigned_indices)
if isempty(assigned)
    @warn "No assigned cases for task" task_id total_cases chunk_size
    exit(0)
end

@info "Worker starting" task_id hostname=gethostname() julia=VERSION threads=Threads.nthreads() cases=length(assigned)
preview_n = min(5, length(assigned))
preview = map(assigned[1:preview_n]) do idx
    T = T_vals[idx]
    gamma_mr, gamma_mc = gamma_from_T(T, a_mc, a_mr)
    (index = idx, T = T, gamma_mr = gamma_mr, gamma_mc = gamma_mc)
end
@info "Task ordered preview" task_id preview

local_mesh_path = haskey(ENV, "SLURM_ARRAY_TASK_ID") ? copy_mesh_to_scratch(mesh_path, "fh_$(name)_") : mesh_path

let
    u0_override = nothing
    task_start = now()

    for (case_counter, idx) in enumerate(assigned)
        T = T_vals[idx]
        gamma_mr, gamma_mc = gamma_from_T(T, a_mc, a_mr)

        if warm_start && isnothing(u0_override) && length(assigned) == 1
            prior_path, prior_T = nearest_smaller_state(state_dir, T)
            if !isnothing(prior_path)
                try
                    u0_override = load_warm_state(prior_path)
                    @info "Loaded disk warm start" T prior_T file=basename(prior_path)
                catch err
                    @warn "Failed to load disk warm start; using default initial guess" T err=sprint(showerror, err)
                    u0_override = nothing
                end
            end
        end

        restart = (!warm_start || isnothing(u0_override)) ? "cold" : "warm"
        @info "Starting T case" progress="$(case_counter)/$(length(assigned))" index=idx T gamma_mc gamma_mr restart

        case_start = now()
        status = "ok"
        iters = missing
        residual = NaN
        runtime_s = NaN
        output_file = ""

        try
            sol, semi = FermiHarmonics.solve(
                local_mesh_path,
                boundary_conditions,
                solve_params,
                gamma_mr,
                gamma_mc;
                max_harmonic=max_harmonic,
                u0_override=(warm_start ? u0_override : nothing),
                name="$(name)_T$(idx)",
            )

            iters = accepted_steps(sol)
            runtime_s = Dates.value(now() - case_start) / 1000.0
            residual = try
                final_residual(sol, semi, solve_params)
            catch err
                @warn "Residual evaluation failed" T err=sprint(showerror, err)
                NaN
            end

            file_params = (
                geometry = "rectangle",
                T = T,
                gamma_mc = gamma_mc,
                gamma_mr = gamma_mr,
                a_mc = a_mc,
                a_mr = a_mr,
                bias = bias,
                p_scatter = p_scatter,
            )
            analysis_path = joinpath(data_dir, DrWatson.savename(file_params; connector="_", sort=true) * ".h5")
            state_path = joinpath(state_dir, "state_" * DrWatson.savename(file_params; connector="_", sort=true) * ".h5")

            save_for_analysis(sol, semi, analysis_path; nvisnodes=nvisnodes)
            save_warm_state(state_path, sol.u[end], T, gamma_mr, gamma_mc)
            output_file = relpath(analysis_path, sweep_dir)

            @info "Case converged" T gamma_mc gamma_mr iterations=iters residual runtime_s output=basename(analysis_path)

            u0_override = warm_start ? copy(sol.u[end]) : nothing
        catch err
            status = "error"
            runtime_s = Dates.value(now() - case_start) / 1000.0
            residual = NaN
            output_file = ""
            u0_override = nothing
            @error "Case failed" T gamma_mc gamma_mr runtime_s exception=(err, catch_backtrace())
        end

        append_summary_row(
            summary_path;
            T=T,
            gamma_mc=gamma_mc,
            gamma_mr=gamma_mr,
            status=status,
            iters=iters,
            residual=residual,
            runtime=runtime_s,
            output_file=output_file,
        )

        GC.gc()
    end

    total_duration = Dates.canonicalize(Dates.CompoundPeriod(now() - task_start))
    total_ms = Dates.value(now() - task_start)
    avg_ms = round(Int, total_ms / max(1, length(assigned)))
    @info "Task complete" task_id cases=length(assigned) total_time=total_duration avg_per_case=Dates.canonicalize(
        Dates.CompoundPeriod(Dates.Millisecond(avg_ms))
    ) summary_file=summary_path
end
