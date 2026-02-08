# run a large parameter sweep for the square bells problem, sweeping over gamma_mr and gamma_mc using SLURM array jobs.
# this script handles both the initial submission (writing shared metadata and submitting the array job) and the worker
# logic. the parameter grid is defined in the script, and each worker picks which cases to run based on its
# SLURM_ARRAY_TASK_ID. results are saved in a shared directory with one file per case.

using Dates
using DrWatson
using FermiHarmonics
using Sockets: gethostname

# ======================================================================================================================
# CLI Arguments
# ======================================================================================================================

function parse_command_line_args(args::Vector{String})
    p_scatter = 1.0
    i = 1
    while i <= length(args)
        arg = args[i]
        value_str = nothing
        if startswith(arg, "--p_scatter=")
            value_str = split(arg, "=", limit=2)[2]
            i += 1
        elseif arg == "--p_scatter"
            i < length(args) || throw(ArgumentError("Missing value after --p_scatter"))
            value_str = args[i + 1]
            i += 2
        else
            throw(ArgumentError("Unknown argument: $arg. Supported: --p_scatter <float>"))
        end

        p_val = tryparse(Float64, value_str)
        p_val === nothing && throw(ArgumentError("Could not parse --p_scatter value '$value_str' as Float64"))
        (0.0 <= p_val <= 1.0) || throw(ArgumentError("--p_scatter must satisfy 0.0 <= p_scatter <= 1.0"))
        p_scatter = p_val
    end

    return (p_scatter = p_scatter,)
end

function p_scatter_tag(p_scatter::Real)
    return "p" * replace(string(Float64(p_scatter)), "." => "p", "-" => "m")
end

function ordered_case_indices(gamma_mr_vals, gamma_mc_vals)
    total_cases = length(gamma_mr_vals) * length(gamma_mc_vals)
    indices = collect(1:total_cases)
    sort!(
        indices;
        by = index_global -> begin
            gamma_mr, gamma_mc = grid_lookup(gamma_mr_vals, gamma_mc_vals, index_global)
            (gamma_mr + gamma_mc, gamma_mr, gamma_mc)
        end,
        rev = true,
    )
    return indices
end

# ======================================================================================================================
# Configuration and Setup
# ======================================================================================================================

# sweep configuration
name = "square_bells_ucsb"
project_root = normpath(joinpath(@__DIR__, ".."))
main_project = normpath(joinpath(project_root, "..", ".."))
mesh_path = joinpath(project_root, "mesh", "square_bells.inp")
results_root = joinpath(project_root, "results")
n_jobs = 250
cases_per_job = 10
mkpath(results_root)

cli_args = parse_command_line_args(collect(ARGS))
bias = 1.0 # change in a0 to drive flow
p_scatter = cli_args.p_scatter # scattering prob. at walls (1.0 = fully diffuse, 0.0 = fully specular)
p_tag = p_scatter_tag(p_scatter)

# slurm parameters
sbatch = Dict(
    :job_name => "$(name)_sweep_$(p_tag)",
    :time => "24:00:00",
    :cpus_per_task => 1,
    :mem_per_cpu => "4G",
    :account => "ucb485_asc2",
    :partition => "amilan",
    :qos => "normal",
    :export => "ALL",
)

# grid of parameters
gamma_mr_vals = 10 .^ range(log10(0.05), log10(100.0), length=50)
gamma_mc_vals = 10 .^ range(log10(0.1), log10(200.0), length=50)
total_cases = length(gamma_mr_vals) * length(gamma_mc_vals)
@assert length(gamma_mr_vals) * length(gamma_mc_vals) == n_jobs * cases_per_job

boundary_conditions = Dict(
    :walls => MaxwellWallBC(p_scatter),
    :contact_top => OhmicContactBC(-bias / 2),
    :contact_bottom => OhmicContactBC(bias / 2),
)

solve_params = SolveParams(;
    min_harmonic = 4,
    max_harmonic_auto = 100,
    polydeg = 3,
    tspan_end = 100.0, # end time for simulation (if not converged earlier)
    residual_tol = 1e-4,
    cfl = 0.5,
    log_every = 500,
)

slurm_metadata = Dict(
    "n_jobs" => n_jobs,
    "cases_per_job" => cases_per_job,
    "job_name" => sbatch[:job_name],
    "time" => sbatch[:time],
    "cpus_per_task" => sbatch[:cpus_per_task],
    "mem_per_cpu" => sbatch[:mem_per_cpu],
    "account" => sbatch[:account],
    "partition" => sbatch[:partition],
    "qos" => sbatch[:qos],
)
sweep_metadata = Dict(
    "bias" => bias,
    "p_scatter" => p_scatter,
)


# ======================================================================================================================
# Submit mode
# ======================================================================================================================

# if not running as part of SLURM array job (initial submission), write shared metadata and submit the job
if !haskey(ENV, "SLURM_ARRAY_TASK_ID")

    @info "Square bells sweep submission" name mesh=basename(mesh_path) bias p_scatter total_cases

    timestamp = Dates.format(now(), "yyyy-mm-dd_HHMMSS")
    sweep_id = "$(name)_sweep_$(p_tag)_$(timestamp)"
    sweep_dir = joinpath(results_root, sweep_id)
    data_dir = joinpath(sweep_dir, "data")
    log_dir = joinpath(sweep_dir, "logs")
    mkpath(data_dir)
    mkpath(log_dir)

    # configure directory for logs 
    sbatch_submit = copy(sbatch)
    sbatch_submit[:output] = joinpath(log_dir, "slurm-%A_%a.out")
    sbatch_submit[:error] = joinpath(log_dir, "slurm-%A_%a.err")

    # submit the sweep job
    submit_sweep!(;
        script=abspath(@__FILE__),
        n_jobs=n_jobs,
        project_dir=main_project,
        script_args=["--p_scatter=$(p_scatter)"],
        env=Dict("FERMI_SWEEP_DIR" => sweep_dir),
        sbatch=sbatch_submit,
    )

    # write metadata (always overwritten on submit) and archive mesh at sweep root
    write_sweep_metadata!(sweep_dir, solve_params, sweep_metadata, slurm_metadata, gamma_mr_vals, gamma_mc_vals)
    archive_mesh!(mesh_path, sweep_dir)

    @info "Submitted" output_dir=sweep_dir
    flush(stdout)
    flush(stderr)
    exit(0)
end

# ======================================================================================================================
# Worker mode
# ======================================================================================================================

# if running as part of SLURM array job, execute the assigned cases

# read metadata to get the parameter grid and other info
job_id = ENV["SLURM_ARRAY_JOB_ID"]
task_id = parse(Int, ENV["SLURM_ARRAY_TASK_ID"])
sweep_dir = get(ENV, "FERMI_SWEEP_DIR", joinpath(results_root, "$(name)_sweep_$(p_tag)_$(job_id)"))
data_dir = joinpath(sweep_dir, "data")
@info "Worker starting" job_id task_id hostname=gethostname() julia=VERSION threads=Threads.nthreads()

# pick which cases to run based on task ID and cases per task, and log the assignment
case_positions, total_cases, n_cases_this_task, _ = select_cases(
    gamma_mr_vals, gamma_mc_vals; cases_per_task=cases_per_job, env=ENV
)
ordered_global_indices = ordered_case_indices(gamma_mr_vals, gamma_mc_vals)
assigned_global_indices = ordered_global_indices[collect(case_positions)]
@info "Task assignment" task_id cases=n_cases_this_task position_range=(first(case_positions), last(case_positions))

preview_n = min(5, length(assigned_global_indices))
preview = map(assigned_global_indices[1:preview_n]) do idx
    gamma_mr, gamma_mc = grid_lookup(gamma_mr_vals, gamma_mc_vals, idx)
    (index = idx, gamma_mr = gamma_mr, gamma_mc = gamma_mc, gamma_total = gamma_mr + gamma_mc)
end
@info "Task ordered preview" task_id preview

# avoid data races by copying the mesh to a local scratch directory for this task
local_mesh_path = copy_mesh_to_scratch(mesh_path, "fh_$(name)_")

# loop over assigned cases, solve, save results, and use warm_starting within each block
let
    u0_override = nothing
    task_start = now()

    for (case_counter, index_global) in enumerate(assigned_global_indices)
        gamma_mr, gamma_mc = grid_lookup(gamma_mr_vals, gamma_mc_vals, index_global)
        gamma_total = gamma_mr + gamma_mc

        # determine cold or warm start
        restart = isnothing(u0_override) ? "cold" : "warm"
        @info "Starting case" progress="$(case_counter)/$(n_cases_this_task)" index=index_global gamma_mr gamma_mc gamma_total restart

        case_start = now()

        # solve the case!
        sol, semi = FermiHarmonics.solve(
            local_mesh_path, boundary_conditions, solve_params, gamma_mr, gamma_mc;
            max_harmonic=:auto,
            u0_override=u0_override,
            name="$(name)_g$(index_global)",
        )
        @info "Case converged" iterations=sol.destats.naccept duration=Dates.canonicalize(
            Dates.CompoundPeriod(now() - case_start)
        )
        
        # save a0, a1, b1 for analysis
        file_params = (bias=bias, p_scatter=p_scatter, gamma_mr=gamma_mr, gamma_mc=gamma_mc)
        small_path = joinpath(data_dir, DrWatson.savename(file_params; connector="_", sort=true) * ".h5")
        save_for_analysis(sol, semi, small_path)
        @info "Saved" path=basename(small_path) small_mb=round(filesize(small_path) / 1e6, digits=2)

        # set up warm start for next case in block
        u0_override = copy(sol.u[end])

        # collect garbage
        GC.gc()
    end

    # some nice logging
    total_duration = Dates.canonicalize(Dates.CompoundPeriod(now() - task_start))
    total_ms = Dates.value(now() - task_start)
    avg_ms = round(Int, total_ms / n_cases_this_task)
    @info "Task complete" task_id cases=n_cases_this_task total_time=total_duration avg_per_case=Dates.canonicalize(
        Dates.CompoundPeriod(Dates.Millisecond(avg_ms))
    )
end
