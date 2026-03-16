# run a 2D parameter sweep for the simple_geometries dogleg problem, sweeping over gamma_mr and gamma_mc using
# SLURM array jobs. this script handles both the initial submission (writing shared metadata and submitting the array
# job) and the worker logic. the parameter grid is defined in the script, and each worker picks which cases to run
# based on its SLURM_ARRAY_TASK_ID. results are saved in a shared directory with one file per case.

using Dates
using DrWatson
using FermiHarmonics
using Sockets: gethostname

parse_env_float(name::AbstractString, default::Float64) =
    haskey(ENV, name) ? parse(Float64, ENV[name]) : default


# ======================================================================================================================
# Configuration and Setup
# ======================================================================================================================

name = "simple_geometries_dogleg_gamma_grid"
project_root = normpath(joinpath(@__DIR__, "..", ".."))
main_project = normpath(joinpath(project_root, "..", ".."))
mesh_path = joinpath(project_root, "meshes", "dogleg", "dogleg.inp")
results_root = joinpath(project_root, "results")
n_jobs = 25
cases_per_job = 20
mkpath(results_root)

bias = 1.0
p_scatter = parse_env_float("FERMI_P_SCATTER", 1.0)

sbatch = Dict(
    :job_name => "$(name)_sweep",
    :time => "24:00:00",
    :cpus_per_task => 1,
    :mem_per_cpu => "4G",
    :account => "ucb485_asc2",
    :partition => "amilan",
    :qos => "normal",
    :export => "ALL",
)

n_gamma_mr = 10
n_gamma_mc = 50
gamma_mr_vals = 10 .^ range(log10(1e-2), log10(1e1), length=n_gamma_mr)
gamma_mc_vals = 10 .^ range(log10(1e-2), log10(1e2), length=n_gamma_mc)
total_cases = length(gamma_mr_vals) * length(gamma_mc_vals)
@assert total_cases == n_jobs * cases_per_job

boundary_conditions = Dict(
    :walls => MaxwellWallBC(p_scatter),
    :source => OhmicContactBC(bias / 2),
    :drain => OhmicContactBC(-bias / 2),
)

solve_params = SolveParams(;
    min_harmonic = 4,
    max_harmonic_auto = 100,
    polydeg = 3,
    tspan_end = 100.0,
    residual_tol = 1e-5,
    residual_mode = :absolute_all,
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
    "gamma_mr_min" => first(gamma_mr_vals),
    "gamma_mr_max" => last(gamma_mr_vals),
    "n_gamma_mr" => n_gamma_mr,
    "gamma_mc_min" => first(gamma_mc_vals),
    "gamma_mc_max" => last(gamma_mc_vals),
    "n_gamma_mc" => n_gamma_mc,
)


# ======================================================================================================================
# Submit mode
# ======================================================================================================================

if !haskey(ENV, "SLURM_ARRAY_TASK_ID")

    @info "Dogleg gamma-grid sweep submission" name mesh=basename(mesh_path) bias p_scatter total_cases n_gamma_mr n_gamma_mc

    timestamp = Dates.format(now(), "yyyy-mm-dd_HHMMSS")
    sweep_id = "$(name)_sweep_$(timestamp)"
    sweep_dir = joinpath(results_root, sweep_id)
    data_dir = joinpath(sweep_dir, "data")
    log_dir = joinpath(sweep_dir, "logs")
    mkpath(data_dir)
    mkpath(log_dir)

    sbatch_submit = copy(sbatch)
    sbatch_submit[:output] = joinpath(log_dir, "slurm-%A_%a.out")
    sbatch_submit[:error] = joinpath(log_dir, "slurm-%A_%a.err")

    submit_sweep!(;
        script=abspath(@__FILE__),
        n_jobs=n_jobs,
        project_dir=main_project,
        env=Dict(
            "FERMI_SWEEP_DIR" => sweep_dir,
            "FERMI_P_SCATTER" => string(p_scatter),
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
# Worker mode
# ======================================================================================================================

job_id = ENV["SLURM_ARRAY_JOB_ID"]
task_id = parse(Int, ENV["SLURM_ARRAY_TASK_ID"])
sweep_dir = get(ENV, "FERMI_SWEEP_DIR", joinpath(results_root, "$(name)_sweep_$(job_id)"))
data_dir = joinpath(sweep_dir, "data")
@info "Worker starting" job_id task_id hostname=gethostname() julia=VERSION threads=Threads.nthreads()

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

local_mesh_path = copy_mesh_to_scratch(mesh_path, "fh_$(name)_")

let
    u0_override = nothing
    task_start = now()

    for (case_counter, index_global) in enumerate(assigned_global_indices)
        gamma_mr, gamma_mc = grid_lookup(gamma_mr_vals, gamma_mc_vals, index_global)
        gamma_total = gamma_mr + gamma_mc

        restart = isnothing(u0_override) ? "cold" : "warm"
        @info "Starting case" progress="$(case_counter)/$(n_cases_this_task)" index=index_global gamma_mr gamma_mc gamma_total restart

        case_start = now()

        sol, semi = FermiHarmonics.solve(
            local_mesh_path, boundary_conditions, solve_params, gamma_mr, gamma_mc;
            max_harmonic=:auto,
            u0_override=u0_override,
            name="$(name)_g$(index_global)",
        )
        @info "Case converged" iterations=sol.destats.naccept duration=Dates.canonicalize(
            Dates.CompoundPeriod(now() - case_start)
        )

        file_params = (bias = bias, p_scatter = p_scatter, gamma_mr = gamma_mr, gamma_mc = gamma_mc)
        small_path = joinpath(data_dir, DrWatson.savename(file_params; connector="_", sort=true) * ".h5")
        save_for_analysis(sol, semi, small_path)
        @info "Saved" path=basename(small_path) small_mb=round(filesize(small_path) / 1e6, digits=2)

        u0_override = copy(sol.u[end])
        GC.gc()
    end

    total_duration = Dates.canonicalize(Dates.CompoundPeriod(now() - task_start))
    total_ms = Dates.value(now() - task_start)
    avg_ms = round(Int, total_ms / n_cases_this_task)
    @info "Task complete" task_id cases=n_cases_this_task total_time=total_duration avg_per_case=Dates.canonicalize(
        Dates.CompoundPeriod(Dates.Millisecond(avg_ms))
    )
end
