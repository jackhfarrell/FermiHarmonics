using ElectronKinetics, Trixi

function tesla_cluster_required_env(primary::AbstractString, fallback::AbstractString="")
    value = get(ENV, primary, "")
    isempty(value) && !isempty(fallback) && (value = get(ENV, fallback, ""))
    isempty(value) && error("Set $(primary) in the environment before submitting this sweep")
    return value
end

function tesla_cluster_submission_sbatch(direction::AbstractString, log_dir::AbstractString;
                                         mem::AbstractString="8G",
                                         partition::AbstractString,
                                         qos::AbstractString,
                                         time_limit::AbstractString="24:00:00")
    return Dict{Symbol, Any}(
        :job_name => "tesla_valve_bias_" * direction,
        :cpus_per_task => 16,
        :mem => mem,
        :partition => partition,
        :qos => qos,
        :time => time_limit,
        :output => joinpath(log_dir, direction * "_%j.out"),
        :error => joinpath(log_dir, direction * "_%j.err"),
    )
end

function tesla_cluster_submission_commands(; output_root::AbstractString=joinpath(@__DIR__, "..", "data",
                                                                                  "tesla_valve_cluster_bias_sweep"),
                                            mem::AbstractString=get(ENV, "TESLA_SWEEP_MEM", "8G"),
                                            dry_run::Bool=false)
    project_root = normpath(joinpath(@__DIR__, "..", "..", ".."))
    log_dir = joinpath(@__DIR__, "slurm_logs")
    mkpath(log_dir)
    partition = get(ENV, "TESLA_SWEEP_PARTITION", get(ENV, "SLURM_PARTITION", "amilan"))
    qos = get(ENV, "TESLA_SWEEP_QOS", "normal")
    time_limit = get(ENV, "TESLA_SWEEP_TIME", "24:00:00")

    commands = Dict{String, Any}()
    for direction in ("forward", "reverse")
        commands[direction] = submit_sweep!(
            script=joinpath("projects", "nonlinear", "scripts", "run_tesla_valve_cluster_bias_sweep.jl"),
            n_jobs=1,
            project_dir=project_root,
            env=Dict(
                "JULIA_NUM_THREADS" => "16",
                "TESLA_DIRECTION" => direction,
                "TESLA_OUTPUT_DIR" => output_root,
            ),
            sbatch=tesla_cluster_submission_sbatch(direction, log_dir;
                                                   mem=mem,
                                                   partition=partition,
                                                   qos=qos,
                                                   time_limit=time_limit),
            array="1-1",
            dry_run=dry_run,
        )
    end
    return commands
end

function submit_tesla_valve_cluster_bias_sweep_main()
    dry_run = get(ENV, "TESLA_SWEEP_DRY_RUN", "0") == "1"
    commands = tesla_cluster_submission_commands(; dry_run=dry_run)
    if dry_run
        for direction in ("forward", "reverse")
            println(commands[direction])
        end
    end
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    submit_tesla_valve_cluster_bias_sweep_main()
end
