# Utilities for running large parameter sweeps on SLURM clusters.  provides functions for submitting a sweep job as a
# SLURM array job and for executing the assigned cases within the array job, including writing shared metadata and
# archiving the mesh for reproducibility.

using TOML
using Dates
using Sockets

solver_metadata_dict(params) = Dict(string(f) => getfield(params, f) for f in fieldnames(typeof(params)))


"""
    submit_sweep!(; 
                  script, 
                  n_jobs, 
                  project_dir=pwd(), 
                  julia_cmd="julia", 
                  script_args=String[], 
                  env=Dict(), 
                  sbatch=Dict(), 
                  array=nothing, 
                  dry_run=false)
Submit a SLURM array job to run a parameter sweep. The `script` will be executed with `julia_cmd` and
`--project=project_dir`. The `sbatch` dict can include any additional SLURM options (e.g. `:time => "01:00:00"`). 
If `array` is provided, it will override the default array specification of `1-n_jobs`. If `dry_run` is true, the 
constructed `Cmd` is returned without submission. Otherwise, the job is submitted and the function returns `nothing`.

Returns:
- `Cmd` when `dry_run=true`.
- `nothing` when submitted.
"""
function submit_sweep!(;
                       script::AbstractString,
                       n_jobs::Integer,
                       project_dir::AbstractString=pwd(),
                       julia_cmd::AbstractString="julia",
                       script_args::Vector{String}=String[],
                       env::Dict{String, String}=Dict{String, String}(),
                       sbatch::Dict{Symbol, Any}=Dict{Symbol, Any}(),
                       array::Union{Nothing, AbstractString}=nothing,
                       dry_run::Bool=false)
    script_path = isabspath(script) ? script : joinpath(project_dir, script)
    isfile(script_path) || error("Script not found: $script_path")

    has_cpus_per_task = haskey(sbatch, :cpus_per_task) &&
                        sbatch[:cpus_per_task] !== nothing &&
                        sbatch[:cpus_per_task] !== false
    has_tres_per_task = haskey(sbatch, :tres_per_task) &&
                        sbatch[:tres_per_task] !== nothing &&
                        sbatch[:tres_per_task] !== false
    if has_cpus_per_task && has_tres_per_task
        error("Specify only one of :cpus_per_task or :tres_per_task in sbatch.")
    end

    array_spec = array === nothing ? string(get(sbatch, :array, "1-$(n_jobs)")) : string(array)
    sbatch_args = String["--array=$(array_spec)"]
    export_val = get(sbatch, :export, nothing)

    if !isempty(env)
        env_pairs = join(["$k=$v" for (k, v) in env], ",")
        base_export = export_val === nothing ? "ALL" : string(export_val)
        export_val = isempty(env_pairs) ? base_export : "$(base_export),$(env_pairs)"
    end

    for (key, value) in sbatch
        key == :array && continue
        key == :export && continue
        flag = "--" * replace(string(key), "_" => "-")
        if value === true
            push!(sbatch_args, flag)
        elseif value === false || value === nothing
            continue
        else
            push!(sbatch_args, "$(flag)=$(value)")
        end
    end
    if export_val !== nothing
        push!(sbatch_args, "--export=$(export_val)")
    end

    wrap_parts = [julia_cmd, "--project=$(project_dir)", script_path]
    append!(wrap_parts, script_args)
    wrap_cmd = join(Base.shell_escape.(wrap_parts), " ")

    cmd = Cmd(["sbatch"; sbatch_args; "--wrap"; wrap_cmd])
    run_env = copy(ENV)
    if has_cpus_per_task
        delete!(run_env, "SBATCH_TRES_PER_TASK")
    elseif has_tres_per_task
        delete!(run_env, "SBATCH_CPUS_PER_TASK")
    end
    @info "Submitting SLURM sweep" script=script_path array=array_spec
    dry_run && return cmd
    run(setenv(cmd, run_env))
    return nothing
end


"""
    write_sweep_metadata!(data_dir::AbstractString, 
                          params, sweep_metadata, slurm_metadata,
                          gamma_mr_vals, gamma_mc_vals)
    
Write sweep metadata to `sweep_metadata.toml`.

Returns:
- `metadata_path::String`: path to written TOML file.
"""
function write_sweep_metadata!(data_dir::AbstractString, 
                               params, sweep_metadata, slurm_metadata, 
                               gamma_mr_vals, gamma_mc_vals)
        # put all metadata together
        metadata = Dict(
        "solver_metadata" => solver_metadata_dict(params),
        "sweep_metadata" => merge(
            Dict(
                "gamma_mr_values" => collect(gamma_mr_vals),
                "gamma_mc_values" => collect(gamma_mc_vals),
            ),
            Dict(pairs(sweep_metadata)),
        ),
        "slurm_metadata" => Dict(pairs(slurm_metadata)),
    )
    metadata_path = joinpath(data_dir, "sweep_metadata.toml")
    mkpath(data_dir)
    open(metadata_path, "w") do io
        TOML.print(io, metadata)
    end
    @info "Wrote sweep metadata" path=metadata_path
    return metadata_path
end


"""
    archive_mesh!(mesh_path::AbstractString, data_dir::AbstractString)

Copy the mesh file to the sweep data directory for reproducibility.

Returns:
- `archive_path::String`: destination path of archived mesh.
"""
function archive_mesh!(mesh_path::AbstractString, data_dir::AbstractString)
    archive_path = joinpath(data_dir, basename(mesh_path))
    @info "Archiving mesh file for reproducibility..."
    cp(mesh_path, archive_path; force=true)
    @info "Archived mesh" source=mesh_path dest=archive_path
    return archive_path
end


"""
    copy_mesh_to_scratch(mesh_path::AbstractString, prefix::AbstractString)

Copy the mesh file to a scratch directory on the compute node and return the local path. This is done so that we 
avoid data races when many jobs are trying to load the meshes a bunch of times.

Returns:
- `local_path::String`: copied mesh path in scratch space.
"""
function copy_mesh_to_scratch(mesh_path::AbstractString, prefix::AbstractString)
    scratch_root = get(ENV, "SLURM_TMPDIR", get(ENV, "TMPDIR", tempdir()))
    scratch = mktempdir(scratch_root; prefix=prefix)
    local_path = joinpath(scratch, basename(mesh_path))
    cp(mesh_path, local_path)
    @info "Mesh copied to scratch" source=mesh_path scratch=local_path
    return local_path
end


"""
    select_cases(gamma_mr_vals, gamma_mc_vals; cases_per_task::Int, env::AbstractDict=ENV)

Determine which cases to run for the current task based on the total parameter grid and the SLURM_ARRAY_TASK_ID. 

Returns:
- `(case_indices, total_cases, n_cases_this_task, task_id)`.
"""
function select_cases(gamma_mr_vals, gamma_mc_vals; cases_per_task::Int, env::AbstractDict=ENV)
    task_id = parse(Int, get(env, "SLURM_ARRAY_TASK_ID", "1"))
    total_cases = length(gamma_mr_vals) * length(gamma_mc_vals)
    base_index = (task_id - 1) * cases_per_task
    last_index = min(base_index + cases_per_task, total_cases)
    case_indices = (base_index + 1):last_index
    return case_indices, total_cases, length(case_indices), task_id
end


"""
    grid_lookup(gamma_mr_vals, gamma_mc_vals, index_global::Int)

Given a global case index, look up the corresponding gamma_mr and gamma_mc values from the parameter grid. Assumes the 
grid is ordered with gamma_mc varying fastest.

Returns:
- `(gamma_mr, gamma_mc)` for `index_global`.
"""
function grid_lookup(gamma_mr_vals, gamma_mc_vals, index_global::Int)
    n_mc = length(gamma_mc_vals)
    index_mr = div(index_global - 1, n_mc) + 1
    index_mc = mod(index_global - 1, n_mc) + 1
    return gamma_mr_vals[index_mr], gamma_mc_vals[index_mc]
end


"""
    ordered_case_indices(gamma_mr_vals, gamma_mc_vals)

Return global case indices ordered from hardest to easiest expected solves,
using descending `(gamma_mr + gamma_mc, gamma_mr, gamma_mc)`.
"""
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
