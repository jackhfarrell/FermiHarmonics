# Helper to submit all simple_geometries sweep scripts configured for 125 jobs each
# (total expected submissions: 625 SLURM array jobs).
#
# Usage:
#   julia --project=. projects/simple_geometries/scripts/sweeps/submit_all_sweeps.jl

project_root = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))
simple_geometries_root = joinpath(project_root, "projects", "simple_geometries")

julia_bin = joinpath(Sys.BINDIR, Base.julia_exename())

sweeps = [
    (name = "rectangle", path = joinpath(simple_geometries_root, "scripts", "sweeps", "rectangle_array.jl"), jobs = 125),
    (name = "diverging_nozzle", path = joinpath(simple_geometries_root, "scripts", "sweeps", "diverging_nozzle_array.jl"), jobs = 125),
    (name = "junction", path = joinpath(simple_geometries_root, "scripts", "junction_array.jl"), jobs = 125),
    (name = "intersection", path = joinpath(simple_geometries_root, "scripts", "intersection_config_1.jl"), jobs = 125),
    (name = "chamber", path = joinpath(simple_geometries_root, "scripts", "chamber_array.jl"), jobs = 125),
]

total_jobs = sum(s.jobs for s in sweeps)
@assert total_jobs == 625 "Expected 625 total jobs, got $(total_jobs)"

for sweep in sweeps
    isfile(sweep.path) || error("Sweep script not found: $(sweep.path)")
    cmd = Cmd([
        julia_bin,
        "--project=$(project_root)",
        sweep.path,
    ])
    @info "Submitting sweep" name=sweep.name expected_jobs=sweep.jobs cmd
    run(cmd)
end

@info "All sweep submissions complete" total_expected_jobs=total_jobs
