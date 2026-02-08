# Lightweight helper to submit two square-bells sweeps with different wall scattering.
#
# Usage:
#   julia --project=. projects/square_bells_ucsb/scripts/submit_two_p_sweeps.jl

project_root = normpath(joinpath(@__DIR__, "..", "..", ".."))
sweep_script = joinpath(project_root, "projects", "square_bells_ucsb", "scripts", "square_bells_array.jl")

function submit_for_p(p_scatter::Float64)
    julia_bin = joinpath(Sys.BINDIR, Base.julia_exename())
    cmd = Cmd([
        julia_bin,
        "--project=$(project_root)",
        sweep_script,
        "--p_scatter=$(p_scatter)",
    ])
    @info "Submitting sweep" p_scatter cmd
    run(cmd)
end

submit_for_p(1.0)
submit_for_p(0.0)
