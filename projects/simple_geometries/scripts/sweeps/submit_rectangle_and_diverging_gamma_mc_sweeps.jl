# submit gamma_mc-only sweeps for rectangle, diverging nozzle, and dogleg with fixed p_scatter and gamma_mr.

project_root = normpath(joinpath(@__DIR__, "..", ".."))
main_project = normpath(joinpath(project_root, "..", ".."))

p_scatter = 1.0
gamma_mr = 1e-2

common_env = Dict(
    "FERMI_P_SCATTER" => string(p_scatter),
    "FERMI_GAMMA_MR" => string(gamma_mr),
)

sweep_scripts = [
    joinpath(@__DIR__, "rectangle_gamma_mc_sweep_array.jl"),
    joinpath(@__DIR__, "diverging_nozzle_gamma_mc_sweep_array.jl"),
    joinpath(@__DIR__, "dogleg_gamma_mc_sweep_array.jl"),
]

for script in sweep_scripts
    println("Submitting $(basename(script)) with FERMI_P_SCATTER=$(p_scatter), FERMI_GAMMA_MR=$(gamma_mr)")
    cmd = `julia --project=$(main_project) $(script)`
    run(setenv(cmd, merge(copy(ENV), common_env)))
end
