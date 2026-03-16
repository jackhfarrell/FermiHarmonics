# submit the 2D dogleg gamma_mr x gamma_mc sweep with an optional wall-scattering override.

project_root = normpath(joinpath(@__DIR__, "..", ".."))
main_project = normpath(joinpath(project_root, "..", ".."))

p_scatter = 1.0

common_env = Dict(
    "FERMI_P_SCATTER" => string(p_scatter),
)

script = joinpath(@__DIR__, "dogleg_gamma_grid_sweep_array.jl")
println("Submitting $(basename(script)) with FERMI_P_SCATTER=$(p_scatter)")
cmd = `julia --project=$(main_project) $(script)`
run(setenv(cmd, merge(copy(ENV), common_env)))
