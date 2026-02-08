# simple script to run a demo with three regimes:
# 1. ohmic/diffusive (gamma_mr=100, gamma_mc=0)
# 2. hydrodynamic (gamma_mr=0, gamma_mc=100)
# 3. ballistic (gamma_mr=0, gamma_mc=0)
# runs are warm-started in sequence. harmonics are chosen automatically with
# min_harmonic=4 and max_harmonic_auto=100.

using Plots
using FermiHarmonics

# ======================================================================================================================
# Configuration
# ======================================================================================================================

function main()
    project_root = normpath(joinpath(@__DIR__, ".."))
    mesh_path = joinpath(project_root, "projects", "square_bells_ucsb", "mesh", "square_bells.inp")
    output_dir = joinpath(@__DIR__, "data")
    mkpath(output_dir)

    bias = 1.0
    p_scatter = 1.0

    boundary_conditions = Dict(
        :walls => MaxwellWallBC(p_scatter),
        :contact_top => OhmicContactBC(-bias / 2),
        :contact_bottom => OhmicContactBC(bias / 2),
    )

    params = SolveParams(;
        min_harmonic = 4,
        max_harmonic_auto = 100,
        polydeg = 3,
        tspan_end = 150.0,
        residual_tol = 1e-4,
        cfl = 0.8,
        log_every = 1000,
    )

    regimes = [
        (name = "diffusive", gamma_mr = 100.0, gamma_mc = 0.0),
        (name = "hydrodynamic", gamma_mr = 0.0, gamma_mc = 100.0),
        (name = "ballistic", gamma_mr = 0.0, gamma_mc = 0.0),
    ]

# ======================================================================================================================
# Warm-Started Sweep
# ======================================================================================================================

    u0 = nothing

    for regime in regimes
        @info "Running regime" name=regime.name gamma_mr=regime.gamma_mr gamma_mc=regime.gamma_mc

        sol, semi = FermiHarmonics.solve(
            mesh_path,
            boundary_conditions,
            params,
            regime.gamma_mr,
            regime.gamma_mc;
            max_harmonic=:auto,
            u0_override=u0,
            visualize=true,
            name=regime.name,
        )

        save_path = joinpath(output_dir, "$(regime.name).h5")
        FermiHarmonics.save_for_analysis(sol, semi, save_path)
        @info "Saved analysis output" path=save_path

        u0 = copy(sol.u[end])
    end

    @info "Demo complete" output_dir=output_dir
    return nothing
end

main()
