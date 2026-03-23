# simple script to run a demo with four regimes:
# 1. ohmic/diffusive (gamma_mr=100, gamma_ee=0)
# 2. hydrodynamic (gamma_mr=0, gamma_ee=100)
# 3. intermediate/mixed, hydro-leaning (gamma_mr=0.05, gamma_ee=60)
# 4. ballistic (gamma_mr=0, gamma_ee=0)
# runs are warm-started in sequence. harmonics are chosen automatically with
# min_harmonic=4 and max_harmonic_auto=100.

using Plots
using ElectronKinetics
using Trixi

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

    config = SolverConfig(;
        min_harmonic = 4,
        max_harmonic_auto = 100,
        polydeg = 3,
        tspan_end = 150.0,
        residual_tol = 1e-4,
        cfl = 0.8,
        log_every = 1000,
    )

    # Keep the "intermediate" demo closer to hydrodynamic by default (gamma_ee >> gamma_mr).
    # Optional overrides:
    #   FERMI_DEMO_INTERMEDIATE_GAMMA_MR
    #   FERMI_DEMO_INTERMEDIATE_GAMMA_EE
    intermediate_gamma_mr = parse(Float64, get(ENV, "FERMI_DEMO_INTERMEDIATE_GAMMA_MR", "0.05"))
    intermediate_gamma_ee = parse(Float64, get(ENV, "FERMI_DEMO_INTERMEDIATE_GAMMA_EE", "60.0"))

    regimes = [
        (name = "diffusive", gamma_mr = 100.0, gamma_ee = 0.0),
        (name = "hydrodynamic", gamma_mr = 0.0, gamma_ee = 100.0),
        (name = "intermediate", gamma_mr = intermediate_gamma_mr, gamma_ee = intermediate_gamma_ee),
        (name = "ballistic", gamma_mr = 0.0, gamma_ee = 0.0),
    ]
    requested_regime = get(ENV, "FERMI_DEMO_REGIME", "")
    if !isempty(requested_regime)
        regimes = filter(r -> r.name == requested_regime, regimes)
        isempty(regimes) && error("Unknown FERMI_DEMO_REGIME=$(requested_regime). Expected one of: diffusive, hydrodynamic, intermediate, ballistic")
    end

# ======================================================================================================================
# Warm-Started Sweep
# ======================================================================================================================

    u0 = nothing

    for regime in regimes
        @info "Running regime" name=regime.name gamma_mr=regime.gamma_mr gamma_ee=regime.gamma_ee

        model = KineticModel2D(
            Isotropic2DFermiSurface(; vF=1.0, nu=1.0, mass=1.0, charge=-1.0),
            HarmonicBasis(:auto),
            IsotropicHarmonicStreaming(),
            LinearBGKCollision(regime.gamma_mr, TwoRateProfile(regime.gamma_ee)),
        )

        sol, semi = solve(
            TrixiProblem(; mesh_path=mesh_path, boundary_conditions=boundary_conditions),
            model,
            config;
            u0_override=u0,
            visualize=true,
            name=regime.name,
        )

        save_path = joinpath(output_dir, "$(regime.name).h5")
        save_for_analysis(sol, semi, save_path)
        @info "Saved analysis output" path=save_path

        u0 = copy(sol.u[end])
    end

    @info "Demo complete" output_dir=output_dir
    return nothing
end

main()
