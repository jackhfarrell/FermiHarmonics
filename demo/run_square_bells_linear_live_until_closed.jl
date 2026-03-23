using Dates
using Pkg

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(PROJECT_ROOT; io=devnull)

if Base.find_package("GLMakie") === nothing
    error("This live demo requires GLMakie in the active environment. Run `julia --project=. -e 'using Pkg; Pkg.add(\"GLMakie\")'` once, then rerun this script.")
end

using ElectronKinetics
import GLMakie
using SciMLBase
using Trixi

function main()
    mesh_path = joinpath(PROJECT_ROOT, "projects", "square_bells_ucsb", "mesh", "square_bells.inp")
    output_dir = joinpath(@__DIR__, "data_live")
    mkpath(output_dir)

    reference = blg_reference_setup()
    bias = 1.0
    gamma_mr = 0.0
    gamma_mc = 100.0
    run_stamp = Dates.format(now(), "yyyy-mm-dd_HHMMSS")
    run_name = "square_bells_linear_live_gamma_mc100_Mauto_" * run_stamp

    boundary_conditions = Dict(
        :walls => MaxwellWallBC(reference.p_scatter),
        :contact_top => OhmicContactBC(-bias / 2),
        :contact_bottom => OhmicContactBC(bias / 2),
    )

    model = KineticModel2D(
        Isotropic2DFermiSurface(; vF=reference.vF, nu=1.0, mass=reference.mass, charge=-1.0),
        HarmonicBasis(:auto),
        IsotropicHarmonicStreaming(),
        LinearBGKCollision(gamma_mr, TwoRateProfile(gamma_mc)),
        reference=reference,
    )

    config = SolverConfig(;
        polydeg=3,
        tspan_end=1.0e6,
        residual_tol=1.0e-14,
        cfl=0.8,
        log_every=200,
        min_harmonic=4,
        max_harmonic_auto=100,
    )

    live_visualization = LiveVisualizationConfig(;
        geometry_mode=:cartesian,
        field=:current_magnitude,
        nvisnodes=200,
        accepted_step_interval=1000,
        min_update_seconds=0.2,
        show_window=true,
    )

    @info "Running square-bells linear live demo" mesh_path gamma_mr gamma_mc harmonic_mode=:auto bias run_name
    @info "Close the live window to stop the solve and save the current state" output_dir

    sol, semi = ElectronKinetics.solve(
        ElectronKinetics.TrixiProblem(; mesh_path=mesh_path, boundary_conditions=boundary_conditions),
        model,
        config;
        live_visualization=live_visualization,
        name=run_name,
    )

    restart_path = joinpath(output_dir, run_name * "_restart.h5")
    analysis_path = joinpath(output_dir, run_name * "_mesh_native.h5")
    ElectronKinetics.save_solution_custom(sol, semi, restart_path)
    ElectronKinetics.save_mesh_native_analysis(sol, semi, analysis_path; refine=4)

    stop_reason_override = hasproperty(sol, :retcode) && getproperty(sol, :retcode) === SciMLBase.ReturnCode.Terminated ? :window_closed : nothing
    status = ElectronKinetics.solve_status(sol, semi, config; stop_reason_override=stop_reason_override)
    @info "Saved square-bells live demo output" stop_reason=status.stop_reason final_time=status.final_time restart_path analysis_path
    return nothing
end

main()
