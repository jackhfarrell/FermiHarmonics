module ElectronKineticsTrixiExt

using ElectronKinetics
using SciMLBase
using Trixi
using LinearAlgebra
using SparseArrays
using StaticArrays
using HDF5
using FFTW
using NLsolve

const ElectronKineticsExt = @__MODULE__

import ElectronKinetics: AbstractAnalyticSurface,
                         AbstractBoundaryCondition,
                         AbstractContactBC,
                         AbstractGridDiscretization,
                         AbstractHarmonicDiscretization,
                         AbstractLinearCollision,
                         AbstractNonlinearAngleCollision,
                         AbstractUserDefinedSurface,
                         AbstractWallBC,
                         LinearTransport,
                         NonlinearParabolicTransport,
                         solve,
                         solve_status,
                         save_solution_custom,
                         save_for_analysis,
                         save_mesh_native_analysis,
                         compute_analysis_grids,
                         current_norm_variables,
                         visualization_callback,
                         evaluate_solution,
                         evaluate_observables,
                         create_live_dashboard,
                         finalize_live_dashboard!,
                         update_live_dashboard!,
                         collision_sources!,
                         Band,
                         BCProjectorCache,
                         create_angle_transport_data,
                         create_nonlinear_transport_data,
                         CustomModeRateProfile,
                         ExactAngleBGKCollision,
                         HarmonicBasis,
                         KineticModel2D,
                         LiveFieldSnapshot,
                         LiveProgressSnapshot,
                         LiveVisualizationConfig,
                         LiveVisualizationSnapshot,
                         live_dashboard_is_open,
                         LinearBGKCollision,
                         LinearCollisionMatrix,
                         MeshBuildConfig,
                         MaxwellWallBC,
                         NonlinearBoundaryFaceData,
                         OddQuarticRateProfile,
                         OhmicContactBC,
                         QuadraticBGKCollision,
                         SolveParams,
                         SolverConfig,
                         TrixiProblem,
                         TwoRateAngleBGKCollision,
                         AngleRateBGKCollision,
                         TwoRateProfile,
                         AngleGrid,
                         Isotropic2DFermiSurface,
                         IsotropicAngleStreaming,
                         IsotropicHarmonicStreaming,
                         band_momentum_weight,
                         boundary_condition_name,
                         coerce_band,
                         collision_electrostatic_coupling,
                         collision_gamma_ee,
                         collision_gamma_mr,
                         collision_mass,
                         collision_mu0,
                         collision_symbol,
                         collision_theta_oversample,
                         cosine_index,
                         estimate_max_harmonic,
                         harmonics_flux!,
                         harmonic_state_nvars,
                         mode_profile,
                         mode_rate,
                         build_collision_matrix,
                         nonlinear_timestep_speed,
                         nonlinear_timing_enabled,
                         profile_reference_rate,
                         record_nonlinear_timing!,
                         residual_progress_fraction,
                         resolve_max_harmonic,
                         resize_multiband_warm_start,
                         resize_warm_start,
                         surface_vF,
                         surface_max_speed,
                         surface_vF_angle,
                         surface_density_of_states,
                         surface_mass,
                         surface_charge,
                         streaming_matrices,
                         transport_symbol,
                         validate,
                         validate_nonlinear_warm_start,
                         zero_state_speed

include("trixi/trixi_equations_internal.jl")
include("trixi/nonlinear_transport.jl")
include("trixi/source_terms.jl")
include("trixi/boundary_conditions.jl")
include("trixi/trixi_interface.jl")
include("trixi/io_utils.jl")
include("trixi/live_visualization_trixi.jl")
include("trixi/trixi_runner.jl")

function legacy_collision_model(
    gamma_mr::Real,
    gamma_ee::Real;
    transport::Symbol,
    collision_model::Union{Nothing, Symbol},
    gamma3::Union{Nothing, Real},
    mu0::Union{Nothing, Real},
    mass::Union{Nothing, Real},
    chi::Real,
)
    if transport === :linear
        isnothing(collision_model) || collision_model in (:linear, :linear_mrt) ||
            throw(ArgumentError("transport=:linear supports only collision_model=:linear_mrt"))
        return LinearBGKCollision(gamma_mr, TwoRateProfile(gamma_ee))
    end

    transport === :parabolic_nonlinear ||
        throw(ArgumentError("unsupported transport=$(transport)"))

    collision_symbol_value = isnothing(collision_model) ? :quadratic_bgk : collision_model
    if collision_symbol_value === :quadratic_bgk
        isnothing(mu0) && throw(ArgumentError("mu0 is required for collision_model=:quadratic_bgk"))
        isnothing(mass) && throw(ArgumentError("mass is required for collision_model=:quadratic_bgk"))
        profile = isnothing(gamma3) ? OddQuarticRateProfile(gamma_ee) : OddQuarticRateProfile(gamma_ee, gamma3)
        return QuadraticBGKCollision(
            gamma_mr,
            profile;
            mu0=mu0,
            mass=mass,
            electrostatic_coupling=chi,
        )
    elseif collision_symbol_value === :exact_bgk
        isnothing(mu0) && throw(ArgumentError("mu0 is required for collision_model=:exact_bgk"))
        isnothing(mass) && throw(ArgumentError("mass is required for collision_model=:exact_bgk"))
        return ExactAngleBGKCollision(;
            gamma_mr=gamma_mr,
            gamma_ee=gamma_ee,
            mu0=mu0,
            mass=mass,
            electrostatic_coupling=chi,
        )
    elseif collision_symbol_value === :two_rate_bgk
        isnothing(mu0) && throw(ArgumentError("mu0 is required for collision_model=:two_rate_bgk"))
        isnothing(mass) && throw(ArgumentError("mass is required for collision_model=:two_rate_bgk"))
        return TwoRateAngleBGKCollision(;
            gamma_mr=gamma_mr,
            gamma_ee=gamma_ee,
            mu0=mu0,
            mass=mass,
            electrostatic_coupling=chi,
        )
    elseif collision_symbol_value === :angle_rate_bgk
        isnothing(mu0) && throw(ArgumentError("mu0 is required for collision_model=:angle_rate_bgk"))
        isnothing(mass) && throw(ArgumentError("mass is required for collision_model=:angle_rate_bgk"))
        profile = isnothing(gamma3) ? OddQuarticRateProfile(gamma_ee) : OddQuarticRateProfile(gamma_ee, gamma3)
        return AngleRateBGKCollision(
            gamma_mr,
            profile;
            mu0=mu0,
            mass=mass,
            electrostatic_coupling=chi,
        )
    end

    throw(ArgumentError("unsupported collision_model=$(collision_symbol_value)"))
end

function legacy_live_visualization(
    config::SolverConfig;
    visualize::Bool,
    visualize_every::Union{Nothing, Integer},
    visualization_mode::Symbol,
    viz_field::Union{Nothing, Symbol}=nothing,
    viz_colormap::Symbol=:magma,
)
    visualize || return nothing
    interval = isnothing(visualize_every) ? config.log_every : Int(visualize_every)
    interval > 0 || throw(ArgumentError("visualize_every must be positive"))
    return LiveVisualizationConfig(;
        accepted_step_interval=interval,
        geometry_mode=visualization_mode,
        field=viz_field,
        colormap=viz_colormap,
        show_window=true,
    )
end

function legacy_surface(mu0::Union{Nothing, Real}, mass::Union{Nothing, Real})
    if isnothing(mu0) || isnothing(mass)
        return Isotropic2DFermiSurface()
    end
    return Isotropic2DFermiSurface(; vF=sqrt(2 * Float64(mu0) / Float64(mass)), nu=1.0, mass=mass, charge=-1.0)
end

function solve(
    mesh_path::AbstractString,
    boundary_conditions::Dict{Symbol, Any},
    params::SolveParams,
    gamma_mr::Real,
    gamma_ee::Real;
    max_harmonic::Union{Integer, Symbol, Nothing}=:auto,
    n_angles::Union{Nothing, Integer}=nothing,
    transport::Symbol=:linear,
    collision_model::Union{Nothing, Symbol}=nothing,
    gamma3::Union{Nothing, Real}=nothing,
    mu0::Union{Nothing, Real}=nothing,
    mass::Union{Nothing, Real}=nothing,
    chi::Real=0.0,
    u0_override::Union{Nothing, AbstractVector}=nothing,
    visualize::Bool=false,
    visualize_every::Union{Nothing, Integer}=nothing,
    visualization_mode::Symbol=:cartesian,
    viz_field::Union{Nothing, Symbol}=nothing,
    viz_colormap::Symbol=:magma,
    mesh_build::MeshBuildConfig=MeshBuildConfig(),
    name::AbstractString="run",
)
    validate(params)
    collision = legacy_collision_model(
        gamma_mr,
        gamma_ee;
        transport=transport,
        collision_model=collision_model,
        gamma3=gamma3,
        mu0=mu0,
        mass=mass,
        chi=chi,
    )
    discretization = if collision isa Union{ExactAngleBGKCollision, TwoRateAngleBGKCollision, AngleRateBGKCollision}
        isnothing(n_angles) && throw(ArgumentError("n_angles is required for collision_model=$(collision_symbol(collision))"))
        AngleGrid(n_angles)
    else
        !isnothing(n_angles) && throw(ArgumentError("n_angles is only supported for angle-based nonlinear collision models"))
        HarmonicBasis(max_harmonic)
    end

    model = KineticModel2D(
        legacy_surface(mu0, mass),
        discretization,
        discretization isa AngleGrid ? IsotropicAngleStreaming() : IsotropicHarmonicStreaming(),
        collision,
    )

    return solve(
        TrixiProblem(; mesh_path=String(mesh_path), boundary_conditions=boundary_conditions, mesh_build=mesh_build),
        model,
        params;
        u0_override=u0_override,
        live_visualization=legacy_live_visualization(params; visualize=visualize, visualize_every=visualize_every, visualization_mode=visualization_mode, viz_field=viz_field, viz_colormap=viz_colormap),
        name=name,
    )
end

function solve(
    mesh_path::AbstractString,
    boundary_conditions::Dict{Symbol, Any},
    params::SolveParams,
    bands::AbstractVector;
    max_harmonic::Union{Integer, Symbol, Nothing}=:auto,
    transport::Symbol=:linear,
    collision_model::Union{Nothing, Symbol}=nothing,
    gamma_drag::Real=0.0,
    gamma3::Union{Nothing, Real}=nothing,
    mu0::Union{Nothing, Real}=nothing,
    mass::Union{Nothing, Real}=nothing,
    chi::Real=0.0,
    u0_override::Union{Nothing, AbstractVector}=nothing,
    visualize::Bool=false,
    visualize_every::Union{Nothing, Integer}=nothing,
    visualization_mode::Symbol=:cartesian,
    mesh_build::MeshBuildConfig=MeshBuildConfig(),
    name::AbstractString="run",
)
    validate(params)
    transport === :linear || throw(ArgumentError("legacy multiband solve currently supports only transport=:linear"))
    isnothing(collision_model) || collision_model in (:linear, :linear_mrt) ||
        throw(ArgumentError("legacy multiband solve supports only collision_model=:linear_mrt"))
    isnothing(gamma3) || throw(ArgumentError("gamma3 is not supported for multiband linear solves"))
    isnothing(mu0) || throw(ArgumentError("mu0 is not supported for multiband linear solves"))
    isnothing(mass) || throw(ArgumentError("mass is not supported for multiband linear solves"))
    chi == 0.0 || throw(ArgumentError("chi is not supported for multiband linear solves"))

    band_specs = Band[coerce_band(band) for band in bands]
    isempty(band_specs) && throw(ArgumentError("multiband solve requires at least one band"))
    first_band = band_specs[1]
    model = KineticModel2D(
        first_band.surface,
        HarmonicBasis(max_harmonic),
        IsotropicHarmonicStreaming(),
        LinearBGKCollision(0.0, TwoRateProfile(0.0));
        bands=band_specs,
        gamma_drag=gamma_drag,
    )

    return solve(
        TrixiProblem(; mesh_path=String(mesh_path), boundary_conditions=boundary_conditions, mesh_build=mesh_build),
        model,
        params;
        u0_override=u0_override,
        live_visualization=legacy_live_visualization(params; visualize=visualize, visualize_every=visualize_every, visualization_mode=visualization_mode),
        name=name,
    )
end

@inline function current_norm_variables(u, equations::FermiHarmonics2D{NVARS}) where {NVARS}
    a1 = length(u) >= 2 ? u[2] : 0.0
    b1 = length(u) >= 3 ? u[3] : 0.0
    j_norm = hypot(a1, b1)
    return SVector{NVARS, Float64}(ntuple(i -> i == 1 ? j_norm : 0.0, NVARS))
end

function Trixi.varnames(::typeof(current_norm_variables), equations::FermiHarmonics2D{NVARS}) where {NVARS}
    return ntuple(i -> i == 1 ? "j_norm" : "_viz_pad_$(i)", NVARS)
end

function visualization_callback(params::SolveParams, semi, name::AbstractString; interval::Int=params.log_every, mode::Symbol=:cartesian)
    if transport_is_nonlinear(semi.equations)
        analysis_path = mode === :mesh_native ? "live_viz_$(name)_mesh_native.h5" : "live_viz_$(name).h5"
        nvisnodes = 120

        return SciMLBase.DiscreteCallback(
            (u, t, integrator) -> integrator.stats.naccept % interval == 0,
            integrator -> begin
                if mode === :mesh_native
                    mesh_data = compute_mesh_native_analysis(integrator.u, semi; refine=6)
                    analysis_write_mesh_native_hdf5(analysis_path, mesh_data, integrator.t, semi.equations)
                elseif mode === :cartesian
                    grids = compute_analysis_grids(integrator.u, semi; nvisnodes=nvisnodes)
                    analysis_write_hdf5(
                        analysis_path,
                        grids.density,
                        grids.a1,
                        grids.b1,
                        grids.jx,
                        grids.jy,
                        grids.x,
                        grids.y,
                        grids.mask,
                        integrator.t,
                        grids.equations,
                        band_grids=grids.bands,
                    )
                else
                    throw(ArgumentError("unsupported visualization_mode=$(mode); use :cartesian or :mesh_native"))
                end
                @info "Updated nonlinear analysis snapshot" path=analysis_path t=round(integrator.t, digits=4)
                nothing
            end;
            save_positions=(false, false),
        )
    end

    variable_names = semi.equations isa MultiBandFermiHarmonics2D ?
        collect(Trixi.varnames(Trixi.cons2cons, semi.equations)) :
        ["a0", "a1", "b1"]
    return Trixi.VisualizationCallback(
        semi;
        interval=interval,
        variable_names=variable_names,
        filename="live_viz_$(name)",
        overwrite=true,
        seriescolor=:magma,
    )
end

end
