module FermiFlowsTrixiExt

using FermiFlows
using SciMLBase
using Trixi
using LinearAlgebra
using SparseArrays
using StaticArrays
using HDF5
using FFTW
using NLsolve

const FermiHarmonics = @__MODULE__

import FermiFlows: solve,
                   solve_status,
                   save_solution_custom,
                   save_for_analysis,
                   save_mesh_native_analysis,
                   evaluate_solution,
                   evaluate_observables,
                   create_live_dashboard,
                   finalize_live_dashboard!,
                   update_live_dashboard!,
                   collision_sources!,
                   BandSpec,
                   BCProjectorCache,
                   CustomModeRateProfile,
                   ExactAngleBGKCollision,
                   HarmonicBasis,
                   KineticModel2D,
                   LiveFieldSnapshot,
                   LiveProgressSnapshot,
                   LiveVisualizationConfig,
                   LiveVisualizationSnapshot,
                   LinearBGKCollision,
                   MaxwellWallBC,
                   NonlinearBoundaryFaceData,
                   OddQuarticRateProfile,
                   OhmicContactBC,
                   QuadraticBGKCollision,
                   SolverConfig,
                   TrixiProblem,
                   TwoRateAngleBGKCollision,
                   TwoRateProfile,
                   band_momentum_weight,
                   boundary_condition_name,
                   coerce_band_spec,
                   collision_electrostatic_coupling,
                   collision_gamma_mc,
                   collision_gamma_mr,
                   collision_mass,
                   collision_mu0,
                   collision_symbol,
                   collision_theta_oversample,
                   cosine_index,
                   create_angle_transport_data,
                   create_nonlinear_transport_data,
                   estimate_max_harmonic,
                   harmonics_flux!,
                   harmonic_state_nvars,
                   mode_profile,
                   mode_rate,
                   nonlinear_timestep_speed,
                   nonlinear_timing_enabled,
                   profile_reference_rate,
                   record_nonlinear_timing!,
                   residual_progress_fraction,
                   resolve_max_harmonic,
                   resize_multiband_warm_start,
                   resize_warm_start,
                   surface_vF,
                   streaming_matrices,
                   transport_symbol,
                   validate,
                   validate_nonlinear_warm_start,
                   zero_state_speed

include("../src/trixi_equations_internal.jl")
include("../src/nonlinear_transport.jl")
include("../src/source_terms.jl")
include("../src/boundary_conditions.jl")
include("../src/trixi_interface.jl")
include("../src/io_utils.jl")
include("../src/live_visualization_trixi.jl")
include("../src/trixi_runner.jl")

end
