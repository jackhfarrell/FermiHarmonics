module ElectronKinetics

using FFTW
using LinearAlgebra
using SparseArrays
using StaticArrays

include("core_api.jl")
include("live_visualization_api.jl")
include("reference_setup.jl")
include("mesh_generation.jl")
include("slurm_utils.jl")

const SolveParams = SolverConfig

export AbstractFermiSurface2D,
       AbstractAngularDiscretization2D,
       AbstractStreamingOperator2D,
       AbstractModeRateProfile,
       AbstractCollisionModel2D,
       Band,
       Isotropic2DFermiSurface,
       EllipticFermiSurface2D,
       GeneralFermiSurface2D,
       HarmonicBasis,
       AngleGrid,
       IsotropicHarmonicStreaming,
       IsotropicAngleStreaming,
       TwoRateProfile,
       OddQuarticRateProfile,
       ConstantModeRateProfile,
       CustomModeRateProfile,
       LinearBGKCollision,
       QuadraticBGKCollision,
       ExactAngleBGKCollision,
       TwoRateAngleBGKCollision,
       KineticModel2D,
       SolverConfig,
       SolveParams,
       LiveVisualizationConfig,
       LiveProgressSnapshot,
       LiveFieldSnapshot,
       LiveVisualizationSnapshot,
       MeshBuildConfig,
       TrixiProblem,
       MaxwellWallBC,
       OhmicContactBC,
       blg_reference_setup,
       estimate_max_harmonic,
       save_solution_custom,
       save_for_analysis,
       save_mesh_native_analysis,
       evaluate_solution,
       evaluate_observables,
       enable_nonlinear_timing!,
       disable_nonlinear_timing!,
       reset_nonlinear_timing!,
       nonlinear_timing_snapshot,
       print_nonlinear_timing_summary,
       solve_status,
       solve,
       mode_rate,
       collision_sources!,
       generate_mesh_from_geo,
       harmonic_state_nvars,
       band_momentum_weight,
       surface_vF,
       surface_max_speed,
       surface_vF_angle,
       surface_density_of_states,
       surface_mass,
       surface_charge,
       streaming_matrices,
       cosine_index,
       sine_index,
       residual_progress_fraction,
       mesh_provenance_attributes,
       resolve_mesh_path,
       submit_sweep!,
       write_sweep_metadata!,
       archive_mesh!,
       copy_mesh_to_scratch,
       select_cases,
       grid_lookup,
       ordered_case_indices

end
