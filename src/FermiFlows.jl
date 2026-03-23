module FermiFlows

using FFTW
using LinearAlgebra
using SparseArrays
using StaticArrays

include("core_api.jl")
include("reference_setup.jl")

export AbstractFermiSurface2D,
       AbstractAngularDiscretization2D,
       AbstractStreamingOperator2D,
       AbstractModeRateProfile,
       AbstractCollisionModel2D,
       BandSpec,
       Isotropic2DFermiSurface,
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
       harmonic_state_nvars,
       band_momentum_weight,
       streaming_matrices,
       cosine_index,
       sine_index

end
