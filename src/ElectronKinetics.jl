"""
    ElectronKinetics

Extensible 2D Fermi-liquid kinetic transport framework.

# Overview

ElectronKinetics provides composable components for space-resolved kinetic transport
simulations on unstructured 2D meshes. It separates physics components from computational
backends through Julia's type system and multiple dispatch.

# Core Architecture

## Type Hierarchies

### Fermi Surfaces (AbstractFermiSurface2D)
Represent quasiparticle properties on the Fermi surface:
- **Analytic surfaces** (AbstractAnalyticSurface)
  - `Isotropic2DFermiSurface` — constant velocity in all directions
  - `EllipticFermiSurface2D` — elliptic anisotropy: vF(θ) = vF0/√(cos²θ + sin²θ/a²)
- **User-defined surfaces** (AbstractUserDefinedSurface)
  - `GeneralFermiSurface2D` — arbitrary velocity function vF(θ)

All surfaces implement:
```
surface_vF(s), surface_max_speed(s), surface_vF_angle(s, θ),
surface_density_of_states(s), surface_mass(s), surface_charge(s)
```

### Angular Discretizations (AbstractAngularDiscretization2D)
Choose how to represent angular dependence:
- `HarmonicBasis` — Fourier harmonic expansion (linear transport, efficient)
- `AngleGrid` — Direct angle sampling (nonlinear transport, exact)

### Collision Models (AbstractCollisionModel2D)
**Linear collisions** (AbstractLinearCollision) — for harmonic basis:
- `LinearBGKCollision` — BGK collision with mode-dependent scattering profile
- `LinearCollisionMatrix` — user-supplied collision matrix

**Nonlinear collisions** (AbstractNonlinearAngleCollision) — for angle grids:
- `QuadraticBGKCollision` — quadratic (parabolic) band approximation
- `ExactAngleBGKCollision` — exact angle-dependent BGK
- `TwoRateAngleBGKCollision` — simplified two-rate model
- `AngleRateBGKCollision` — mode-rate filtered angle BGK

### Mode Rate Profiles (AbstractModeRateProfile)
Define m-dependent scattering rates for harmonic modes:
- **Built-in profiles** (AbstractBuiltInProfile)
  - `TwoRateProfile` — piecewise: 0 for m<2, γ_ee for m≥2
  - `OddQuarticRateProfile` — quartic scaling for odd modes
  - `ConstantModeRateProfile` — constant for all m≥2
- **User-defined profiles** (AbstractUserProfile)
  - `CustomModeRateProfile` — arbitrary rate(m) function

### Boundary Conditions (AbstractBoundaryCondition)
Define behavior at domain boundaries:
- **Wall boundaries** (AbstractWallBC)
  - `MaxwellWallBC` — specular + diffuse reflection (p_scatter ∈ [0,1])
- **Contact boundaries** (AbstractContactBC)
  - `OhmicContactBC` — voltage-controlled carrier injection/extraction (bias parameter)
  - `CurrentContactBC` — current-controlled contact (target flux parameter)

# Usage Example

```julia
using ElectronKinetics

# Define carrier properties (e.g., graphene electrons)
surface = Isotropic2DFermiSurface(
    vF = 1.0,           # Fermi velocity
    nu = 1.0,           # density of states
    mass = 1.0,         # effective mass
    charge = -1.0       # electron charge
)

# Choose angular discretization
discretization = HarmonicBasis(:auto)  # Auto-estimate max harmonic

# Define collision model (linear BGK with e-e scattering for m≥2)
collision = LinearBGKCollision(
    0.1,                           # momentum-relaxing rate
    TwoRateProfile(0.4)            # e-e rate
)

# Assemble kinetic model
model = KineticModel2D(
    surface,
    discretization,
    IsotropicHarmonicStreaming(),
    collision
)
```

# Extensibility

The type hierarchy enables adding new physics:

**New collision model?** Just inherit from appropriate abstract type:
```julia
struct MyCollision <: AbstractLinearCollision
    gamma_mr::Float64
    profile::AbstractModeRateProfile
end
```

**New Fermi surface?** Implement the surface interface:
```julia
struct MySurface <: AbstractAnalyticSurface
    # ... fields ...
end
surface_vF(s::MySurface) = s.vF
surface_vF_angle(s::MySurface, θ) = ...
# ... etc
```

**New BC?** Add to appropriate abstract class:
```julia
mutable struct MyBC <: AbstractWallBC
    p_scatter::Float64
    tol::Float64
    cache::BCProjectorCache
end
```

# References

- Project documentation: https://fermiharmonics.jackhfarrell.com
- Julia docs: https://docs.julialang.org/
- Physics background: See project README for theory
"""
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
       AbstractAnalyticSurface,
       AbstractUserDefinedSurface,
       AbstractAngularDiscretization2D,
       AbstractStreamingOperator2D,
       AbstractModeRateProfile,
       AbstractBuiltInProfile,
       AbstractUserProfile,
       AbstractCollisionModel2D,
       AbstractLinearCollision,
       AbstractNonlinearAngleCollision,
       AbstractBoundaryCondition,
       AbstractWallBC,
       AbstractContactBC,
       AbstractTransportMode,
       LinearTransport,
       NonlinearParabolicTransport,
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
       MagneticField2D,
       LinearBGKCollision,
       LinearCollisionMatrix,
       QuadraticBGKCollision,
       ExactAngleBGKCollision,
       TwoRateAngleBGKCollision,
       AngleRateBGKCollision,
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
       CurrentContactBC,
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
       build_collision_matrix,
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
