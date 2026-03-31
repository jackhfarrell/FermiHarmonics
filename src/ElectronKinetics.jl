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

# Extensibility Guide

The type hierarchy enables adding custom physics without modifying core code.

## Adding a Custom Collision Model

For linear harmonic transport:
```julia
struct MyLinearCollision <: AbstractLinearCollision
    gamma_mr::Float64
    # Add your custom parameters here
end

# Constructor with validation
function MyLinearCollision(gamma_mr::Real; my_param::Real)
    return MyLinearCollision(
        Float64(gamma_mr) >= 0 ? Float64(gamma_mr) : error("gamma_mr must be ≥ 0")
    )
end
```

For nonlinear angle transport:
```julia
struct MyNonlinearCollision <: AbstractNonlinearAngleCollision
    gamma_mr::Float64
    mu0::Float64
    mass::Float64
    electrostatic_coupling::Float64
end
```

## Adding a Custom Fermi Surface

For analytic surfaces with known formulas:
```julia
struct MyAnalyticSurface <: AbstractAnalyticSurface
    vF0::Float64
    anisotropy_param::Float64
    nu::Float64
    mass::Float64
    charge::Float64
end

# Implement the required interface (6 methods):
surface_vF(s::MyAnalyticSurface) = s.vF0
surface_max_speed(s::MyAnalyticSurface) = compute_max_speed(s)
surface_vF_angle(s::MyAnalyticSurface, θ) = s.vF0 * my_formula(θ, s.anisotropy_param)
surface_density_of_states(s::MyAnalyticSurface) = s.nu
surface_mass(s::MyAnalyticSurface) = s.mass
surface_charge(s::MyAnalyticSurface) = s.charge
```

For user-defined functions:
```julia
struct MyCustomSurface <: AbstractUserDefinedSurface
    vF_func::Function  # θ -> vF
    max_vF::Float64    # CFL bound
    nu::Float64
    mass::Float64
    charge::Float64
end
```

## Adding a Custom Boundary Condition

For wall-type boundaries:
```julia
mutable struct MyWallBC <: AbstractWallBC
    p_scatter::Float64
    my_parameter::Float64
    tol::Float64
    cache::BCProjectorCache
end
```

For contact-type boundaries:
```julia
mutable struct MyContactBC <: AbstractContactBC
    p_ohmic_absorb::Float64
    my_control_param::Float64
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

# ============================================================================
# Core Physics API — types and solvers users directly instantiate
# ============================================================================
export
    # Fermi Surfaces
    Isotropic2DFermiSurface,
    EllipticFermiSurface2D,
    GeneralFermiSurface2D,

    # Angular Discretizations
    HarmonicBasis,
    AngleGrid,

    # Streaming Operators
    IsotropicHarmonicStreaming,
    IsotropicAngleStreaming,

    # Collision Models (Linear)
    LinearBGKCollision,
    LinearCollisionMatrix,

    # Collision Models (Nonlinear)
    QuadraticBGKCollision,
    ExactAngleBGKCollision,
    TwoRateAngleBGKCollision,
    AngleRateBGKCollision,

    # Mode Rate Profiles
    TwoRateProfile,
    OddQuarticRateProfile,
    ConstantModeRateProfile,
    CustomModeRateProfile,

    # Boundary Conditions
    MaxwellWallBC,
    OhmicContactBC,
    CurrentContactBC,

    # Core Models
    Band,
    KineticModel2D,
    MagneticField2D,

    # Configuration
    SolverConfig,
    SolveParams,
    LiveVisualizationConfig,
    MeshBuildConfig,
    TrixiProblem,

    # Main Solver Interface
    solve,
    solve_status

# ============================================================================
# Reference Utilities
# ============================================================================
export
    blg_reference_setup

# ============================================================================
# Analysis and Output
# ============================================================================
export
    save_solution_custom,
    save_for_analysis,
    save_mesh_native_analysis,
    evaluate_solution,
    evaluate_observables,
    LiveProgressSnapshot,
    LiveFieldSnapshot,
    LiveVisualizationSnapshot

# ============================================================================
# Mesh Generation
# ============================================================================
export
    generate_mesh_from_geo,
    mesh_provenance_attributes,
    resolve_mesh_path

# ============================================================================
# HPC/Sweep Utilities
# ============================================================================
export
    submit_sweep!,
    write_sweep_metadata!,
    archive_mesh!,
    copy_mesh_to_scratch,
    select_cases,
    grid_lookup,
    ordered_case_indices

# ============================================================================
# Advanced API — For custom implementations and extensions
# ============================================================================
export
    # Abstract types for extensibility
    AbstractFermiSurface2D,
    AbstractAnalyticSurface,
    AbstractUserDefinedSurface,
    AbstractAngularDiscretization2D,
    AbstractHarmonicDiscretization,
    AbstractGridDiscretization,
    AbstractStreamingOperator2D,
    AbstractIsotropicStreaming,
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

    # Advanced helper functions (for custom implementations)
    mode_rate,
    build_collision_matrix,
    estimate_max_harmonic,
    streaming_matrices,
    harmonic_state_nvars,
    band_momentum_weight,

    # Surface interface (implement these for custom surfaces)
    surface_vF,
    surface_max_speed,
    surface_vF_angle,
    surface_density_of_states,
    surface_mass,
    surface_charge,

    # Low-level utilities
    cosine_index,
    sine_index,
    residual_progress_fraction,
    collision_sources!

end
