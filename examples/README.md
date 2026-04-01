# ElectronKinetics Examples

This directory contains runnable examples demonstrating common use cases of ElectronKinetics.

## Getting Started

### Running Examples

To run an example from the project root directory:

```julia
using ElectronKinetics

include("examples/01_linear_harmonic_transport.jl")
```

Or from the examples directory:

```julia
julia --project=.. 01_linear_harmonic_transport.jl
```

## Examples Overview

### 1. Linear Harmonic Basis Transport (`01_linear_harmonic_transport.jl`)

**What**: The most common use case — linear kinetic transport using efficient harmonic basis expansion.

**When to use**:
- Linear response regime (small perturbations)
- Memory-efficient representation needed
- Well-separated timescales (weak damping)
- Most routine kinetic transport problems

**Key concepts**:
- Fermi surface definition
- Harmonic basis with automatic mode estimation
- Linear BGK collision model with mode-dependent profiles
- Model assembly and component dispatch

### 2. Nonlinear Angle Grid Transport (`02_nonlinear_angle_transport.jl`)

**What**: Nonlinear kinetic transport using direct angle grid sampling.

**When to use**:
- Strong (nonlinear) driving fields
- Detailed angle-dependent physics
- Harmonic basis insufficient (saturation effects)
- Need exact angular dependence

**Key concepts**:
- Angle grid discretization
- Parabolic band collision models
- Nonlinear BGK collision variants
- Angular resolution trade-offs

### 3. Multiband Linear Transport (`03_multiband_transport.jl`)

**What**: Coupled kinetic transport for multiple carrier species.

**When to use**:
- Compensated semimetals
- Two-dimensional materials with multiple valleys
- Graphene bilayers or heterostructures
- Systems where band coupling is important

**Key concepts**:
- Band definition with individual surfaces
- Harmonic basis for multiple bands
- Inter-band drag coupling
- Multi-carrier equilibration

### 4. Custom Fermi Surface (`04_custom_fermi_surface.jl`)

**What**: Defining custom Fermi surfaces with user-supplied functions.

**When to use**:
- Non-standard band structures
- Testing theoretical predictions
- Arbitrary velocity functions
- Comparing different dispersion relations

**Key concepts**:
- Built-in analytic surfaces (elliptic)
- User-defined functions
- CFL stability bounds
- Custom model assembly

## Decision Tree: Which Example to Use?

```
Start with a kinetic transport problem
│
├─ Is it multiband (e.g., electrons + holes)?
│  └─ YES → Example 3 (multiband_transport.jl)
│
├─ Is it linear response (small driving)?
│  ├─ YES → Example 1 (linear_harmonic_transport.jl) ← MOST COMMON
│  └─ NO  → Need nonlinear?
│     └─ YES → Example 2 (nonlinear_angle_transport.jl)
│
└─ Do you need a custom Fermi surface?
   └─ YES → Example 4 (custom_fermi_surface.jl)
```

## Live Visualization

Each example accepts `--live` to enable the GLMakie live dashboard. This requires
`GLMakie` in your environment. Use `--mesh-only` (implies `--live`) to preview the
mesh and exit before the solve starts. Mesh coarseness is controlled by
`MeshBuildConfig(mesh_scale=...)` when using `.geo` inputs.

```bash
julia --project=. examples/01_linear_harmonic_transport.jl --live
```

```bash
julia --project=. examples/01_linear_harmonic_transport.jl --live --mesh-only
```

Note: the mesh preview shows the quad mesh from the `.inp`, while the live solver
dashboard visualizes a triangulated mesh-native grid.

## Next Steps

1. **Read the module docstring**: `?ElectronKinetics` for architecture overview
2. **Check type documentation**: `?LinearBGKCollision`, `?HarmonicBasis`, etc.
3. **Look at test cases**: `test/` directory for additional usage patterns
4. **Read the physics**: Project documentation at <https://fermiharmonics.jackhfarrell.com>

## Common Patterns

### Changing Resolution
```julia
# Finer harmonic basis
HarmonicBasis(15)  # Fixed 15 modes instead of :auto

# Finer angle grid
AngleGrid(64)      # 64 angles instead of 32
```

### Changing Scattering Rates
```julia
# Stronger momentum relaxation
LinearBGKCollision(0.5, TwoRateProfile(0.4))

# Weaker e-e scattering
LinearBGKCollision(0.1, TwoRateProfile(0.1))

# Custom rate profile
LinearBGKCollision(0.1, CustomModeRateProfile(m -> 0.01 * m^2))
```

### Comparing Collision Models
```julia
# Linear: efficient, linear response
collision_linear = LinearBGKCollision(0.1, TwoRateProfile(0.4))

# Nonlinear: exact angle, strong fields
collision_nonlinear = QuadraticBGKCollision(0.05, OddQuarticRateProfile(0.3); mu0=0.5, mass=1.0)

# Two-rate angle model: simplified angle dependence
collision_angle = TwoRateAngleBGKCollision(gamma_mr=0.05, gamma_ee=0.3; mu0=0.5, mass=1.0)
```

## Performance Tips

- **Harmonic basis**: Fast for weak damping (γ << 1)
- **Angle grid**: More accurate for strong fields
- **Auto-estimation**: `:auto` is usually good; adjust if needed
- **Multiband**: Scales as (N_bands × N_harmonics) in memory
- **Custom surfaces**: General functions slower than built-in types

## Troubleshooting

**Error**: "HarmonicBasis requires linear collision model"
- **Fix**: Use `LinearBGKCollision` or `LinearCollisionMatrix`, not `QuadraticBGKCollision`

**Error**: "AngleGrid requires nonlinear collision model"
- **Fix**: Use `QuadraticBGKCollision`, `ExactAngleBGKCollision`, etc., not `LinearBGKCollision`

**Error**: "collision matrix must be odd"
- **Fix**: Use `LinearCollisionMatrix(matrix)` where `size(matrix, 1)` is odd (n = 1 + 2M)

**Model slow**: Too many harmonics
- **Fix**: Try `HarmonicBasis(10)` instead of `:auto`, or use angle grid for nonlinear

## Questions or Issues?

- **Theory questions**: Check project docs at <https://fermiharmonics.jackhfarrell.com>
- **Julia questions**: <https://docs.julialang.org/>
- **Bug reports**: Create an issue with a minimal example
