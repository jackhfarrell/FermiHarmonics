# Quick Start Guide

Get up and running with ElectronKinetics in 5 minutes.

## Installation

```julia
using Pkg
Pkg.add("ElectronKinetics")
```

## The 30-Second Happy Path

```julia
using ElectronKinetics

# 1. Define your Fermi surface
surface = Isotropic2DFermiSurface(vF=1.0, nu=1.0, mass=1.0)

# 2. Choose angular discretization
discretization = HarmonicBasis(:auto)  # Auto-estimate modes

# 3. Define scattering
collision = LinearBGKCollision(0.1, TwoRateProfile(0.4))

# 4. Assemble model
model = KineticModel2D(
    surface,
    discretization,
    IsotropicHarmonicStreaming(),
    collision
)

println("Model ready!")
```

That's it! Your model is ready to use with a solver.

## What Each Line Does

### Define the Fermi Surface

```julia
surface = Isotropic2DFermiSurface(
    vF = 1.0,      # Fermi velocity
    nu = 1.0,      # Density of states
    mass = 1.0,    # Effective mass
    charge = -1.0  # Carrier charge (negative for electrons)
)
```

**Alternatives**:
- `EllipticFermiSurface2D(vF0=1.0, aspect=2.0, ...)` — Anisotropic
- `GeneralFermiSurface2D(θ -> 1.0, ...)` — Custom function

### Choose Angular Discretization

```julia
HarmonicBasis(:auto)        # Auto-estimate (RECOMMENDED)
HarmonicBasis(10)           # Fixed 10 modes
AngleGrid(32)               # Nonlinear? Use this instead
```

**When to use each**:
- **HarmonicBasis**: Linear transport, efficient, most cases
- **AngleGrid**: Nonlinear response, strong fields

### Define Collision (Scattering)

```julia
# Linear transport (harmonic basis)
LinearBGKCollision(
    0.1,                        # γ_mr: drag rate
    TwoRateProfile(0.4)         # γ_ee: e-e rate (m ≥ 2)
)

# Nonlinear transport (angle grid)
QuadraticBGKCollision(
    0.1,                        # γ_mr
    OddQuarticRateProfile(0.4); # γ_ee profile
    mu0 = 0.5,                  # Band bottom
    mass = 1.0                  # Band mass
)
```

### Assemble the Model

```julia
model = KineticModel2D(
    surface,                        # Your Fermi surface
    discretization,                 # Harmonic basis or angle grid
    IsotropicHarmonicStreaming(),   # Streaming (fixed by discretization)
    collision                       # Collision model
)
```

The streaming operator is determined by discretization:
- HarmonicBasis → `IsotropicHarmonicStreaming()`
- AngleGrid → `IsotropicAngleStreaming()`

## Common Tasks

### Task 1: Change Scattering Rate

```julia
# Stronger drag
collision = LinearBGKCollision(0.5, TwoRateProfile(0.4))

# Weaker drag
collision = LinearBGKCollision(0.05, TwoRateProfile(0.2))

# Different e-e rate profile
collision = LinearBGKCollision(0.1, OddQuarticRateProfile(0.4))
```

### Task 2: Use Anisotropic Surface

```julia
# Elliptic anisotropy (built-in)
surface = EllipticFermiSurface2D(
    vF0 = 1.0,      # Isotropic Fermi velocity
    aspect = 2.0    # vF_y / vF_x ratio
)

# Custom function
surface = GeneralFermiSurface2D(
    θ -> 1.0 * (1.0 + 0.3 * cos(2θ)),  # vF(θ)
    max_vF = 1.3,  # Upper bound (for CFL)
    nu = 1.0,
    mass = 1.0
)
```

### Task 3: Multiple Carriers (Multiband)

```julia
band_e = Band(:electrons, surface, gamma_mr=0.1, gamma_ee=0.05)
band_h = Band(:holes, surface, gamma_mr=0.08, gamma_ee=0.04)

model = KineticModel2D(
    surface,
    HarmonicBasis(:auto),
    IsotropicHarmonicStreaming(),
    LinearBGKCollision(0.1),
    bands = [band_e, band_h],
    gamma_drag = 0.02  # Inter-band coupling
)
```

### Task 4: Nonlinear Transport

```julia
# Must use angle grid + nonlinear collision
model = KineticModel2D(
    surface,
    AngleGrid(32),                   # Angle grid (not harmonic)
    IsotropicAngleStreaming(),       # Angle streaming
    QuadraticBGKCollision(           # Nonlinear collision
        0.05,
        OddQuarticRateProfile(0.3);
        mu0 = 0.5,
        mass = 1.0
    )
)
```

## Decision Flowchart

```
┌─────────────────────────────────────┐
│ What's your physics problem?        │
└─────────────────────────────────────┘
         │
         ├─ "Linear response"
         │  └─→ Use HarmonicBasis(:auto) + LinearBGKCollision ✓
         │
         ├─ "Nonlinear effects"
         │  └─→ Use AngleGrid(32) + QuadraticBGKCollision
         │
         ├─ "Multiple carriers"
         │  └─→ Add bands=[band1, band2] and gamma_drag=...
         │
         ├─ "Anisotropic band"
         │  └─→ Use EllipticFermiSurface2D or GeneralFermiSurface2D
         │
         └─ "Custom model?"
            └─→ See EXTENSION_GUIDE.md
```

## Next: Solving on a Mesh

Once you have your model, solve on a mesh with Trixi:

```julia
using ElectronKinetics

# Your model
model = KineticModel2D(
    Isotropic2DFermiSurface(vF=1.0),
    HarmonicBasis(:auto),
    IsotropicHarmonicStreaming(),
    LinearBGKCollision(0.1, TwoRateProfile(0.4))
)

# Define problem
problem = TrixiProblem(
    mesh_path = "your_mesh.msh",        # Path to Gmsh file
    boundary_conditions = Dict(
        :top => MaxwellWallBC(0.5),     # Specular/diffuse
        :bottom => OhmicContactBC(0.0)  # Voltage-driven contact
    )
)

# Configuration
config = SolverConfig(
    max_harmonic = 60,      # Limit on harmonics
    tspan_end = 100.0,      # End time
    cfl = 0.8,              # CFL stability parameter
    log_every = 500         # Print progress every N steps
)

# Solve!
solution = solve(problem, model, config)
```

## Getting Help

### Check the Docs

```julia
# Module overview
?ElectronKinetics

# Specific types
?LinearBGKCollision
?HarmonicBasis
?KineticModel2D

# Advanced topics
?AbstractLinearCollision
?mode_rate
```

### Examples

See `examples/` directory:
- `01_linear_harmonic_transport.jl` — Common case
- `02_nonlinear_angle_transport.jl` — Nonlinear example
- `03_multiband_transport.jl` — Multiple carriers
- `04_custom_fermi_surface.jl` — Custom surfaces

### Guides

- **PERFORMANCE.md** — Optimize your model
- **EXTENSION_GUIDE.md** — Custom components
- **examples/README.md** — Decision tree

## Troubleshooting

### Error: "HarmonicBasis requires linear collision"

**Problem**: You used `QuadraticBGKCollision` with `HarmonicBasis`

**Fix**: Use `AngleGrid` instead:
```julia
# WRONG:
model = KineticModel2D(surface, HarmonicBasis(10), ..., QuadraticBGKCollision(...))

# RIGHT:
model = KineticModel2D(surface, AngleGrid(32), ..., QuadraticBGKCollision(...))
```

### Error: "AngleGrid requires nonlinear collision"

**Problem**: You used `LinearBGKCollision` with `AngleGrid`

**Fix**: Use `QuadraticBGKCollision` or similar:
```julia
# WRONG:
model = KineticModel2D(surface, AngleGrid(32), ..., LinearBGKCollision(...))

# RIGHT:
model = KineticModel2D(surface, AngleGrid(32), ..., QuadraticBGKCollision(...))
```

### Error about mu0 or mass required

**Problem**: Nonlinear collision needs band parameters

**Fix**: Provide `mu0` and `mass` keywords:
```julia
QuadraticBGKCollision(
    0.1, OddQuarticRateProfile(0.4);
    mu0 = 0.5,       # Add these!
    mass = 1.0
)
```

### Model assembly succeeds but solve fails

**Problem**: Usually mesh-related or BC issue

**Fix**:
1. Check mesh file exists: `mesh_path = "..."`
2. Check boundary condition symbols match mesh tags
3. See solver documentation for `TrixiProblem`

## Performance Tips

1. **Start simple**: Use `HarmonicBasis(:auto)` first
2. **Small mesh**: Test with few elements before full simulation
3. **Profile carefully**: Use `SolverConfig(log_every=100)` to track progress
4. **Reduce modes if slow**: Try `HarmonicBasis(10)` if `:auto` is too many

## Minimal Working Example

```julia
using ElectronKinetics

# Simplest possible model
surface = Isotropic2DFermiSurface()
model = KineticModel2D(
    surface,
    HarmonicBasis(5),
    IsotropicHarmonicStreaming(),
    LinearBGKCollision(0.1)
)

println("Success! Model: $(typeof(model))")
```

This creates a minimal model suitable for testing and prototyping.

## What's Next?

1. ✅ You can now assemble models
2. 📖 Read examples for specific physics cases
3. ⚙️ Check PERFORMANCE.md to optimize
4. 🔧 See EXTENSION_GUIDE.md for custom components
5. 🌐 Visit docs at <https://fermiharmonics.jackhfarrell.com> for theory
