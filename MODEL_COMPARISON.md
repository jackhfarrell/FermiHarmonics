# Model Comparison: Choosing the Right Components

This guide helps you choose between different physics components based on your problem.

## Angular Discretization: HarmonicBasis vs AngleGrid

| Feature | HarmonicBasis | AngleGrid |
|---------|---------------|-----------|
| **Use Case** | Linear response, weak damping | Nonlinear transport, strong fields |
| **Speed** | ⚡ Fast | 🐢 Slower |
| **Memory** | 💾 Low | 💾 Medium-High |
| **Accuracy** | High for linear | Exact for arbitrary angles |
| **Requires** | `LinearBGKCollision` | `QuadraticBGKCollision` variants |
| **Modes/Angles** | 5-60 modes (small M) | 32-128 angles (larger) |
| **Time stepping** | Standard RK methods | May need smaller dt |
| **Example** | `HarmonicBasis(:auto)` | `AngleGrid(32)` |

### When to Choose HarmonicBasis

✅ **Use HarmonicBasis if**:
- Small perturbations (linear regime)
- Weak scattering (γ < 1)
- You need speed/efficiency
- Memory is limited
- Most routine simulations

❌ **Don't use if**:
- Strong driving fields (nonlinear)
- You need exact angle dependence
- Saturation effects important

### When to Choose AngleGrid

✅ **Use AngleGrid if**:
- Large driving fields
- Need nonlinear response
- Detailed angle dependence
- Saturation/heating effects
- Have memory/time budget

❌ **Don't use if**:
- Only linear response needed
- Memory is very limited
- Simulation must be fast

## Collision Models: Which One?

### Linear Collision Models (with HarmonicBasis)

| Model | Parameters | When to Use | Notes |
|-------|-----------|-------------|-------|
| `LinearBGKCollision` | `gamma_mr`, `profile` | **Most common** | BGK with mode-dependent rates |
| `LinearCollisionMatrix` | Explicit matrix | Custom rates | Pre-computed collision matrix |

**LinearBGKCollision Details**:
```julia
LinearBGKCollision(
    gamma_mr::Float64,              # Momentum-relaxing rate
    profile::AbstractModeRateProfile # Mode-dependent e-e rate
)
```

**Choose Profile**:

| Profile | Formula | When to Use |
|---------|---------|-------------|
| `TwoRateProfile(γ_ee)` | γ_ee(m≥2), 0 (m<2) | **Most common**, weak e-e |
| `OddQuarticRateProfile(γ_ee)` | γ_ee·m² (odd m) | Quartic band, nonequilibrium |
| `ConstantModeRateProfile(γ)` | γ (all m≥2) | Uniform scattering |
| `CustomModeRateProfile(f)` | Custom f(m) | Arbitrary rate function |

### Nonlinear Collision Models (with AngleGrid)

| Model | Parameters | When to Use |
|-------|-----------|-------------|
| `QuadraticBGKCollision` | `gamma_mr`, `profile`, `mu0`, `mass` | **Standard nonlinear** |
| `ExactAngleBGKCollision` | `gamma_mr`, `gamma_ee`, `mu0`, `mass` | Exact angle-dependent rates |
| `TwoRateAngleBGKCollision` | `gamma_mr`, `gamma_ee`, `mu0`, `mass` | Simplified two-rate model |
| `AngleRateBGKCollision` | `gamma_mr`, `profile`, `mu0`, `mass` | Mode-rate filtered |

**QuadraticBGKCollision Details** (most common):
```julia
QuadraticBGKCollision(
    gamma_mr::Float64,                      # Momentum relaxation
    profile::AbstractModeRateProfile;       # Mode-dependent rates
    mu0::Float64,                           # Band bottom energy
    mass::Float64,                          # Band mass
    electrostatic_coupling::Float64=0.0,    # e-e interaction
    theta_oversample::Int=1                 # FFT oversampling
)
```

## Fermi Surface: Shape and Isotropy

| Surface | Isotropy | Formula | Memory | Speed | When to Use |
|---------|----------|---------|--------|-------|------------|
| `Isotropic2DFermiSurface` | ✅ Isotropic | vF = const | Minimal | ⚡ Fast | **Default choice** |
| `EllipticFermiSurface2D` | ⚠️ Elliptic | vF(θ) = vF₀/√(...) | Minimal | 🔥 Fast | 2D anisotropy |
| `GeneralFermiSurface2D` | ❌ Arbitrary | vF(θ) = user fn | Minimal | 🐢 Slower | Custom bands |

### When to Choose Each

**Isotropic2DFermiSurface**:
- Graphene, Dirac materials
- Circular Fermi surface
- No preferred direction
- **Use by default**

**EllipticFermiSurface2D**:
- Elliptic Fermi pockets
- Different vF in x vs y
- Uniaxial anisotropy
- Transition metal dichalcogenides

**GeneralFermiSurface2D**:
- Hexagonal Fermi surfaces
- Triangular/square lattice
- Custom dispersion relations
- Comparing theory vs experiment

## Complete Decision Tables

### Linear Transport (HarmonicBasis)

Choose this combination for most problems:

```julia
# CASE 1: Simple linear transport (MOST COMMON)
model = KineticModel2D(
    Isotropic2DFermiSurface(vF=1.0, nu=1.0),
    HarmonicBasis(:auto),
    IsotropicHarmonicStreaming(),
    LinearBGKCollision(0.1, TwoRateProfile(0.4))
)

# CASE 2: Weak e-e scattering (clean limit)
model = KineticModel2D(
    Isotropic2DFermiSurface(vF=1.0),
    HarmonicBasis(:auto),
    IsotropicHarmonicStreaming(),
    LinearBGKCollision(0.1, TwoRateProfile(0.05))  # Weak γ_ee
)

# CASE 3: Anisotropic band
model = KineticModel2D(
    EllipticFermiSurface2D(vF0=1.0, aspect=2.0),
    HarmonicBasis(:auto),
    IsotropicHarmonicStreaming(),
    LinearBGKCollision(0.1, TwoRateProfile(0.4))
)

# CASE 4: Custom Fermi surface
model = KineticModel2D(
    GeneralFermiSurface2D(θ -> 1.0 + 0.1*cos(4θ), max_vF=1.1, nu=1.0),
    HarmonicBasis(:auto),
    IsotropicHarmonicStreaming(),
    LinearBGKCollision(0.1, TwoRateProfile(0.4))
)

# CASE 5: Multi-carrier system
model = KineticModel2D(
    surface,
    HarmonicBasis(:auto),
    IsotropicHarmonicStreaming(),
    LinearBGKCollision(0.1),
    bands=[electron_band, hole_band],
    gamma_drag=0.02
)
```

### Nonlinear Transport (AngleGrid)

Use for strong fields and nonlinear effects:

```julia
# CASE 1: Standard nonlinear transport
model = KineticModel2D(
    Isotropic2DFermiSurface(vF=1.0),
    AngleGrid(32),
    IsotropicAngleStreaming(),
    QuadraticBGKCollision(
        0.05, OddQuarticRateProfile(0.3);
        mu0=0.5, mass=1.0
    )
)

# CASE 2: High-resolution (more angles)
model = KineticModel2D(
    surface,
    AngleGrid(64),  # Higher resolution
    IsotropicAngleStreaming(),
    QuadraticBGKCollision(0.05, OddQuarticRateProfile(0.3); mu0=0.5, mass=1.0)
)

# CASE 3: Exact angle-dependent rates
model = KineticModel2D(
    surface,
    AngleGrid(32),
    IsotropicAngleStreaming(),
    ExactAngleBGKCollision(gamma_mr=0.05, gamma_ee=0.3; mu0=0.5, mass=1.0)
)

# CASE 4: With e-e coupling
model = KineticModel2D(
    surface,
    AngleGrid(32),
    IsotropicAngleStreaming(),
    QuadraticBGKCollision(
        0.05, OddQuarticRateProfile(0.3);
        mu0=0.5, mass=1.0,
        electrostatic_coupling=0.1  # Add e-e interaction
    )
)
```

## Comparison: Scattering Rates

### Effect of γ_mr (momentum relaxation)

```julia
# Weak drag (clean limit)
LinearBGKCollision(0.01, TwoRateProfile(0.05))

# Moderate drag
LinearBGKCollision(0.1, TwoRateProfile(0.4))

# Strong drag (dirty limit)
LinearBGKCollision(1.0, TwoRateProfile(2.0))
```

**Impact**:
- Small γ_mr: Long scattering time, sharp features
- Large γ_mr: Fast relaxation, broadened features

### Effect of γ_ee (electron-electron)

```julia
# Weak e-e (low density, high T)
TwoRateProfile(0.1)

# Moderate e-e
TwoRateProfile(0.4)

# Strong e-e (high density, low T)
TwoRateProfile(1.0)
```

**Impact**:
- Damps higher harmonics (m ≥ 2)
- Weak e-e: Many harmonics needed
- Strong e-e: Fewer harmonics sufficient

## Quick Reference: Which Model for...?

| Problem | Surface | Discretization | Collision |
|---------|---------|-----------------|-----------|
| Graphene transport | Isotropic | HarmonicBasis(:auto) | LinearBGKCollision |
| Bilayer graphene | Isotropic | HarmonicBasis(:auto) | LinearBGKCollision + bands |
| TMD valley | Isotropic | HarmonicBasis(:auto) | LinearBGKCollision |
| Strong field | Isotropic | AngleGrid(32) | QuadraticBGKCollision |
| Anisotropic material | EllipticFermi | HarmonicBasis(:auto) | LinearBGKCollision |
| Arbitrary band | GeneralFermi | HarmonicBasis(:auto) | LinearBGKCollision |
| Nonlinear Hall | Isotropic | AngleGrid(32) | QuadraticBGKCollision |
| Clean limit | Isotropic | HarmonicBasis(20) | LinearBGKCollision(0.01) |
| Dirty limit | Isotropic | HarmonicBasis(5) | LinearBGKCollision(1.0) |

## Memory/Speed Trade-offs

### To make model FASTER:

1. Reduce modes: `HarmonicBasis(10)` instead of `:auto`
2. Use isotropic: `Isotropic2DFermiSurface` not `GeneralFermiSurface2D`
3. Reduce angles: `AngleGrid(32)` instead of `64`
4. Higher scattering: More modes damped naturally

### To make model MORE ACCURATE:

1. Increase modes: `HarmonicBasis(30)` for weak damping
2. Use anisotropic surface if needed
3. More angles: `AngleGrid(64)` for fine angle structure
4. Finer mesh: More spatial resolution

### To reduce MEMORY:

1. Lower modes: `HarmonicBasis(5)` minimum
2. Fewer bands: Combine if possible
3. Angle grid: Fixed θ resolution
4. Single-carrier: Avoid multiband

## Validation Checklist

Before running a solve:

- [ ] Surface properties reasonable? (vF > 0, nu > 0, mass > 0)
- [ ] Harmonics or angles chosen? (not both)
- [ ] HarmonicBasis ↔ LinearBGKCollision? (or AngleGrid ↔ nonlinear)
- [ ] Nonlinear collision has mu0, mass?
- [ ] Band names unique (multiband)?
- [ ] max_vF tight upper bound (custom surface)?
- [ ] Boundary conditions match mesh tags?

## Performance Expectations

Typical timings (desktop CPU, 2D mesh ~10k elements):

| Configuration | Time per step | Memory |
|---|---|---|
| Harmonic(10) + Linear | 1-5 ms | 50 MB |
| Harmonic(20) + Linear | 5-20 ms | 150 MB |
| Harmonic(30) + Linear | 20-50 ms | 300 MB |
| AngleGrid(32) + Quad | 10-30 ms | 200 MB |
| Multiband(2) + Harmonic(10) | 2-10 ms | 100 MB |

Scale as: Time ~ (M² or N_angles)² for spatial discretization.

## More Information

- **Quick-start**: QUICKSTART.md
- **Examples**: examples/ directory
- **Performance**: PERFORMANCE.md
- **Extensions**: EXTENSION_GUIDE.md
- **Theory**: <https://fermiharmonics.jackhfarrell.com>
