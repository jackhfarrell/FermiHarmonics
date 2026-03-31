# Frequently Asked Questions

Common questions and answers about ElectronKinetics.

## Installation & Setup

### Q1: How do I install ElectronKinetics?

```julia
using Pkg
Pkg.add("ElectronKinetics")
```

Or with specific version:
```julia
Pkg.add("ElectronKinetics@0.1")
```

For development:
```julia
Pkg.develop(path="/path/to/FermiHarmonics")
```

### Q2: What are the dependencies?

Core dependencies:
- Julia 1.9+
- FFTW.jl
- LinearAlgebra (standard library)
- SparseArrays (standard library)
- StaticArrays.jl

For solving on meshes:
- Trixi.jl (loaded via extension)
- SciMLBase.jl
- OrdinaryDiffEq.jl
- HDF5.jl

## Choosing Components

### Q3: Should I use HarmonicBasis or AngleGrid?

**Use HarmonicBasis if**:
- Linear response (small perturbations)
- Weak scattering (γ < 1)
- Need fast computation
- Memory is limited
- **Most cases** ← default choice

**Use AngleGrid if**:
- Nonlinear response (strong fields)
- Saturation/heating effects
- Need exact angle dependence
- Have memory/time budget

See MODEL_COMPARISON.md for detailed comparison.

### Q4: What collision model should I use?

**For linear transport (HarmonicBasis)**:
- `LinearBGKCollision` ← **Most common**
- `LinearCollisionMatrix` — only if you have pre-computed matrix

**For nonlinear transport (AngleGrid)**:
- `QuadraticBGKCollision` ← **Default nonlinear**
- `ExactAngleBGKCollision` — exact angle dependence
- `TwoRateAngleBGKCollision` — simplified model
- `AngleRateBGKCollision` — mode-rate filtered

See MODEL_COMPARISON.md for when to use each.

### Q5: When do I need to specify mu0 and mass?

**You need mu0 and mass when**:
- Using `QuadraticBGKCollision` or other nonlinear collision
- Modeling parabolic band approximation
- Working with angle grid discretization

**You don't need them for**:
- `LinearBGKCollision` (they're optional)
- Linear transport generally

Example:
```julia
# Linear: mu0 and mass optional
LinearBGKCollision(0.1, TwoRateProfile(0.4))

# Nonlinear: mu0 and mass REQUIRED
QuadraticBGKCollision(
    0.1, OddQuarticRateProfile(0.4);
    mu0=0.5,  # Required!
    mass=1.0  # Required!
)
```

### Q6: What's the difference between TwoRateProfile and OddQuarticRateProfile?

| Profile | Rate Function | When to Use |
|---------|---------------|------------|
| `TwoRateProfile(γ_ee)` | γ_ee for m≥2, 0 for m<2 | **Most common**, default choice |
| `OddQuarticRateProfile(γ_ee)` | γ_ee·m² for odd m | Nonequilibrium, heating |
| `ConstantModeRateProfile(γ)` | γ for all m≥2 | Uniform damping |
| `CustomModeRateProfile(f)` | Custom function f(m) | Arbitrary physics |

In most cases, start with `TwoRateProfile`.

### Q7: How do I model electron-hole plasma (multiband)?

Define bands with their properties, then pass to model:

```julia
electron = Band(
    :e,
    Isotropic2DFermiSurface(vF=1.0, nu=0.5, mass=0.5, charge=-1.0),
    gamma_mr=0.1,
    gamma_ee=0.05
)

hole = Band(
    :h,
    Isotropic2DFermiSurface(vF=0.8, nu=0.5, mass=0.6, charge=1.0),
    gamma_mr=0.08,
    gamma_ee=0.04
)

model = KineticModel2D(
    electron.surface,
    HarmonicBasis(:auto),
    IsotropicHarmonicStreaming(),
    LinearBGKCollision(0.1),
    bands=[electron, hole],
    gamma_drag=0.02  # Inter-band coupling
)
```

Band names must be unique.

## Model Assembly & Validation

### Q8: Why does my model fail with "HarmonicBasis requires linear collision"?

You mixed incompatible components:

```julia
# WRONG: QuadraticBGKCollision with HarmonicBasis
model = KineticModel2D(
    surface,
    HarmonicBasis(10),  # Linear-only!
    IsotropicHarmonicStreaming(),
    QuadraticBGKCollision(...)  # Nonlinear collision
)
```

**Fix**: Use AngleGrid instead:
```julia
# RIGHT: Nonlinear collision needs AngleGrid
model = KineticModel2D(
    surface,
    AngleGrid(32),  # Angle grid for nonlinear
    IsotropicAngleStreaming(),
    QuadraticBGKCollision(...)
)
```

### Q9: Why does my model fail with "AngleGrid requires nonlinear collision"?

You mixed incompatible components:

```julia
# WRONG: LinearBGKCollision with AngleGrid
model = KineticModel2D(
    surface,
    AngleGrid(32),  # Nonlinear-only!
    IsotropicAngleStreaming(),
    LinearBGKCollision(...)  # Linear collision
)
```

**Fix**: Use HarmonicBasis instead:
```julia
# RIGHT: Linear collision needs HarmonicBasis
model = KineticModel2D(
    surface,
    HarmonicBasis(:auto),  # Harmonic basis for linear
    IsotropicHarmonicStreaming(),
    LinearBGKCollision(...)
)
```

The rule:
- HarmonicBasis ↔ LinearBGKCollision / LinearCollisionMatrix
- AngleGrid ↔ QuadraticBGKCollision / ExactAngleBGKCollision / etc.

### Q10: What does max_vF mean in GeneralFermiSurface2D?

`max_vF` is the tight **upper bound** on vF(θ):

```julia
surface = GeneralFermiSurface2D(
    θ -> 1.0 * (1.0 + 0.3 * cos(2θ)),  # vF(θ)
    max_vF = 1.3,  # Must satisfy: max_vF ≥ vF(θ) for all θ
    nu = 1.0
)
```

**Why it matters**: Used for CFL stability. Too loose → unstable. Too tight → wasteful.

**How to choose**:
```julia
# Compute actual max
vF_func(θ) = 1.0 * (1.0 + 0.3 * cos(2θ))
θ_test = range(0, 2π, 1000)
v_test = vF_func.(θ_test)
max_vF = maximum(v_test) * 1.01  # Add 1% buffer

surface = GeneralFermiSurface2D(vF_func, max_vF=max_vF, ...)
```

## Physical Parameters

### Q11: What units should I use?

ElectronKinetics is **dimensionless**. Choose any consistent system:

**Option 1: Energy units (typical)**
```julia
# γ in units of energy/ℏ
# vF in units of ℏk_F/m
# mu0, mass in same energy units
# Time: ℏ/energy_unit
LinearBGKCollision(0.1, TwoRateProfile(0.4))  # γ = 0.1 meV (if energy units are meV)
```

**Option 2: Length-time units**
```julia
# Velocities in units of distance/time
# Rates in units of 1/time
# Time: in chosen time units
```

**Just be consistent!**

### Q12: How do I choose scattering rates γ_mr and γ_ee?

From transport measurements or theory:

**From experiment**:
```julia
# Measure momentum relaxation time τ_mr
# γ_mr = 1/τ_mr = 1/τ_momentum_relaxation
LinearBGKCollision(gamma_mr=1/tau_mr, ...)
```

**From theory**:
```julia
# Coulomb scattering: γ ∝ T²
# Impurity scattering: γ = const
# Phonon scattering: γ ∝ T³
# Your model here...
```

**Typical ranges**:
- Clean: γ_mr ~ 0.01-0.1
- Intermediate: γ_mr ~ 0.1-1.0
- Dirty: γ_mr > 1.0

For e-e: `gamma_ee` typically 0.1-1.0 times γ_mr.

### Q13: What does electrostatic_coupling (chi) do?

Coupling parameter for electron-electron interaction:

```julia
QuadraticBGKCollision(
    0.1, OddQuarticRateProfile(0.4);
    mu0=0.5,
    mass=1.0,
    electrostatic_coupling=0.1  # χ parameter
)
```

**Physical meaning**:
- χ = 0: No e-e coupling (default)
- χ > 0: Screening effects
- χ ~ 1: Strong screening
- χ >> 1: Very strong coupling

Affects nonlinear transport and heating. Usually set to 0 unless studying e-e effects.

## Running Simulations

### Q14: Why does my simulation crash or give wrong results?

Common issues and fixes:

**1. Model-mesh mismatch**
```julia
# Problem: Boundary condition symbols don't match mesh
problem = TrixiProblem(
    mesh_path="mesh.msh",
    boundary_conditions=Dict(
        :top => MaxwellWallBC(0.5)      # Must match mesh tag names!
    )
)

# Fix: Check mesh tags in Gmsh
# gmsh meshfile.geo -info
```

**2. Timestep too large**
```julia
# CFL condition: dt < CFL * dx / v_max
config = SolverConfig(
    cfl=0.8,  # Default; reduce to 0.4 if unstable
    tspan_end=100.0
)
```

**3. Too few harmonics**
```julia
# Check: Are modes damped at cutoff?
discretization = HarmonicBasis(20)  # Increase if needed

# Or auto-estimate
discretization = HarmonicBasis(:auto)
```

**4. Boundary condition incompatible**
```julia
# Check dimension compatibility
# Wall BC needs position on boundary
# Contact BC needs voltage/current specification
```

### Q15: How do I know if my simulation converged?

Check progress output:

```julia
config = SolverConfig(
    log_every=500,      # Print every 500 steps
    residual_tol=1e-5   # Convergence tolerance
)

# Watch the `residual_progress_fraction` value
# Should decrease monotonically
# Should reach ~1.0 at end
```

**In solve output**, look for:
- Residual decreasing
- No NaNs or Infs
- Run time reasonable
- Final state makes physical sense

## Troubleshooting & Debugging

### Q16: How do I profile my code for bottlenecks?

```julia
using BenchmarkTools, ProfileView

# Benchmark key functions
@benchmark streaming_matrices(10, $surface)
@benchmark mode_rate($profile, 5)
@benchmark collision_gamma_mr($collision)

# Profile a solve
@profview solve(problem, model, config)
```

See PERFORMANCE.md for detailed profiling guide.

### Q17: How do I get better error messages?

ElectronKinetics uses structured `PhysicsError` type:

```julia
try
    bad_collision = LinearBGKCollision(-0.1)  # Negative rate!
catch err
    if err isa PhysicsError
        println("Problem: $(err.problem)")
        println("Details: $(err.details)")
        println("Suggestion: $(err.suggestion)")
    end
end
```

All validation uses this structured format to explain:
1. What went wrong
2. What was provided
3. How to fix it
4. Where to find more info

### Q18: How do I test custom components?

```julia
using InteractiveUtils

# Type stability check (no allocation warnings)
@code_warntype surface_vF(my_surface)

# Test mode rate function
profile = MyModeRateProfile(0.4, 0.1)
@test mode_rate(profile, 1) == 0.0
@test mode_rate(profile, 2) > 0.0

# Test model assembly
model = KineticModel2D(my_surface, HarmonicBasis(10), ..., collision)
@test typeof(model) <: KineticModel2D
```

See EXTENSION_GUIDE.md for more testing patterns.

## Advanced Questions

### Q19: Can I use weak dependencies to extend?

Yes! ElectronKinetics uses extension mechanism:

```julia
# In ext/MyExtension.jl
module MyExtensionModule
using ElectronKinetics
using MyDependency

# Add your custom implementations here
struct MyCollision <: AbstractLinearCollision
    # ...
end
end
```

Then declare in Project.toml:
```toml
[weakdeps]
MyDependency = "uuid-here"

[extensions]
MyExtensionModule = "MyDependency"
```

See EXTENSION_GUIDE.md for full examples.

### Q20: How do I contribute improvements?

1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Submit a pull request with clear description
5. Ensure all tests pass

See CLAUDE.md for development guidelines.

## Performance Questions

### Q21: My simulation is too slow. How do I optimize?

1. **Reduce harmonics**:
   ```julia
   HarmonicBasis(10)  # Instead of :auto
   ```

2. **Use faster surface**:
   ```julia
   Isotropic2DFermiSurface  # Fastest
   EllipticFermiSurface2D   # Medium
   GeneralFermiSurface2D    # Slowest
   ```

3. **Coarser mesh** (for testing):
   ```julia
   # First test on small mesh
   # Then scale up
   ```

4. **Larger CFL**:
   ```julia
   config = SolverConfig(cfl=0.9)  # More aggressive
   ```

See PERFORMANCE.md for detailed optimization guide.

### Q22: How much memory does my model use?

Rough estimates:

```
HarmonicBasis(M):        5 + 10*M MB
AngleGrid(nθ):           10 + 10*nθ MB
Multiband(N):            × N MB
Custom surface:          no extra memory
```

For N=4 bands, M=20 harmonics: ~200 MB total.

For detailed breakdown, see PERFORMANCE.md.

## Conceptual Questions

### Q23: What's the difference between elastic and inelastic scattering?

In ElectronKinetics:
- **Elastic**: Momentum-relaxing (γ_mr)
- **Inelastic**: Energy-relaxing (included via nonlinear collision)

The library focuses on momentum-relaxing processes. Energy-relaxing processes affect temperature evolution in some models.

### Q24: Why is there a difference between γ_mr and γ_ee?

**γ_mr** (momentum-relaxing):
- Drops carrier out of specific direction
- Decays first harmonic fast
- Dominates transport

**γ_ee** (electron-electron):
- Damps higher harmonics (m ≥ 2)
- Less drastic than γ_mr
- Couple through mode-dependent profile

See physics literature for detailed explanation.

### Q25: What does "harmonic basis" really mean?

Fourier expansion of angular dependence:

```
f(θ) ≈ a₀ + Σ[a_m cos(mθ) + b_m sin(mθ)]
```

- a₀: Isotropic part
- a₁, b₁: Dipole (linear response)
- a₂, b₂: Quadrupole (e-e effects)
- Higher m: Progressively damped

HarmonicBasis(M) truncates at m = M.

## Getting More Help

- **Examples**: See `examples/` directory
- **Guides**: QUICKSTART.md, MODEL_COMPARISON.md, EXTENSION_GUIDE.md
- **Theory**: <https://fermiharmonics.jackhfarrell.com>
- **Julia docs**: <https://docs.julialang.org/>
- **Type docs**: `?YourType` in Julia REPL

## Still Have Questions?

Check:
1. **Relevant guide**: Look in docs/ or top-level md files
2. **Examples**: Run relevant example and modify
3. **API docs**: `?TypeName` in Julia
4. **Source code**: Well-commented, readable
5. **GitHub issues**: Check existing issues or open new one

Remember: All questions are valid! Feel free to ask for clarification.
