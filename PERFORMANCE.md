# Performance and Type Stability Guide

This document covers performance characteristics and optimization strategies for ElectronKinetics.

## Type Stability Analysis

The library is designed for **type stability** — the Julia compiler can infer concrete types for all variables at compile time, enabling specialization and elimination of runtime branches.

### Core Design Principles

1. **Abstract Dispatch**: Type decisions happen at compile time (model construction), not runtime
2. **Concrete Structures**: All physics types are concrete (not unions) when fully specified
3. **Parametric Types**: Type parameters carry physics information (e.g., `FermiHarmonics2D{NVARS, LinearTransport, ...}`)

### Type-Stable Patterns

#### ✅ Good: Dispatch at construction time

```julia
# Type fully determined at construction
model = KineticModel2D(surface, HarmonicBasis(10), IsotropicHarmonicStreaming(), collision)
# typeof(model) is concrete, all derived types are concrete
```

#### ✅ Good: Generic functions on abstract types

```julia
function my_solver(model::KineticModel2D)
    # Dispatch resolves once at call time
    _solve_dispatch(model.discretization, model)
end

function _solve_dispatch(::HarmonicBasis, model::KineticModel2D)
    # Specialized code for harmonic basis
end
```

#### ❌ Avoid: Type unions in hot paths

```julia
# BAD: Union type forces runtime checking
function bad_solver(collision::Union{LinearBGKCollision, QuadraticBGKCollision})
    # Compiler cannot specialize: runtime type check required
end

# GOOD: Use dispatch instead
function good_solver(collision::LinearBGKCollision)
    # ...
end

function good_solver(collision::QuadraticBGKCollision)
    # ...
end
```

## Performance Characteristics

### Memory Usage

| Component | Memory/Item | Notes |
|-----------|-----------|-------|
| HarmonicBasis(M) | `(1 + 2M) × F64` per element | Linear in max harmonic |
| AngleGrid(nθ) | `nθ × F64` per element | Fixed angle resolution |
| Multiband N | `N × harmonic_vars` | Scales linearly with bands |
| Streaming matrix | `(1+2M)² × F64` | ~10-20 KB for M=20 |

### Speed Ranking

1. **Isotropic2DFermiSurface** (fastest)
   - Direct formulas, no branches
   - Streaming matrices via closed-form
   - Cost: O(M²) for M harmonics

2. **EllipticFermiSurface2D** (fast)
   - Parametric formula
   - Streaming matrices via 1D integration
   - Cost: O(M² × quadrature points)

3. **GeneralFermiSurface2D** (slower)
   - Arbitrary function evaluation
   - Requires 2D integration for streaming
   - Cost: O(M² × nquad²)

4. **AngleGrid** (memory-intensive)
   - Direct sampling: no integration
   - FFT-based harmonic transforms
   - Cost: O(nθ log nθ) per transform

### Solver Scaling

For typical problems:

- **Linear transport, M=10**: ~1 ms per timestep on modern CPU
- **Linear transport, M=30**: ~10 ms per timestep
- **Angle grid, nθ=32**: ~5-10 ms per timestep (nonlinear)
- **Multiband N=2, M=10**: ~2 ms per timestep (2 bands)
- **Multiband N=3, M=10**: ~3 ms per timestep (3 bands)

## Optimization Strategies

### 1. Choose Right Discretization

**If you have**:
- [x] Linear regime (small perturbations)
- [x] Time/memory constraints
- [x] Well-separated scales

→ **Use HarmonicBasis(:auto)**

**If you need**:
- [x] Nonlinear response (strong fields)
- [x] Exact angle dependence
- [x] Can afford more memory

→ **Use AngleGrid(64)**

### 2. Harmonic Basis Optimization

```julia
# Auto-estimation (recommended)
HarmonicBasis(:auto)  # Uses collision rates to estimate

# Manual control when auto is wrong
HarmonicBasis(10)     # Fixed 10 modes
HarmonicBasis(20)     # Finer resolution (slower)

# Rough rule of thumb:
# M_max ≈ 10 - 20 for weak scattering (γ < 0.1)
# M_max ≈ 5 - 10 for moderate scattering (0.1 < γ < 1)
# M_max ≈ 2 - 5 for strong scattering (γ > 1)
```

### 3. Surface Choice

```julia
# Fastest (preferred if applicable)
surface = Isotropic2DFermiSurface(vF=1.0, ...)

# Nearly as fast (anisotropic systems)
surface = EllipticFermiSurface2D(vF0=1.0, aspect=2.0, ...)

# Slower (arbitrary band structures)
surface = GeneralFermiSurface2D(θ -> 1.0 + 0.1*cos(2θ), ...)
```

### 4. Memory Efficiency

```julia
# Reduce harmonics if possible
HarmonicBasis(10)  # Less memory than 20

# For multiband, use minimal number of bands
KineticModel2D(..., bands=[band1, band2])  # Not 10 bands if 2 suffice

# Angle grid: use minimum angles needed for accuracy
AngleGrid(32)   # Usually sufficient
AngleGrid(64)   # Only if high accuracy needed
AngleGrid(128)  # Rarely necessary
```

## Profiling and Debugging

### Check Type Stability

```julia
using InteractiveUtils

# Create model and collision
model = KineticModel2D(surface, HarmonicBasis(10), IsotropicHarmonicStreaming(), collision)

# Type-check key functions (look for Type Instability warnings)
@code_warntype streaming_matrices(10, surface)
@code_warntype mode_rate(profile, 5)
@code_warntype collision_gamma_mr(collision)
```

### Profile Allocations

```julia
using ProfileView

# Monitor allocations during solve
model = KineticModel2D(...)
problem = TrixiProblem(mesh_path="mesh.msh", boundary_conditions=...)
config = SolverConfig(tspan_end=10.0)

# Run under allocation profiler
@profview sol = solve(problem, model, config)
# Identify which functions allocate excessively
```

### Benchmark Components

```julia
using BenchmarkTools

# Benchmark streaming matrix computation
surface = Isotropic2DFermiSurface(vF=1.0)
profile = TwoRateProfile(0.4)

@benchmark streaming_matrices(10, $surface)
@benchmark mode_rate($profile, 5)
@benchmark estimate_max_harmonic(0.1, 0.4)
```

## Common Performance Issues and Fixes

### Issue 1: Auto-estimation creates too many harmonics

**Symptoms**: Model slow, uses lots of memory

**Fix**:
```julia
# Instead of
discretization = HarmonicBasis(:auto)

# Use
discretization = HarmonicBasis(10)  # Explicit cap
```

### Issue 2: General Fermi surface too slow

**Symptoms**: Streaming matrix computation is slow

**Problem**: Custom vF(θ) requires 2D numerical integration

**Fix**:
```julia
# If you have analytic formula
surface = EllipticFermiSurface2D(...)  # Use analytic if possible

# Or provide simpler vF function
surface = GeneralFermiSurface2D(θ -> 1.0, ...)  # Simpler function
```

### Issue 3: Angle grid too slow or memory-heavy

**Symptoms**: Solve takes long time, memory usage high

**Fix**:
```julia
# Reduce angle resolution if possible
AngleGrid(32)   # Instead of 128

# Or switch to harmonic basis if problem allows
HarmonicBasis(:auto)  # For linear regime
```

### Issue 4: Multiband scaling is poor

**Symptoms**: Time increases linearly with number of bands

**Expected behavior**: Correct (O(N_bands) scaling)

**Optimization**:
```julia
# Reduce harmonics for multiband
model = KineticModel2D(...,
    discretization=HarmonicBasis(5),  # Lower for multiband
    bands=[b1, b2, b3]
)
```

## Architecture Decisions

### Why Abstract Dispatch?

Type unions (e.g., `Union{LinearBGKCollision, QuadraticBGKCollision}`) in hot paths prevent compiler specialization.

**Abstract dispatch** allows specialization at compile time:

```julia
# Method 1: LinearBGKCollision path
function compute_collision(c::LinearBGKCollision, ...)
    # Specialized code for linear
end

# Method 2: QuadraticBGKCollision path
function compute_collision(c::QuadraticBGKCollision, ...)
    # Specialized code for quadratic
end

# Called with concrete type → compiler picks right method at compile time
compute_collision(linear_collision, ...)  # Calls method 1
compute_collision(quadratic_collision, ...)  # Calls method 2
```

### Why Dense Streaming Matrices?

Streaming matrices are small (10-50 × 10-50) and used densely, so:
- Dense storage is optimal
- No benefit from sparsity
- BLAS operations are fast on dense

## Tuning Guide

Start with **Example 1** (linear harmonic), then optimize:

1. **Measure**: Use `SolverConfig(log_every=100)` to track iterations
2. **Profile**: Run with smaller mesh if slow
3. **Adjust**: If slow, try reducing harmonics or switching to angle grid if needed
4. **Validate**: Check against reference if available

## References

- Julia performance tips: <https://docs.julialang.org/en/v1/manual/performance-tips/>
- Type stability: <https://docs.julialang.org/en/v1/manual/performance-tips/#Type-stability>
- Profiling: <https://docs.julialang.org/en/v1/stdlib/Profile/>
- BenchmarkTools: <https://github.com/JuliaCI/BenchmarkTools.jl>
