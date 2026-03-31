# Extension Guide: Creating Custom Physics Components

This guide explains how to extend ElectronKinetics with custom implementations of physics components. The library is designed for extensibility through abstract dispatch — you can implement your own Fermi surfaces, collision models, and boundary conditions without modifying core code.

## Core Principle: Abstract Dispatch

ElectronKinetics uses **abstract type hierarchies** to enable composition:

```
AbstractFermiSurface2D
├── AbstractAnalyticSurface
│   ├── Isotropic2DFermiSurface
│   ├── EllipticFermiSurface2D
│   └── YourCustomSurface
└── AbstractUserDefinedSurface
    └── GeneralFermiSurface2D

AbstractCollisionModel2D
├── AbstractLinearCollision
│   ├── LinearBGKCollision
│   ├── LinearCollisionMatrix
│   └── YourLinearCollision
└── AbstractNonlinearAngleCollision
    ├── QuadraticBGKCollision
    └── YourNonlinearCollision
```

When you define a custom type and implement its interface, existing code dispatches automatically to your implementation.

## Extension 1: Custom Fermi Surface

### Template: Analytic Surface with Formula

Use this for surfaces with closed-form velocity functions.

```julia
struct MyAnalyticSurface <: AbstractAnalyticSurface
    vF0::Float64                      # Reference velocity
    my_anisotropy_parameter::Float64  # Your physics parameter
    nu::Float64                       # Density of states
    mass::Float64                     # Effective mass
    charge::Float64                   # Carrier charge
end

# Constructor with validation
function MyAnalyticSurface(
    vF0::Real,
    my_param::Real;
    nu::Real=1.0,
    mass::Real=1.0,
    charge::Real=-1.0
)
    vF0 > 0 || throw(PhysicsError("vF0 must be positive", "got $vF0", "", ""))
    my_param > 0 || throw(PhysicsError("my_param must be positive", "got $my_param", "", ""))
    return MyAnalyticSurface(Float64(vF0), Float64(my_param), Float64(nu), Float64(mass), Float64(charge))
end

# Implement the 6 required interface methods
@inline surface_vF(s::MyAnalyticSurface) = s.vF0

@inline surface_max_speed(s::MyAnalyticSurface) = s.vF0  # or compute_max_speed(s) if anisotropic

@inline function surface_vF_angle(s::MyAnalyticSurface, θ::Float64)
    # Return vF(θ) for your surface
    # Example: Trigonal warping
    s.vF0 * (1.0 + s.my_anisotropy_parameter * cos(3θ))
end

@inline surface_density_of_states(s::MyAnalyticSurface) = s.nu

@inline surface_mass(s::MyAnalyticSurface) = s.mass

@inline surface_charge(s::MyAnalyticSurface) = s.charge
```

### Template: User-Defined Function

Use this if you have an arbitrary function vF(θ) or don't know it analytically.

```julia
struct MyGeneralSurface <: AbstractUserDefinedSurface
    vF_func::Function         # Must be: θ::Float64 -> vF::Float64
    max_vF::Float64          # MUST be tight upper bound on |vF(θ)|
    nu::Float64
    mass::Float64
    charge::Float64
end

# Or use the built-in GeneralFermiSurface2D
surface = GeneralFermiSurface2D(
    θ -> 1.0 * (1.0 + 0.2 * cos(4θ)),  # vF(θ) function
    name = :square_modulated,
    max_vF = 1.2,                        # Upper bound (critical for CFL)
    nu = 1.0,
    mass = 1.0,
    charge = -1.0
)
```

### Integration: Using Your Surface

```julia
# Your custom surface works with all discretizations and solvers
model = KineticModel2D(
    MyAnalyticSurface(1.0, 0.3),  # Your custom surface
    HarmonicBasis(:auto),
    IsotropicHarmonicStreaming(),
    LinearBGKCollision(0.1)
)
```

## Extension 2: Custom Collision Model

### Template: Linear Collision (Harmonic Basis)

Use this for collision models compatible with harmonic basis.

```julia
struct MyLinearCollision <: AbstractLinearCollision
    gamma_mr::Float64                   # Momentum-relaxing rate
    my_physical_parameter::Float64      # Your parameter
    # Other fields as needed
end

# Constructor with validation
function MyLinearCollision(
    gamma_mr::Real,
    my_param::Real
)
    gamma_mr >= 0 || throw(PhysicsError("gamma_mr must be >= 0", "got $gamma_mr", "", ""))
    my_param > 0 || throw(PhysicsError("my_param must be positive", "got $my_param", "", ""))
    return MyLinearCollision(Float64(gamma_mr), Float64(my_param))
end

# Required accessor functions
@inline collision_gamma_mr(c::MyLinearCollision) = c.gamma_mr

# If your model uses mode-dependent rates, implement this
@inline function mode_profile(c::MyLinearCollision)
    # Return an AbstractModeRateProfile
    # Option 1: Return existing profile
    TwoRateProfile(0.4)

    # Option 2: Create custom profile (see Extension 3)
    # CustomModeRateProfile(m -> your_rate_function(m))
end
```

### Template: Nonlinear Collision (Angle Grid)

Use this for angle grid discretization with nonlinear effects.

```julia
struct MyNonlinearCollision <: AbstractNonlinearAngleCollision
    gamma_mr::Float64                   # Momentum-relaxing rate
    gamma_ee::Float64                   # Electron-electron rate
    mu0::Float64                        # Band bottom energy
    mass::Float64                       # Band mass
    electrostatic_coupling::Float64     # e-e interaction (chi)
    my_param::Float64                   # Your physical parameter
end

# Constructor with validation
function MyNonlinearCollision(;
    gamma_mr::Real,
    gamma_ee::Real,
    mu0::Real,
    mass::Real,
    electrostatic_coupling::Real=0.0,
    my_param::Real=1.0
)
    # Validate all parameters
    gamma_mr >= 0 || throw(PhysicsError("gamma_mr must be >= 0", "got $gamma_mr", "", ""))
    gamma_ee >= 0 || throw(PhysicsError("gamma_ee must be >= 0", "got $gamma_ee", "", ""))
    mu0 > 0 || throw(PhysicsError("mu0 must be positive", "got $mu0", "", ""))
    mass > 0 || throw(PhysicsError("mass must be positive", "got $mass", "", ""))

    return MyNonlinearCollision(
        Float64(gamma_mr),
        Float64(gamma_ee),
        Float64(mu0),
        Float64(mass),
        Float64(electrostatic_coupling),
        Float64(my_param)
    )
end

# Required accessors
@inline collision_gamma_mr(c::MyNonlinearCollision) = c.gamma_mr
@inline collision_gamma_ee(c::MyNonlinearCollision) = c.gamma_ee
@inline collision_mu0(c::MyNonlinearCollision) = c.mu0
@inline collision_mass(c::MyNonlinearCollision) = c.mass
@inline collision_electrostatic_coupling(c::MyNonlinearCollision) = c.electrostatic_coupling
```

## Extension 3: Custom Mode Rate Profile

Use this for custom mode-dependent scattering.

```julia
struct MyModeRateProfile <: AbstractModeRateProfile
    gamma_ref::Float64
    my_param::Float64
end

# Constructor
function MyModeRateProfile(gamma_ref::Real, my_param::Real)
    gamma_ref >= 0 || throw(PhysicsError("gamma_ref must be >= 0", "got $gamma_ref", "", ""))
    return MyModeRateProfile(Float64(gamma_ref), Float64(my_param))
end

# Implement the rate function
@inline function mode_rate(p::MyModeRateProfile, m::Int)
    if m < 2
        0.0  # Usually: no mode damping for m < 2
    else
        # Example: Quartic scaling
        p.gamma_ref * m^2 * (1.0 + p.my_param * m)
    end
end

# Optional: Reference rate (used by solver for time stepping)
@inline profile_reference_rate(p::MyModeRateProfile) = p.gamma_ref
```

## Extension 4: Custom Boundary Condition

### Template: Wall-Type Boundary

Use for specular/diffuse reflection boundaries.

```julia
mutable struct MyWallBC <: AbstractWallBC
    p_scatter::Float64              # Scattering probability [0, 1]
    my_physical_parameter::Float64  # Your parameter
    tol::Float64                    # Numerical tolerance
    cache::BCProjectorCache         # Thread-safe workspace
end

# Constructor
function MyWallBC(
    p_scatter::Real;
    my_param::Real=1.0,
    tol::Real=0.0
)
    p_scatter >= 0 && p_scatter <= 1 || throw(
        PhysicsError("p_scatter must be in [0, 1]", "got $p_scatter", "", "")
    )
    return MyWallBC(
        Float64(p_scatter),
        Float64(my_param),
        Float64(tol),
        BCProjectorCache()
    )
end
```

### Template: Contact-Type Boundary

Use for voltage-controlled or current-controlled contacts.

```julia
mutable struct MyContactBC <: AbstractContactBC
    p_ohmic_absorb::Float64         # Absorption probability
    my_control_param::Float64       # Voltage or current control
    tol::Float64                    # Numerical tolerance
    cache::BCProjectorCache         # Thread-safe workspace
end

# Constructor
function MyContactBC(
    my_control_param::Real;
    p_ohmic_absorb::Real=1.0,
    tol::Real=0.0
)
    return MyContactBC(
        Float64(p_ohmic_absorb),
        Float64(my_control_param),
        Float64(tol),
        BCProjectorCache()
    )
end
```

## Testing Your Extension

### Basic Assembly Test

```julia
using ElectronKinetics

# Create your custom component
my_surface = MyAnalyticSurface(1.0, 0.3)
my_collision = MyLinearCollision(0.1, 0.4)

# Assemble a model
model = KineticModel2D(
    my_surface,
    HarmonicBasis(10),
    IsotropicHarmonicStreaming(),
    my_collision
)

# Check model properties
println("Model type: $(typeof(model))")
println("Surface vF: $(surface_vF(my_surface))")
println("Collision γ_mr: $(collision_gamma_mr(my_collision))")
```

### Type Stability Test

```julia
using InteractiveUtils

# Check that your surface is type-stable
@code_warntype surface_vF(my_surface)
@code_warntype surface_vF_angle(my_surface, 0.5)

# Check collision properties
@code_warntype collision_gamma_mr(my_collision)
@code_warntype mode_rate(mode_profile(my_collision), 5)
```

### Integration with Solver

Once your component is tested, use it with the solver:

```julia
problem = TrixiProblem(
    mesh_path = "your_mesh.msh",
    boundary_conditions = Dict(:domain => MyWallBC(0.5))
)

config = SolverConfig(max_harmonic=20, tspan_end=100.0)

# Solve with custom components
sol = solve(problem, model, config)
```

## Performance Tips for Extensions

1. **Keep functions simple and allocate-free**
   - Avoid creating arrays in hot-path functions
   - Use scalar operations when possible

2. **Use @inline for small functions**
   ```julia
   @inline surface_vF(s::MyAnalyticSurface) = s.vF0
   ```

3. **Avoid type unions**
   ```julia
   # GOOD: Dispatch on concrete types
   compute(c::LinearBGKCollision) = ...
   compute(c::QuadraticBGKCollision) = ...

   # AVOID: Union types in hot paths
   compute(c::Union{Linear, Quadratic}) = ...
   ```

4. **Pre-validate parameters in constructor**
   ```julia
   function MyCollision(gamma::Real)
       gamma >= 0 || throw(PhysicsError(...))  # Fail fast
       return MyCollision(Float64(gamma))
   end
   ```

## Common Patterns

### Pattern 1: Parameterized Surfaces

```julia
struct MyParametricSurface{N} <: AbstractAnalyticSurface
    params::NTuple{N, Float64}
    nu::Float64
    mass::Float64
    charge::Float64
end

# vF(θ) can depend on tuple of parameters
@inline function surface_vF_angle(s::MyParametricSurface, θ)
    # Use s.params[1], s.params[2], etc.
end
```

### Pattern 2: Tabulated Rates

```julia
struct MyTabularProfile <: AbstractModeRateProfile
    rates::Vector{Float64}  # rates[m] for mode m
end

@inline function mode_rate(p::MyTabularProfile, m::Int)
    (m >= 1 && m <= length(p.rates)) ? p.rates[m] : 0.0
end
```

### Pattern 3: Physics-Dependent Behavior

```julia
struct MyAdaptiveCollision <: AbstractLinearCollision
    gamma_mr::Float64
    use_screening::Bool  # Feature flag
end

@inline function mode_profile(c::MyAdaptiveCollision)
    if c.use_screening
        TwoRateProfile(0.4)
    else
        ConstantModeRateProfile(0.2)
    end
end
```

## Troubleshooting Extensions

### Error: "No method matching..."

**Cause**: Missing required interface function

**Fix**: Implement all required methods for your abstract type. Check the documentation for which methods are required.

### Error: "PhysicsError: X requires Y"

**Cause**: Component incompatibility at model assembly time

**Fix**: Check that your components match:
- HarmonicBasis requires AbstractLinearCollision
- AngleGrid requires AbstractNonlinearAngleCollision

### Type Instability in @code_warntype

**Cause**: Function returns non-concrete type

**Fix**: Ensure all branches return same concrete type:
```julia
# GOOD: Both branches return Float64
@inline function my_rate(p::MyProfile, m::Int)
    m < 2 ? 0.0 : p.gamma * m^2
end

# AVOID: Different types in branches
@inline function my_rate(p::MyProfile, m::Int)
    m < 2 ? nothing : p.gamma * m^2  # Type is Union!
end
```

## References

- Abstract types: <https://docs.julialang.org/en/v1/manual/types/#abstract-types>
- Multiple dispatch: <https://docs.julialang.org/en/v1/manual/methods/>
- Type stability: <https://docs.julialang.org/en/v1/manual/performance-tips/#Type-stability>
- ElectronKinetics API: See module docstring and `?YourType`
