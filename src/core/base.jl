"""
    AbstractFermiSurface2D

Fermi surface representation for 2D kinetic transport.

All subtypes must implement the surface interface:
- `fermi_velocity(surface)` — reference Fermi velocity
- `max_speed(surface)` — maximum velocity (for CFL stability)
- `fermi_velocity_angle(surface, θ)` — velocity at angle θ
- `density_of_states(surface)` — density of states
- `mass(surface)` — effective mass
- `charge(surface)` — carrier charge

# Subtypes
- `AbstractAnalyticSurface` — analytic formula (isotropic, elliptic)
- `AbstractUserDefinedSurface` — user-supplied function

# Example
```julia
# Isotropic Fermi surface (constant velocity)
surface = Isotropic2DFermiSurface(fermi_velocity=1.0, nu=1.0, mass=1.0)

# Elliptic anisotropy
surface = EllipticFermiSurface2D(fermi_velocity0=1.0, aspect=2.0, nu=1.0, mass=1.0)
```
"""
abstract type AbstractFermiSurface2D end

"""
    AbstractAnalyticSurface <: AbstractFermiSurface2D

Fermi surfaces with known analytic formulas.
These can use optimized streaming matrix computation.

# Subtypes
- `Isotropic2DFermiSurface` — constant velocity
- `EllipticFermiSurface2D` — elliptic anisotropy
"""
abstract type AbstractAnalyticSurface <: AbstractFermiSurface2D end

"""
    AbstractUserDefinedSurface <: AbstractFermiSurface2D

Fermi surfaces defined by user-supplied functions.
These require numerical quadrature for streaming matrices.

Subtypes:
- `GeneralFermiSurface2D` — arbitrary fermi_velocity(θ) function
"""
abstract type AbstractUserDefinedSurface <: AbstractFermiSurface2D end
abstract type AbstractAngularDiscretization2D end
abstract type AbstractStreamingOperator2D end
abstract type AbstractModeRateProfile end
abstract type AbstractCollisionModel2D end

"""
    AbstractLinearCollision <: AbstractCollisionModel2D

Linear transport collision models for harmonic expansion basis.
All subtypes have a `gamma_mr` field (momentum-relaxing scattering rate).
Collision behavior: τ⁻¹_total = gamma_mr + gamma_ee(mode).

Subtypes:
- `LinearBGKCollision` — BGK collision with mode-dependent profile
- `LinearCollisionMatrix` — User-supplied collision matrix
"""
abstract type AbstractLinearCollision <: AbstractCollisionModel2D end

"""
    AbstractNonlinearAngleCollision <: AbstractCollisionModel2D

Nonlinear parabolic-band collision models for angle-grid discretization.
All subtypes have: `gamma_mr`, `mu0`, `mass`, `electrostatic_coupling`.

Subtypes:
- `QuadraticBGKCollision` — Quadratic band BGK collision
- `ExactAngleBGKCollision` — Exact angle-dependent BGK
- `TwoRateAngleBGKCollision` — Two-rate angle BGK
- `AngleRateBGKCollision` — Mode-rate filtered angle BGK
"""
abstract type AbstractNonlinearAngleCollision <: AbstractCollisionModel2D end

# ============================================================================
# Error Handling and Validation
# ============================================================================

"""
    PhysicsError <: Exception

Structured error type for physics validation issues.

Fields:
- `problem::String` — what went wrong
- `details::String` — specific information (got what?)
- `suggestion::String` — how to fix it
- `context::String` — reference (where in docs?)
"""
struct PhysicsError <: Exception
    problem::String
    details::String
    suggestion::String
    context::String
end

PhysicsError(problem::String) = PhysicsError(problem, "", "", "")

function Base.showerror(io::IO, e::PhysicsError)
    print(io, "\n╔ PhysicsError: ", e.problem)
    if !isempty(e.details)
        print(io, "\n║\n║ Details: ", e.details)
    end
    if !isempty(e.suggestion)
        print(io, "\n║ Suggestion: ", e.suggestion)
    end
    if !isempty(e.context)
        print(io, "\n║ Reference: ", e.context)
    end
    print(io, "\n╚")
end

struct MagneticField2D
    omega_c::Float64
end

function _require_nonneg(name::AbstractString, value::Real)
    val = Float64(value)
    if !(val >= 0.0)
        throw(PhysicsError(
            "$name must be non-negative",
            "got $value",
            "use a non-negative value",
            "physics constraint"
        ))
    end
    return val
end

function _require_pos(name::AbstractString, value::Real)
    val = Float64(value)
    if !(val > 0.0)
        throw(PhysicsError(
            "$name must be positive",
            "got $value",
            "use a positive value (> 0)",
            "physics constraint"
        ))
    end
    return val
end

function MagneticField2D(omega_c::Real)
    omega_value = Float64(omega_c)
    isfinite(omega_value) || throw(ArgumentError("omega_c must be finite"))
    return MagneticField2D(omega_value)
end

"""
    AbstractTransportMode

Phantom type used as a type parameter of `FermiHarmonics2D` to distinguish
linear harmonic transport from nonlinear parabolic-band transport at the type
level, enabling dispatch without runtime branches.
"""
abstract type AbstractTransportMode end

"""Linear harmonic Boltzmann transport."""
struct LinearTransport <: AbstractTransportMode end

"""Nonlinear parabolic-band transport (quadratic flux, BGK collision)."""
struct NonlinearParabolicTransport <: AbstractTransportMode end

"""
    Band{S<:AbstractFermiSurface2D}

Parameters for one carrier species in a multiband linear transport model.
The surface holds all quasiparticle properties (fermi_velocity, nu, mass, charge, shape);
the band adds momentum-relaxing and momentum-conserving scattering rates.
"""
struct Band{S<:AbstractFermiSurface2D}
    name    :: Symbol
    surface :: S
    gamma_mr :: Float64
    gamma_ee :: Float64
end

function Band(surface::S; name, gamma_mr::Real, gamma_ee::Real) where {S<:AbstractFermiSurface2D}
    return Band{S}(
        Symbol(name),
        surface,
        _require_nonneg("gamma_mr", gamma_mr),
        _require_nonneg("gamma_ee", gamma_ee),
    )
end

@inline coerce_band(band::Band) = band

"""
    Isotropic2DFermiSurface

Isotropic Fermi surface with constant velocity in all directions.

# Parameters
- `fermi_velocity::Float64` — Fermi velocity (must be > 0)
- `nu::Float64` — density of states (must be > 0)
- `mass::Float64` — effective mass (must be > 0)
- `charge::Float64` — carrier charge (default: -1.0 for electrons)

# When to use
- Systems without directional anisotropy (e.g., graphene, conventional metals)
- Testing and development (simplest case)
- As reference for comparing with anisotropic surfaces

# Comparison
- `IsotropicFermiSurface` — constant fermi_velocity (simplest)
- `EllipticFermiSurface2D` — elliptic anisotropy (fermi_velocity depends on angle)
- `GeneralFermiSurface2D` — arbitrary custom function

# Example
```julia
# Electrons in graphene (dimensionless units)
surface = Isotropic2DFermiSurface(
    fermi_velocity = 1.0,           # Fermi velocity
    nu = 1.0,           # density of states
    mass = 1.0,         # effective mass (relative to electron)
    charge = -1.0       # electron charge
)

# Use in model
model = KineticModel2D(surface, HarmonicBasis(:auto), ...)
```
"""
struct Isotropic2DFermiSurface <: AbstractAnalyticSurface
    name::Symbol
    fermi_velocity::Float64
    nu::Float64
    mass::Float64
    charge::Float64
end

function Isotropic2DFermiSurface(;
    name=:isotropic_2d,
    fermi_velocity::Real=1.0,
    nu::Real=1.0,
    mass::Real=1.0,
    charge::Real=-1.0,
)
    fermi_velocity_value = Float64(fermi_velocity)
    nu_value = Float64(nu)
    mass_value = Float64(mass)
    charge_value = Float64(charge)
    fermi_velocity_value > 0.0 || throw(ArgumentError("surface fermi_velocity must be > 0"))
    nu_value > 0.0 || throw(ArgumentError("surface nu must be > 0"))
    mass_value > 0.0 || throw(ArgumentError("surface mass must be > 0"))
    return Isotropic2DFermiSurface(Symbol(name), fermi_velocity_value, nu_value, mass_value, charge_value)
end

# ------------------------------------------------------------------
# AbstractFermiSurface2D interface — all subtypes must implement these
# ------------------------------------------------------------------
@inline fermi_velocity(s::Isotropic2DFermiSurface) = s.fermi_velocity
@inline max_speed(s::Isotropic2DFermiSurface) = s.fermi_velocity
@inline fermi_velocity_angle(s::Isotropic2DFermiSurface, ::Float64) = s.fermi_velocity
@inline density_of_states(s::Isotropic2DFermiSurface) = s.nu
@inline mass(s::Isotropic2DFermiSurface) = s.mass
@inline charge(s::Isotropic2DFermiSurface) = s.charge

# ------------------------------------------------------------------
# EllipticFermiSurface2D — elliptic Fermi surface
#   fermi_velocity(θ) = fermi_velocity0 / sqrt(cos²θ + sin²θ / aspect²)
#   aspect = b/a ratio; 1.0 → isotropic; <1 → compressed along y
# ------------------------------------------------------------------
struct EllipticFermiSurface2D <: AbstractAnalyticSurface
    name    :: Symbol
    fermi_velocity0     :: Float64   # speed at θ=0 (x-axis semi-axis)
    aspect  :: Float64   # b/a
    nu      :: Float64
    mass    :: Float64
    charge  :: Float64
    max_fermi_velocity  :: Float64   # precomputed CFL bound = fermi_velocity0 * max(1, 1/aspect)
end

function EllipticFermiSurface2D(;
    name = :elliptic_2d,
    fermi_velocity0::Real,
    aspect::Real,
    nu::Real,
    mass::Real,
    charge::Real,
)
    fermi_velocity0_value   = Float64(fermi_velocity0)
    asp_v   = Float64(aspect)
    nu_v    = Float64(nu)
    mass_v  = Float64(mass)
    fermi_velocity0_value  > 0.0 || throw(ArgumentError("fermi_velocity0 must be > 0"))
    asp_v  > 0.0 || throw(ArgumentError("aspect must be > 0"))
    nu_v   > 0.0 || throw(ArgumentError("nu must be > 0"))
    mass_v > 0.0 || throw(ArgumentError("mass must be > 0"))
    return EllipticFermiSurface2D(
        Symbol(name), fermi_velocity0_value, asp_v, nu_v, mass_v, Float64(charge),
        fermi_velocity0_value * max(1.0, asp_v),   # max fermi_velocity at θ=π/2 when aspect>1 (denominator=1/aspect)
    )
end

@inline fermi_velocity(s::EllipticFermiSurface2D)          = s.fermi_velocity0
@inline max_speed(s::EllipticFermiSurface2D)    = s.max_fermi_velocity
@inline fermi_velocity_angle(s::EllipticFermiSurface2D, θ::Float64) =
    s.fermi_velocity0 / hypot(cos(θ), sin(θ) / s.aspect)
@inline density_of_states(s::EllipticFermiSurface2D) = s.nu
@inline mass(s::EllipticFermiSurface2D)         = s.mass
@inline charge(s::EllipticFermiSurface2D)       = s.charge

# ------------------------------------------------------------------
# GeneralFermiSurface2D{F} — user-supplied fermi_velocity(θ) function
#   max_fermi_velocity is required (upper bound on |fermi_velocity(θ)| for CFL)
# ------------------------------------------------------------------
struct GeneralFermiSurface2D{F} <: AbstractUserDefinedSurface
    name    :: Symbol
    fermi_velocity_func :: F         # fermi_velocity_func(θ::Float64)::Float64
    max_fermi_velocity  :: Float64   # user-supplied CFL bound
    nu      :: Float64
    mass    :: Float64
    charge  :: Float64
end

function GeneralFermiSurface2D(
    fermi_velocity_func;
    name = :general_2d,
    max_fermi_velocity::Real,
    nu::Real,
    mass::Real,
    charge::Real,
)
    max_fermi_velocity_value = Float64(max_fermi_velocity)
    nu_v     = Float64(nu)
    mass_v   = Float64(mass)
    max_fermi_velocity_value > 0.0 || throw(ArgumentError("max_fermi_velocity must be > 0"))
    nu_v     > 0.0 || throw(ArgumentError("nu must be > 0"))
    mass_v   > 0.0 || throw(ArgumentError("mass must be > 0"))
    return GeneralFermiSurface2D{typeof(fermi_velocity_func)}(
        Symbol(name), fermi_velocity_func, max_fermi_velocity_value, nu_v, mass_v, Float64(charge),
    )
end

# fermi_velocity returns the CFL-relevant maximum for a general surface
@inline fermi_velocity(s::GeneralFermiSurface2D)          = s.max_fermi_velocity
@inline max_speed(s::GeneralFermiSurface2D)   = s.max_fermi_velocity
@inline fermi_velocity_angle(s::GeneralFermiSurface2D, θ::Float64) = Float64(s.fermi_velocity_func(θ))
@inline density_of_states(s::GeneralFermiSurface2D) = s.nu
@inline mass(s::GeneralFermiSurface2D)        = s.mass
@inline charge(s::GeneralFermiSurface2D)      = s.charge

"""
    AbstractHarmonicDiscretization <: AbstractAngularDiscretization2D

Angular discretization using Fourier harmonic basis expansion.
Efficient for linear transport; represents f(θ) ≈ Σ a_m·cos(mθ) + b_m·sin(mθ).

Subtypes:
- `HarmonicBasis` — primary implementation (configurable max harmonic)

Characteristics:
- Memory efficient (1 + 2M variables for M modes)
- Fast linear transport
- No direct angle discretization needed
"""
abstract type AbstractHarmonicDiscretization <: AbstractAngularDiscretization2D end

"""
    AbstractGridDiscretization <: AbstractAngularDiscretization2D

Angular discretization using direct angle grid sampling.
Required for nonlinear transport; evaluates f(θ) at discrete angle points.

Subtypes:
- `AngleGrid` — direct sampling at N equally-spaced angles

Characteristics:
- Direct representation: f_n = f(θ_n)
- Can represent arbitrary (nonlinear) angular dependence
- Requires FFT for harmonic transformation
- Uses more memory (M angles × number of variables)
"""
abstract type AbstractGridDiscretization <: AbstractAngularDiscretization2D end

"""
    HarmonicBasis

Fourier harmonic basis expansion for angular dependence.

# Parameters
- `max_harmonic::Union{Int, Symbol, Nothing}` — Maximum harmonic mode
  - Integer: explicit mode count (1, 2, 3, ...)
  - `:auto`: auto-estimate from collision parameters
  - `nothing`: not yet determined

# When to use
- Linear response regime (small perturbations)
- Efficient representation (memory-proportional to max_harmonic)
- Harmonic basis is fastest for weak damping
- **Recommended default** for most applications

# Comparison
- `HarmonicBasis(:auto)` — Auto-estimate modes (recommended)
- `HarmonicBasis(10)` — Fixed 10 modes
- `AngleGrid(32)` — Direct angle sampling (for nonlinear)

# Example
```julia
# Auto-estimate based on collision rates
basis = HarmonicBasis(:auto)

# Or manually specify
basis = HarmonicBasis(20)  # Use 20 harmonic modes

# In model
model = KineticModel2D(
    surface,
    HarmonicBasis(:auto),
    IsotropicHarmonicStreaming(),
    LinearBGKCollision(0.1)
)
```
"""
struct HarmonicBasis <: AbstractHarmonicDiscretization
    max_harmonic::Union{Int, Symbol, Nothing}
end

function HarmonicBasis(max_harmonic::Union{Integer, Symbol, Nothing})
    if max_harmonic isa Integer
        max_harmonic_int = Int(max_harmonic)
        max_harmonic_int >= 1 || throw(ArgumentError("max_harmonic must be >= 1"))
        return HarmonicBasis(max_harmonic_int)
    elseif max_harmonic === :auto || isnothing(max_harmonic)
        return HarmonicBasis(max_harmonic)
    else
        throw(ArgumentError("HarmonicBasis expects an integer, :auto, or nothing (got $(repr(max_harmonic)))"))
    end
end

"""
    AngleGrid

Direct angle discretization with N equally-spaced samples.

# Parameters
- `theta_count::Int` — Number of angle points (must be even, ≥ 8)

# When to use
- Nonlinear transport regime (essential for nonlinear collision models)
- Need exact angular resolution (not harmonic expansion)
- Working with QuadraticBGKCollision or angle-dependent models

# Comparison
- `HarmonicBasis(:auto)` — Fourier modes (for linear)
- `AngleGrid(32)` — 32 direct angles (for nonlinear)

# Example
```julia
# 32 equally-spaced angles on [0, 2π)
basis = AngleGrid(32)

# In nonlinear model
model = KineticModel2D(
    surface,
    AngleGrid(64),  # 64 angles for high accuracy
    IsotropicAngleStreaming(),
    QuadraticBGKCollision(0.1; mu0=1.0, mass=1.0)
)
```
"""
struct AngleGrid <: AbstractGridDiscretization
    theta_count::Int
end

function AngleGrid(theta_count::Integer)
    ntheta = Int(theta_count)
    if ntheta < 8
        throw(ArgumentError("AngleGrid requires at least 8 angles (got $ntheta)"))
    end
    if !iseven(ntheta)
        throw(ArgumentError("AngleGrid requires even number of angles (got $ntheta)"))
    end
    return AngleGrid(ntheta)
end

"""
    AbstractIsotropicStreaming <: AbstractStreamingOperator2D

Streaming operators for isotropic Fermi surfaces (constant fermi_velocity in all directions).
Can use optimized matrix formulas.

Subtypes:
- `IsotropicHarmonicStreaming` — for harmonic basis
- `IsotropicAngleStreaming` — for angle grids

Future extension point for anisotropic streaming operators.
"""
abstract type AbstractIsotropicStreaming <: AbstractStreamingOperator2D end

struct IsotropicHarmonicStreaming <: AbstractIsotropicStreaming end
struct IsotropicAngleStreaming <: AbstractIsotropicStreaming end

"""
    AbstractBuiltInProfile <: AbstractModeRateProfile

Built-in mode-dependent scattering rate profiles.
These have parameters determined at construction time.

Subtypes:
- `TwoRateProfile` — step function (γ=0 for m<2, γ=γ_ee for m≥2)
- `OddQuarticRateProfile` — quartic rate for odd modes
- `ConstantModeRateProfile` — constant rate for all m≥2
"""
abstract type AbstractBuiltInProfile <: AbstractModeRateProfile end

"""
    AbstractUserProfile <: AbstractModeRateProfile

User-defined mode-dependent scattering rate profiles.

Subtypes:
- `CustomModeRateProfile` — arbitrary rate function mode_rate(m)
"""
abstract type AbstractUserProfile <: AbstractModeRateProfile end

struct TwoRateProfile <: AbstractBuiltInProfile
    gamma_ee::Float64
end

TwoRateProfile(gamma_ee::Real) = begin
    TwoRateProfile(_require_nonneg("gamma_ee", gamma_ee))
end

struct OddQuarticRateProfile <: AbstractBuiltInProfile
    gamma_ee::Float64
    gamma3::Float64
end

function OddQuarticRateProfile(gamma_ee::Real, gamma3::Real=gamma_ee)
    return OddQuarticRateProfile(
        _require_nonneg("gamma_ee", gamma_ee),
        _require_nonneg("gamma3", gamma3),
    )
end

struct ConstantModeRateProfile <: AbstractBuiltInProfile
    gamma::Float64
end

function ConstantModeRateProfile(gamma::Real)
    return ConstantModeRateProfile(_require_nonneg("gamma", gamma))
end

struct CustomModeRateProfile{F} <: AbstractUserProfile
    rate::F
end

@inline mode_rate(profile::TwoRateProfile, m::Int) = m >= 2 ? profile.gamma_ee : 0.0
@inline mode_rate(profile::OddQuarticRateProfile, m::Int) =
    m < 2 ? 0.0 : (iseven(m) ? profile.gamma_ee : min(profile.gamma_ee, profile.gamma3 * m^4))
@inline mode_rate(profile::ConstantModeRateProfile, m::Int) = m >= 2 ? profile.gamma : 0.0
@inline mode_rate(profile::CustomModeRateProfile, m::Int) = m >= 2 ? Float64(profile.rate(m)) : 0.0

"""
    LinearBGKCollision

Linear BGK collision model for harmonic basis transport.

# Description
Relaxation-time approximation (RTA) collision operator with mode-dependent scattering rates:
- Monopole (ρ) and dipole (v) modes: damped at γ_mr
- Higher harmonics: damped at γ_mr + γ_ee(m) where m is harmonic order

# Parameters
- `gamma_mr::Float64` — momentum-relaxing scattering rate
- `profile::AbstractModeRateProfile` — mode-dependent e-e scattering (default: zero)

# When to use
- Linear response regime (small perturbations around equilibrium)
- Harmonic basis discretization (efficient expansion)
- Single-band or multi-band systems
- Need analytic or efficient solutions

# When NOT to use
- Nonlinear transport regime → use QuadraticBGKCollision
- Very low damping with need for exact angle dependence → use ExactAngleBGKCollision

# Example
```julia
# Simple momentum-relaxing only
collision = LinearBGKCollision(0.1)

# With e-e scattering (different for even/odd modes)
collision = LinearBGKCollision(0.1, OddQuarticRateProfile(0.4, 0.1))

# As part of kinetic model
model = KineticModel2D(
    surface,
    HarmonicBasis(:auto),
    IsotropicHarmonicStreaming(),
    collision
)
```
"""
struct LinearBGKCollision{P<:AbstractModeRateProfile} <: AbstractLinearCollision
    gamma_mr::Float64
    profile::P
end

function LinearBGKCollision(gamma_mr::Real, profile::AbstractModeRateProfile=TwoRateProfile(0.0))
    gamma_mr_value = _require_nonneg("gamma_mr", gamma_mr)
    return LinearBGKCollision{typeof(profile)}(gamma_mr_value, profile)
end

struct LinearCollisionMatrix <: AbstractLinearCollision
    gamma_mr::Float64
    gamma_ee::Float64
    gamma3::Float64
    max_harmonic::Int
    matrix::Matrix{Float64}
end

function LinearCollisionMatrix(
    matrix::AbstractMatrix{<:Real};
    gamma_mr::Real=0.0,
    gamma_ee::Real=0.0,
    gamma3::Real=gamma_ee,
)
    size(matrix, 1) == size(matrix, 2) ||
        throw(ArgumentError("collision matrix must be square"))
    n = size(matrix, 1)
    isodd(n) || throw(ArgumentError("collision matrix size must be odd (n = 1 + 2*M)"))
    max_harmonic = (n - 1) ÷ 2
    max_harmonic >= 1 || throw(ArgumentError("collision matrix size must be >= 3"))
    return LinearCollisionMatrix(
        _require_nonneg("gamma_mr", gamma_mr),
        _require_nonneg("gamma_ee", gamma_ee),
        _require_nonneg("gamma3", gamma3),
        Int(max_harmonic),
        Matrix{Float64}(matrix),
    )
end

struct QuadraticBGKCollision{P<:AbstractModeRateProfile} <: AbstractNonlinearAngleCollision
    gamma_mr::Float64
    profile::P
    mu0::Float64
    mass::Float64
    electrostatic_coupling::Float64
    theta_oversample::Int
end

function QuadraticBGKCollision(
    gamma_mr::Real,
    profile::AbstractModeRateProfile=OddQuarticRateProfile(0.0);
    mu0::Real,
    mass::Real,
    electrostatic_coupling::Real=0.0,
    theta_oversample::Integer=1,
)
    gamma_mr_value = _require_nonneg("gamma_mr", gamma_mr)
    mu0_value = _require_pos("mu0", mu0)
    mass_value = _require_pos("mass", mass)
    theta_oversample_value = Int(theta_oversample)
    theta_oversample_value >= 1 || throw(ArgumentError("theta_oversample must be >= 1"))
    return QuadraticBGKCollision{typeof(profile)}(
        gamma_mr_value,
        profile,
        mu0_value,
        mass_value,
        Float64(electrostatic_coupling),
        theta_oversample_value,
    )
end

struct ExactAngleBGKCollision <: AbstractNonlinearAngleCollision
    gamma_mr::Float64
    gamma_ee::Float64
    mu0::Float64
    mass::Float64
    electrostatic_coupling::Float64
end

function ExactAngleBGKCollision(;
    gamma_mr::Real,
    gamma_ee::Real,
    mu0::Real,
    mass::Real,
    electrostatic_coupling::Real=0.0,
)
    gamma_mr_value = _require_nonneg("gamma_mr", gamma_mr)
    gamma_ee_value = _require_nonneg("gamma_ee", gamma_ee)
    mu0_value = _require_pos("mu0", mu0)
    mass_value = _require_pos("mass", mass)
    return ExactAngleBGKCollision(
        gamma_mr_value,
        gamma_ee_value,
        mu0_value,
        mass_value,
        Float64(electrostatic_coupling),
    )
end

struct TwoRateAngleBGKCollision <: AbstractNonlinearAngleCollision
    gamma_mr::Float64
    gamma_ee::Float64
    mu0::Float64
    mass::Float64
    electrostatic_coupling::Float64
end

function TwoRateAngleBGKCollision(;
    gamma_mr::Real,
    gamma_ee::Real,
    mu0::Real,
    mass::Real,
    electrostatic_coupling::Real=0.0,
)
    gamma_mr_value = _require_nonneg("gamma_mr", gamma_mr)
    gamma_ee_value = _require_nonneg("gamma_ee", gamma_ee)
    mu0_value = _require_pos("mu0", mu0)
    mass_value = _require_pos("mass", mass)
    return TwoRateAngleBGKCollision(
        gamma_mr_value,
        gamma_ee_value,
        mu0_value,
        mass_value,
        Float64(electrostatic_coupling),
    )
end

struct AngleRateBGKCollision{P<:AbstractModeRateProfile} <: AbstractNonlinearAngleCollision
    gamma_mr::Float64
    profile::P
    mu0::Float64
    mass::Float64
    electrostatic_coupling::Float64
    filter_tail_modes::Int
    filter_max_multiplier::Float64
    filter_shape::Symbol
end

function AngleRateBGKCollision(
    gamma_mr::Real,
    profile::AbstractModeRateProfile=OddQuarticRateProfile(0.0);
    mu0::Real,
    mass::Real,
    electrostatic_coupling::Real=0.0,
    filter_tail_modes::Integer=10,
    filter_max_multiplier::Real=20.0,
    filter_shape::Symbol=:cosine,
)
    gamma_mr_value = _require_nonneg("gamma_mr", gamma_mr)
    mu0_value = _require_pos("mu0", mu0)
    mass_value = _require_pos("mass", mass)
    filter_tail_value = Int(filter_tail_modes)
    filter_max_value = Float64(filter_max_multiplier)
    filter_tail_value >= 0 || throw(ArgumentError("filter_tail_modes must be >= 0"))
    filter_max_value >= 1.0 || throw(ArgumentError("filter_max_multiplier must be >= 1"))
    filter_shape in (:cosine, :linear, :power) || throw(ArgumentError("filter_shape must be :cosine, :linear, or :power"))
    return AngleRateBGKCollision{typeof(profile)}(
        gamma_mr_value,
        profile,
        mu0_value,
        mass_value,
        Float64(electrostatic_coupling),
        filter_tail_value,
        filter_max_value,
        filter_shape,
    )
end

"""
    KineticModel2D

Kinetic transport model assembling all physics components.

# Parameters
- `surface::AbstractFermiSurface2D` — Fermi surface properties
- `discretization::AbstractAngularDiscretization2D` — Angular basis (HarmonicBasis or AngleGrid)
- `streaming::AbstractStreamingOperator2D` — Streaming operator
- `collision::AbstractCollisionModel2D` — Collision model
- `bands::Vector{Band}` — Multi-band specs (default: single band)
- `gamma_drag::Float64` — Inter-band drag coupling (default: 0.0)
- `magnetic_field::Union{Nothing, MagneticField2D}` — Magnetic field (default: none)

# Validation
The constructor validates compatibility:
- HarmonicBasis requires linear collision + isotropic streaming
- AngleGrid requires nonlinear collision + isotropic streaming
- Multiband requires HarmonicBasis + LinearBGKCollision
- Magnetic field requires single-band linear harmonic model

# Example: Linear Transport
```julia
# Single-band linear transport
model = KineticModel2D(
    surface = Isotropic2DFermiSurface(fermi_velocity=1.0),
    discretization = HarmonicBasis(:auto),
    streaming = IsotropicHarmonicStreaming(),
    collision = LinearBGKCollision(0.1, TwoRateProfile(0.4))
)

sol = solve(problem, model, SolverConfig())
```

# Example: Nonlinear Transport
```julia
# Nonlinear parabolic band transport
model = KineticModel2D(
    surface = Isotropic2DFermiSurface(fermi_velocity=1.0),
    discretization = AngleGrid(32),  # 32 angles
    streaming = IsotropicAngleStreaming(),
    collision = QuadraticBGKCollision(0.1; mu0=1.0, mass=1.0)
)

sol = solve(problem, model, SolverConfig())
```

# Example: Multi-band Transport
```julia
# Two-carrier system (electrons and holes)
bands = [
    Band(Isotropic2DFermiSurface(fermi_velocity=1.0); name=:electrons, gamma_mr=0.1, gamma_ee=0.2),
    Band(Isotropic2DFermiSurface(fermi_velocity=0.8); name=:holes, gamma_mr=0.1, gamma_ee=0.3)
]

model = KineticModel2D(
    surface = Isotropic2DFermiSurface(),  # Reference surface
    discretization = HarmonicBasis(:auto),
    streaming = IsotropicHarmonicStreaming(),
    collision = LinearBGKCollision(0.1),
    bands = bands,
    gamma_drag = 0.05  # Drag coupling
)
```
"""
Base.@kwdef struct KineticModel2D{
    TS<:AbstractFermiSurface2D,
    TD<:AbstractAngularDiscretization2D,
    TStream<:AbstractStreamingOperator2D,
    TCollision<:AbstractCollisionModel2D,
    TBands<:AbstractVector{<:Band},
}
    surface::TS
    discretization::TD
    streaming::TStream
    collision::TCollision
    bands::TBands = Band[]
    gamma_drag::Float64 = 0.0
    magnetic_field::Union{Nothing, MagneticField2D} = nothing
    reference = nothing
end

function build_collision_matrix(
    max_harmonic::Integer;
    gamma_mr::Real,
    gamma_ee::Real,
    gamma3::Real=gamma_ee,
    profile::Union{Nothing, AbstractModeRateProfile}=nothing,
)
    M = Int(max_harmonic)
    M >= 1 || throw(ArgumentError("max_harmonic must be >= 1"))
    gamma_mr_value = _require_nonneg("gamma_mr", gamma_mr)
    gamma_ee_value = _require_nonneg("gamma_ee", gamma_ee)
    gamma3_value = _require_nonneg("gamma3", gamma3)
    rate_profile = isnothing(profile) ? OddQuarticRateProfile(gamma_ee_value, gamma3_value) : profile
    n = 1 + 2 * M
    C = zeros(Float64, n, n)
    @inbounds begin
        C[1, 1] = 0.0
        if n >= 3
            C[cosine_index(1), cosine_index(1)] = -gamma_mr_value
            C[sine_index(1), sine_index(1)] = -gamma_mr_value
        end
        for m in 2:M
            gamma_mode = gamma_mr_value + mode_rate(rate_profile, m)
            C[cosine_index(m), cosine_index(m)] = -gamma_mode
            C[sine_index(m), sine_index(m)] = -gamma_mode
        end
    end
    return C
end

# Validation dispatch for model compatibility
function _validate_model_compatibility(
    discretization::AbstractHarmonicDiscretization,
    streaming::AbstractStreamingOperator2D,
    collision::AbstractCollisionModel2D
)
    if !(streaming isa AbstractIsotropicStreaming)
        throw(PhysicsError(
            "HarmonicBasis requires isotropic streaming",
            "got $(typeof(streaming))",
            "use IsotropicHarmonicStreaming()",
            ""
        ))
    end
    if !(collision isa AbstractLinearCollision)
        throw(PhysicsError(
            "HarmonicBasis incompatible with collision model",
            "got $(typeof(collision))",
            "use LinearBGKCollision, LinearCollisionMatrix, or QuadraticBGKCollision",
            "see @ref LinearBGKCollision"
        ))
    end
end

function _validate_model_compatibility(
    discretization::AbstractGridDiscretization,
    streaming::AbstractStreamingOperator2D,
    collision::AbstractCollisionModel2D
)
    if !(streaming isa AbstractIsotropicStreaming)
        throw(PhysicsError(
            "AngleGrid requires isotropic streaming",
            "got $(typeof(streaming))",
            "use IsotropicAngleStreaming()",
            ""
        ))
    end
    if !(collision isa AbstractNonlinearAngleCollision)
        throw(PhysicsError(
            "AngleGrid requires nonlinear collision model",
            "got $(typeof(collision))",
            "use QuadraticBGKCollision, ExactAngleBGKCollision, TwoRateAngleBGKCollision, or AngleRateBGKCollision",
            "see @ref QuadraticBGKCollision"
        ))
    end
end

function KineticModel2D(
    surface::AbstractFermiSurface2D,
    discretization::AbstractAngularDiscretization2D,
    streaming::AbstractStreamingOperator2D,
    collision::AbstractCollisionModel2D;
    bands=Band[],
    gamma_drag::Real=0.0,
    magnetic_field::Union{Nothing, MagneticField2D}=nothing,
    reference=nothing,
)
    band_specs = Band[coerce_band(band) for band in bands]
    gamma_drag_value = Float64(gamma_drag)
    gamma_drag_value >= 0.0 || throw(ArgumentError("gamma_drag must be >= 0"))

    # Dispatch-based validation replaces explicit type checks
    _validate_model_compatibility(discretization, streaming, collision)

    if !isempty(band_specs)
        collision isa AbstractLinearCollision ||
            throw(ArgumentError("multiband support requires linear collision (got $(typeof(collision)))"))
        discretization isa AbstractHarmonicDiscretization ||
            throw(ArgumentError("multiband support requires harmonic basis discretization"))
        length(unique(b.name for b in band_specs)) == length(band_specs) ||
            throw(ArgumentError("band names must be unique"))
        # No N-band limit — drag coupling supported for N=2; N≥3 uses mean-field formula
    end

    if !isnothing(magnetic_field)
        discretization isa AbstractHarmonicDiscretization &&
            collision isa AbstractLinearCollision &&
            isempty(band_specs) ||
            throw(ArgumentError("magnetic_field currently supported only for single-band linear harmonics"))
    end

    return KineticModel2D(
        surface,
        discretization,
        streaming,
        collision,
        band_specs,
        gamma_drag_value,
        magnetic_field,
        reference,
    )
end

Base.@kwdef struct SolverConfig
    max_harmonic::Int = 60
    min_harmonic::Int = 4
    max_harmonic_auto::Int = 100
    polydeg::Int = 3
    tspan_end::Float64 = 100.0
    residual_tol::Float64 = 1e-5
    residual_reltol::Float64 = 0.0
    cfl::Float64 = 0.8
    log_every::Int = 500
end

function validate(config::SolverConfig)
    config.max_harmonic >= 1 || throw(ArgumentError("config.max_harmonic must be >= 1"))
    config.min_harmonic >= 1 || throw(ArgumentError("config.min_harmonic must be >= 1"))
    config.max_harmonic_auto >= config.min_harmonic ||
        throw(ArgumentError("config.max_harmonic_auto must be >= config.min_harmonic"))
    config.polydeg >= 1 || throw(ArgumentError("config.polydeg must be >= 1"))
    config.tspan_end > 0 || throw(ArgumentError("config.tspan_end must be > 0"))
    config.residual_tol >= 0 || throw(ArgumentError("config.residual_tol must be >= 0"))
    config.residual_reltol >= 0 || throw(ArgumentError("config.residual_reltol must be >= 0"))
    config.residual_tol + config.residual_reltol > 0 ||
        throw(ArgumentError("at least one of residual_tol or residual_reltol must be > 0"))
    config.cfl > 0 || throw(ArgumentError("config.cfl must be > 0"))
    config.log_every >= 1 || throw(ArgumentError("config.log_every must be >= 1"))
    return config
end

"""
    MeshBuildConfig(; ...)

Mesh generation controls when building `.inp` meshes from `.geo` sources.

Defaults:
- `algorithm = 8` (quad-friendly Gmsh algorithm for quasi-structured quads)
- `recombine_all = true` (prefer quadrilateral elements)
- `mesh_scale = 3.0` (coarser meshes; passes through to `Mesh.CharacteristicLengthFactor`)
"""
Base.@kwdef struct MeshBuildConfig
    recombine_all::Bool = true
    algorithm::Int = 8
    save_groups_of_nodes::Bool = true
    mesh_scale::Float64 = 3.0
    output_mode::Symbol = :temporary
    output_dir::Union{Nothing, String} = nothing
    prefix::String = "ek_mesh_"
    gmsh_options::Dict{String, Float64} = Dict{String, Float64}()
end

function validate(config::MeshBuildConfig)
    config.output_mode in (:temporary, :persistent) ||
        throw(ArgumentError("mesh_build.output_mode must be :temporary or :persistent"))
    config.mesh_scale > 0 || throw(ArgumentError("mesh_build.mesh_scale must be > 0"))
    isempty(config.prefix) && throw(ArgumentError("mesh_build.prefix must not be empty"))
    return config
end

struct TrixiProblem
    mesh_path::Union{Nothing, String}
    geometry_path::Union{Nothing, String}
    boundary_conditions::Dict{Symbol, Any}
    mesh_build::MeshBuildConfig
end

function TrixiProblem(;
                      mesh_path::Union{Nothing, AbstractString}=nothing,
                      geometry_path::Union{Nothing, AbstractString}=nothing,
                      boundary_conditions::AbstractDict{Symbol, <:Any},
                      mesh_build::MeshBuildConfig=MeshBuildConfig())
    normalized_mesh_path = isnothing(mesh_path) ? nothing : String(mesh_path)
    normalized_geometry_path = isnothing(geometry_path) ? nothing : String(geometry_path)
    normalized_boundary_conditions = Dict{Symbol, Any}(boundary_conditions)
    problem = TrixiProblem(
        normalized_mesh_path,
        normalized_geometry_path,
        normalized_boundary_conditions,
        mesh_build,
    )
    return validate(problem)
end

function validate(problem::TrixiProblem)
    has_mesh_path = !isnothing(problem.mesh_path)
    has_geometry_path = !isnothing(problem.geometry_path)
    xor(has_mesh_path, has_geometry_path) ||
        throw(ArgumentError("TrixiProblem requires exactly one of mesh_path or geometry_path"))
    validate(problem.mesh_build)
    return problem
end

const AUTO_HARMONIC_GAMMA_HIGH = 300.0

function estimate_max_harmonic(
    gamma_mr::Real,
    gamma_ee::Real;
    min_harmonic::Integer = 4,
    max_harmonic::Integer = 100,
)::Int
    gamma_mr < 0 && throw(ArgumentError("gamma_mr must be >= 0"))
    gamma_ee < 0 && throw(ArgumentError("gamma_ee must be >= 0"))

    min_h = Int(min_harmonic)
    max_h = Int(max_harmonic)
    min_h >= 1 || throw(ArgumentError("min_harmonic must be >= 1"))
    max_h >= min_h || throw(ArgumentError("max_harmonic must be >= min_harmonic"))

    gamma_total = Float64(gamma_mr) + Float64(gamma_ee)
    gamma_total <= 0 && return max_h
    gamma_total >= AUTO_HARMONIC_GAMMA_HIGH && return min_h

    frac = log1p(gamma_total) / log1p(AUTO_HARMONIC_GAMMA_HIGH)
    estimate = max_h - (max_h - min_h) * frac
    return clamp(ceil(Int, estimate), min_h, max_h)
end

function resolve_max_harmonic(max_harmonic_kw, config::SolverConfig, gamma_mr::Real, gamma_ee::Real)
    if max_harmonic_kw isa Integer
        M = Int(max_harmonic_kw)
        M >= 1 || throw(ArgumentError("max_harmonic must be >= 1"))
        return M, :manual
    end

    if max_harmonic_kw === :auto || isnothing(max_harmonic_kw)
        M = estimate_max_harmonic(
            gamma_mr,
            gamma_ee;
            min_harmonic=config.min_harmonic,
            max_harmonic=config.max_harmonic_auto,
        )
        return M, :auto
    end

    throw(ArgumentError("max_harmonic must be an Integer, :auto, or nothing"))
end

@inline function _warm_start_same(u0_override::AbstractVector, target_nvars::Integer)
    return (
        u0 = collect(Float64, u0_override),
        mode = :same,
        source_nvars = Int(target_nvars),
        target_nvars = Int(target_nvars),
    )
end

function resize_warm_start(
    u0_override::AbstractVector,
    target_u0::AbstractVector,
    target_nvars::Integer,
)
    target_nvars_int = Int(target_nvars)
    target_len = length(target_u0)
    target_len % target_nvars_int == 0 ||
        throw(ArgumentError("Target state length $target_len is incompatible with nvars=$target_nvars_int"))
    target_block = target_len ÷ target_nvars_int
    source_len = length(u0_override)

    if source_len == target_len
        return _warm_start_same(u0_override, target_nvars_int)
    end

    if source_len % target_block == 0
        source_nvars = source_len ÷ target_block
        source_state = reshape(collect(Float64, u0_override), source_nvars, target_block)
        target_state = zeros(Float64, target_nvars_int, target_block)
        ncopy = min(source_nvars, target_nvars_int)
        @views target_state[1:ncopy, :] .= source_state[1:ncopy, :]
        mode = source_nvars < target_nvars_int ? :padded : :truncated
        return (
            u0 = vec(target_state),
            mode = mode,
            source_nvars = source_nvars,
            target_nvars = target_nvars_int,
        )
    end

    resized = zeros(Float64, target_len)
    ncopy = min(source_len, target_len)
    @inbounds resized[1:ncopy] .= u0_override[1:ncopy]
    mode = source_len < target_len ? :padded_flat : :truncated_flat
    return (
        u0 = resized,
        mode = mode,
        source_nvars = nothing,
        target_nvars = target_nvars_int,
    )
end

function validate_nonlinear_warm_start(
    u0_override::AbstractVector,
    target_u0::AbstractVector,
    target_nvars::Integer,
)
    source_len = length(u0_override)
    target_len = length(target_u0)
    source_len == target_len || throw(ArgumentError(
        "nonlinear warm start length $source_len does not match target length $target_len for n_angles=$(Int(target_nvars))",
    ))
    return _warm_start_same(u0_override, target_nvars)
end

function resize_multiband_warm_start(
    u0_override::AbstractVector,
    target_u0::AbstractVector,
    target_nbands::Integer,
    target_band_nvars::Integer,
)
    target_nvars = Int(target_nbands) * Int(target_band_nvars)
    target_len = length(target_u0)
    target_len % target_nvars == 0 ||
        throw(ArgumentError("Target state length $target_len is incompatible with nvars=$target_nvars"))
    target_block = target_len ÷ target_nvars
    source_len = length(u0_override)

    if source_len == target_len
        return _warm_start_same(u0_override, target_nvars)
    end

    source_len % target_block == 0 || throw(ArgumentError(
        "warm start length $source_len is incompatible with target discretization for multiband resize",
    ))
    source_nvars = source_len ÷ target_block
    source_nbands = Int(target_nbands)
    source_nvars % source_nbands == 0 || throw(ArgumentError(
        "warm start nvars=$source_nvars is incompatible with source_nbands=$source_nbands",
    ))

    source_band_nvars = source_nvars ÷ source_nbands
    target_band_nvars_int = Int(target_band_nvars)
    source_state = reshape(collect(Float64, u0_override), source_nvars, target_block)
    target_state = zeros(Float64, target_nvars, target_block)

    @inbounds for block_index in 1:target_block, band_index in 1:source_nbands
        source_offset = (band_index - 1) * source_band_nvars
        target_offset = (band_index - 1) * target_band_nvars_int
        ncopy = min(source_band_nvars, target_band_nvars_int)
        target_state[(target_offset + 1):(target_offset + ncopy), block_index] .=
            source_state[(source_offset + 1):(source_offset + ncopy), block_index]
    end

    mode = source_band_nvars < target_band_nvars_int ? :padded_multiband : :truncated_multiband
    return (
        u0 = vec(target_state),
        mode = mode,
        source_nvars = source_nvars,
        target_nvars = target_nvars,
    )
end
