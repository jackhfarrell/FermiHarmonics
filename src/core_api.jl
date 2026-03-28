abstract type AbstractFermiSurface2D end
abstract type AbstractAngularDiscretization2D end
abstract type AbstractStreamingOperator2D end
abstract type AbstractModeRateProfile end
abstract type AbstractCollisionModel2D end
struct MagneticField2D
    omega_c::Float64
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
The surface holds all quasiparticle properties (vF, nu, mass, charge, shape);
the band adds momentum-relaxing and momentum-conserving scattering rates.
"""
struct Band{S<:AbstractFermiSurface2D}
    name    :: Symbol
    surface :: S
    gamma_mr :: Float64
    gamma_mc :: Float64
end

function Band(surface::S; name, gamma_mr::Real, gamma_mc::Real) where {S<:AbstractFermiSurface2D}
    gamma_mr >= 0 || throw(ArgumentError("gamma_mr must be >= 0"))
    gamma_mc >= 0 || throw(ArgumentError("gamma_mc must be >= 0"))
    return Band{S}(Symbol(name), surface, Float64(gamma_mr), Float64(gamma_mc))
end

@inline coerce_band(band::Band) = band

struct Isotropic2DFermiSurface <: AbstractFermiSurface2D
    name::Symbol
    vF::Float64
    nu::Float64
    mass::Float64
    charge::Float64
end

function Isotropic2DFermiSurface(;
    name=:isotropic_2d,
    vF::Real=1.0,
    nu::Real=1.0,
    mass::Real=1.0,
    charge::Real=-1.0,
)
    vF_value = Float64(vF)
    nu_value = Float64(nu)
    mass_value = Float64(mass)
    charge_value = Float64(charge)
    vF_value > 0.0 || throw(ArgumentError("surface vF must be > 0"))
    nu_value > 0.0 || throw(ArgumentError("surface nu must be > 0"))
    mass_value > 0.0 || throw(ArgumentError("surface mass must be > 0"))
    return Isotropic2DFermiSurface(Symbol(name), vF_value, nu_value, mass_value, charge_value)
end

# ------------------------------------------------------------------
# AbstractFermiSurface2D interface — all subtypes must implement these
# ------------------------------------------------------------------
@inline surface_vF(s::Isotropic2DFermiSurface) = s.vF
@inline surface_max_speed(s::Isotropic2DFermiSurface) = s.vF
@inline surface_vF_angle(s::Isotropic2DFermiSurface, ::Float64) = s.vF
@inline surface_density_of_states(s::Isotropic2DFermiSurface) = s.nu
@inline surface_mass(s::Isotropic2DFermiSurface) = s.mass
@inline surface_charge(s::Isotropic2DFermiSurface) = s.charge

# ------------------------------------------------------------------
# EllipticFermiSurface2D — elliptic Fermi surface
#   vF(θ) = vF0 / sqrt(cos²θ + sin²θ / aspect²)
#   aspect = b/a ratio; 1.0 → isotropic; <1 → compressed along y
# ------------------------------------------------------------------
struct EllipticFermiSurface2D <: AbstractFermiSurface2D
    name    :: Symbol
    vF0     :: Float64   # speed at θ=0 (x-axis semi-axis)
    aspect  :: Float64   # b/a
    nu      :: Float64
    mass    :: Float64
    charge  :: Float64
    max_vF  :: Float64   # precomputed CFL bound = vF0 * max(1, 1/aspect)
end

function EllipticFermiSurface2D(;
    name = :elliptic_2d,
    vF0::Real,
    aspect::Real,
    nu::Real,
    mass::Real,
    charge::Real,
)
    vF0_v   = Float64(vF0)
    asp_v   = Float64(aspect)
    nu_v    = Float64(nu)
    mass_v  = Float64(mass)
    vF0_v  > 0.0 || throw(ArgumentError("vF0 must be > 0"))
    asp_v  > 0.0 || throw(ArgumentError("aspect must be > 0"))
    nu_v   > 0.0 || throw(ArgumentError("nu must be > 0"))
    mass_v > 0.0 || throw(ArgumentError("mass must be > 0"))
    return EllipticFermiSurface2D(
        Symbol(name), vF0_v, asp_v, nu_v, mass_v, Float64(charge),
        vF0_v * max(1.0, asp_v),   # max vF at θ=π/2 when aspect>1 (denominator=1/aspect)
    )
end

@inline surface_vF(s::EllipticFermiSurface2D)          = s.vF0
@inline surface_max_speed(s::EllipticFermiSurface2D)    = s.max_vF
@inline surface_vF_angle(s::EllipticFermiSurface2D, θ::Float64) =
    s.vF0 / hypot(cos(θ), sin(θ) / s.aspect)
@inline surface_density_of_states(s::EllipticFermiSurface2D) = s.nu
@inline surface_mass(s::EllipticFermiSurface2D)         = s.mass
@inline surface_charge(s::EllipticFermiSurface2D)       = s.charge

# ------------------------------------------------------------------
# GeneralFermiSurface2D{F} — user-supplied vF(θ) function
#   max_vF is required (upper bound on |vF(θ)| for CFL)
# ------------------------------------------------------------------
struct GeneralFermiSurface2D{F} <: AbstractFermiSurface2D
    name    :: Symbol
    vF_func :: F         # vF_func(θ::Float64)::Float64
    max_vF  :: Float64   # user-supplied CFL bound
    nu      :: Float64
    mass    :: Float64
    charge  :: Float64
end

function GeneralFermiSurface2D(
    vF_func;
    name = :general_2d,
    max_vF::Real,
    nu::Real,
    mass::Real,
    charge::Real,
)
    max_vF_v = Float64(max_vF)
    nu_v     = Float64(nu)
    mass_v   = Float64(mass)
    max_vF_v > 0.0 || throw(ArgumentError("max_vF must be > 0"))
    nu_v     > 0.0 || throw(ArgumentError("nu must be > 0"))
    mass_v   > 0.0 || throw(ArgumentError("mass must be > 0"))
    return GeneralFermiSurface2D{typeof(vF_func)}(
        Symbol(name), vF_func, max_vF_v, nu_v, mass_v, Float64(charge),
    )
end

# surface_vF returns the CFL-relevant maximum for a general surface
@inline surface_vF(s::GeneralFermiSurface2D)          = s.max_vF
@inline surface_max_speed(s::GeneralFermiSurface2D)   = s.max_vF
@inline surface_vF_angle(s::GeneralFermiSurface2D, θ::Float64) = Float64(s.vF_func(θ))
@inline surface_density_of_states(s::GeneralFermiSurface2D) = s.nu
@inline surface_mass(s::GeneralFermiSurface2D)        = s.mass
@inline surface_charge(s::GeneralFermiSurface2D)      = s.charge

struct HarmonicBasis <: AbstractAngularDiscretization2D
    max_harmonic::Union{Int, Symbol, Nothing}
end

function HarmonicBasis(max_harmonic::Union{Integer, Symbol, Nothing})
    if max_harmonic isa Integer
        Int(max_harmonic) >= 1 || throw(ArgumentError("max_harmonic must be >= 1"))
    elseif !(max_harmonic === :auto || isnothing(max_harmonic))
        throw(ArgumentError("HarmonicBasis expects an integer, :auto, or nothing"))
    end
    return HarmonicBasis(max_harmonic isa Integer ? Int(max_harmonic) : max_harmonic)
end

struct AngleGrid <: AbstractAngularDiscretization2D
    theta_count::Int
end

function AngleGrid(theta_count::Integer)
    ntheta = Int(theta_count)
    ntheta >= 8 || throw(ArgumentError("n_angles must be >= 8"))
    iseven(ntheta) || throw(ArgumentError("n_angles must be even"))
    return AngleGrid(ntheta)
end

struct IsotropicHarmonicStreaming <: AbstractStreamingOperator2D end
struct IsotropicAngleStreaming <: AbstractStreamingOperator2D end

struct TwoRateProfile <: AbstractModeRateProfile
    gamma_mc::Float64
end

TwoRateProfile(gamma_mc::Real) = begin
    gamma_value = Float64(gamma_mc)
    gamma_value >= 0.0 || throw(ArgumentError("gamma_mc must be >= 0"))
    TwoRateProfile(gamma_value)
end

struct OddQuarticRateProfile <: AbstractModeRateProfile
    gamma_mc::Float64
    gamma3::Float64
end

function OddQuarticRateProfile(gamma_mc::Real, gamma3::Real=gamma_mc)
    gamma_mc_value = Float64(gamma_mc)
    gamma3_value = Float64(gamma3)
    gamma_mc_value >= 0.0 || throw(ArgumentError("gamma_mc must be >= 0"))
    gamma3_value >= 0.0 || throw(ArgumentError("gamma3 must be >= 0"))
    return OddQuarticRateProfile(gamma_mc_value, gamma3_value)
end

struct ConstantModeRateProfile <: AbstractModeRateProfile
    gamma::Float64
end

function ConstantModeRateProfile(gamma::Real)
    gamma_value = Float64(gamma)
    gamma_value >= 0.0 || throw(ArgumentError("gamma must be >= 0"))
    return ConstantModeRateProfile(gamma_value)
end

struct CustomModeRateProfile{F} <: AbstractModeRateProfile
    rate::F
end

@inline mode_rate(profile::TwoRateProfile, m::Int) = m >= 2 ? profile.gamma_mc : 0.0
@inline mode_rate(profile::OddQuarticRateProfile, m::Int) =
    m < 2 ? 0.0 : (iseven(m) ? profile.gamma_mc : min(profile.gamma_mc, profile.gamma3 * m^4))
@inline mode_rate(profile::ConstantModeRateProfile, m::Int) = m >= 2 ? profile.gamma : 0.0
@inline mode_rate(profile::CustomModeRateProfile, m::Int) = m >= 2 ? Float64(profile.rate(m)) : 0.0

struct LinearBGKCollision{P<:AbstractModeRateProfile} <: AbstractCollisionModel2D
    gamma_mr::Float64
    profile::P
end

function LinearBGKCollision(gamma_mr::Real, profile::AbstractModeRateProfile=TwoRateProfile(0.0))
    gamma_mr_value = Float64(gamma_mr)
    gamma_mr_value >= 0.0 || throw(ArgumentError("gamma_mr must be >= 0"))
    return LinearBGKCollision{typeof(profile)}(gamma_mr_value, profile)
end

struct QuadraticBGKCollision{P<:AbstractModeRateProfile} <: AbstractCollisionModel2D
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
    gamma_mr_value = Float64(gamma_mr)
    mu0_value = Float64(mu0)
    mass_value = Float64(mass)
    theta_oversample_value = Int(theta_oversample)
    gamma_mr_value >= 0.0 || throw(ArgumentError("gamma_mr must be >= 0"))
    mu0_value > 0.0 || throw(ArgumentError("mu0 must be > 0"))
    mass_value > 0.0 || throw(ArgumentError("mass must be > 0"))
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

struct ExactAngleBGKCollision <: AbstractCollisionModel2D
    gamma_mr::Float64
    gamma_mc::Float64
    mu0::Float64
    mass::Float64
    electrostatic_coupling::Float64
end

function ExactAngleBGKCollision(;
    gamma_mr::Real,
    gamma_mc::Real,
    mu0::Real,
    mass::Real,
    electrostatic_coupling::Real=0.0,
)
    gamma_mr_value = Float64(gamma_mr)
    gamma_mc_value = Float64(gamma_mc)
    mu0_value = Float64(mu0)
    mass_value = Float64(mass)
    gamma_mr_value >= 0.0 || throw(ArgumentError("gamma_mr must be >= 0"))
    gamma_mc_value >= 0.0 || throw(ArgumentError("gamma_mc must be >= 0"))
    mu0_value > 0.0 || throw(ArgumentError("mu0 must be > 0"))
    mass_value > 0.0 || throw(ArgumentError("mass must be > 0"))
    return ExactAngleBGKCollision(
        gamma_mr_value,
        gamma_mc_value,
        mu0_value,
        mass_value,
        Float64(electrostatic_coupling),
    )
end

struct TwoRateAngleBGKCollision <: AbstractCollisionModel2D
    gamma_mr::Float64
    gamma_mc::Float64
    mu0::Float64
    mass::Float64
    electrostatic_coupling::Float64
end

function TwoRateAngleBGKCollision(;
    gamma_mr::Real,
    gamma_mc::Real,
    mu0::Real,
    mass::Real,
    electrostatic_coupling::Real=0.0,
)
    gamma_mr_value = Float64(gamma_mr)
    gamma_mc_value = Float64(gamma_mc)
    mu0_value = Float64(mu0)
    mass_value = Float64(mass)
    gamma_mr_value >= 0.0 || throw(ArgumentError("gamma_mr must be >= 0"))
    gamma_mc_value >= 0.0 || throw(ArgumentError("gamma_mc must be >= 0"))
    mu0_value > 0.0 || throw(ArgumentError("mu0 must be > 0"))
    mass_value > 0.0 || throw(ArgumentError("mass must be > 0"))
    return TwoRateAngleBGKCollision(
        gamma_mr_value,
        gamma_mc_value,
        mu0_value,
        mass_value,
        Float64(electrostatic_coupling),
    )
end

struct AngleRateBGKCollision{P<:AbstractModeRateProfile} <: AbstractCollisionModel2D
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
    gamma_mr_value = Float64(gamma_mr)
    mu0_value = Float64(mu0)
    mass_value = Float64(mass)
    filter_tail_value = Int(filter_tail_modes)
    filter_max_value = Float64(filter_max_multiplier)
    gamma_mr_value >= 0.0 || throw(ArgumentError("gamma_mr must be >= 0"))
    mu0_value > 0.0 || throw(ArgumentError("mu0 must be > 0"))
    mass_value > 0.0 || throw(ArgumentError("mass must be > 0"))
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

    if discretization isa HarmonicBasis
        streaming isa IsotropicHarmonicStreaming ||
            throw(ArgumentError("HarmonicBasis currently requires IsotropicHarmonicStreaming"))
        collision isa Union{LinearBGKCollision, QuadraticBGKCollision} ||
            throw(ArgumentError("HarmonicBasis currently supports LinearBGKCollision or QuadraticBGKCollision"))
    elseif discretization isa AngleGrid
        streaming isa IsotropicAngleStreaming ||
            throw(ArgumentError("AngleGrid currently requires IsotropicAngleStreaming"))
        collision isa Union{ExactAngleBGKCollision, TwoRateAngleBGKCollision, AngleRateBGKCollision} ||
            throw(ArgumentError("AngleGrid currently supports ExactAngleBGKCollision, TwoRateAngleBGKCollision, or AngleRateBGKCollision"))
    else
        throw(ArgumentError("unsupported angular discretization $(typeof(discretization))"))
    end

    if !isempty(band_specs)
        collision isa LinearBGKCollision ||
            throw(ArgumentError("multiband support currently requires LinearBGKCollision"))
        discretization isa HarmonicBasis ||
            throw(ArgumentError("multiband support currently requires HarmonicBasis"))
        length(unique(b.name for b in band_specs)) == length(band_specs) ||
            throw(ArgumentError("band names must be unique"))
        # No N-band limit — drag coupling supported for N=2; N≥3 uses mean-field formula
    end

    if !isnothing(magnetic_field)
        discretization isa HarmonicBasis &&
            collision isa LinearBGKCollision &&
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

Base.@kwdef struct MeshBuildConfig
    recombine_all::Bool = true
    algorithm::Int = 8
    save_groups_of_nodes::Bool = true
    output_mode::Symbol = :temporary
    output_dir::Union{Nothing, String} = nothing
    prefix::String = "ek_mesh_"
    gmsh_options::Dict{String, Float64} = Dict{String, Float64}()
end

function validate(config::MeshBuildConfig)
    config.output_mode in (:temporary, :persistent) ||
        throw(ArgumentError("mesh_build.output_mode must be :temporary or :persistent"))
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
    gamma_mc::Real;
    min_harmonic::Integer = 4,
    max_harmonic::Integer = 100,
)::Int
    gamma_mr < 0 && throw(ArgumentError("gamma_mr must be >= 0"))
    gamma_mc < 0 && throw(ArgumentError("gamma_mc must be >= 0"))

    min_h = Int(min_harmonic)
    max_h = Int(max_harmonic)
    min_h >= 1 || throw(ArgumentError("min_harmonic must be >= 1"))
    max_h >= min_h || throw(ArgumentError("max_harmonic must be >= min_harmonic"))

    gamma_total = Float64(gamma_mr) + Float64(gamma_mc)
    gamma_total <= 0 && return max_h
    gamma_total >= AUTO_HARMONIC_GAMMA_HIGH && return min_h

    frac = log1p(gamma_total) / log1p(AUTO_HARMONIC_GAMMA_HIGH)
    estimate = max_h - (max_h - min_h) * frac
    return clamp(ceil(Int, estimate), min_h, max_h)
end

function resolve_max_harmonic(max_harmonic_kw, config::SolverConfig, gamma_mr::Real, gamma_mc::Real)
    if max_harmonic_kw isa Integer
        M = Int(max_harmonic_kw)
        M >= 1 || throw(ArgumentError("max_harmonic must be >= 1"))
        return M, :manual
    end

    if max_harmonic_kw === :auto || isnothing(max_harmonic_kw)
        M = estimate_max_harmonic(
            gamma_mr,
            gamma_mc;
            min_harmonic=config.min_harmonic,
            max_harmonic=config.max_harmonic_auto,
        )
        return M, :auto
    end

    throw(ArgumentError("max_harmonic must be an Integer, :auto, or nothing"))
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
        return (
            u0 = collect(Float64, u0_override),
            mode = :same,
            source_nvars = target_nvars_int,
            target_nvars = target_nvars_int,
        )
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
    return (
        u0 = collect(Float64, u0_override),
        mode = :same,
        source_nvars = Int(target_nvars),
        target_nvars = Int(target_nvars),
    )
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
        return (
            u0 = collect(Float64, u0_override),
            mode = :same,
            source_nvars = target_nvars,
            target_nvars = target_nvars,
        )
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

mutable struct AngleThreadCache{PF, PI}
    spectrum::Vector{ComplexF64}
    scratch_spectrum::Vector{ComplexF64}
    real_buffer::Vector{Float64}
    fft_plan::PF
    ifft_plan::PI
end

struct AngleTransportData{TC<:AngleThreadCache}
    theta_count::Int
    theta::Vector{Float64}
    cos_theta::Vector{Float64}
    sin_theta::Vector{Float64}
    weight::Float64
    thread_caches::Vector{TC}
end

mutable struct NonlinearThreadCache{PF, PI}
    spectrum::Vector{ComplexF64}
    samples::Vector{ComplexF64}
    scratch_samples::Vector{ComplexF64}
    work_samples::Vector{ComplexF64}
    gradx_samples::Vector{ComplexF64}
    grady_samples::Vector{ComplexF64}
    theta_derivative_samples::Vector{ComplexF64}
    real_work::Vector{Float64}
    real_scratch::Vector{Float64}
    fft_plan::PF
    ifft_plan::PI
end

struct NonlinearTransportData{TC<:NonlinearThreadCache}
    theta_count::Int
    theta::Vector{Float64}
    cos_theta::Vector{Float64}
    sin_theta::Vector{Float64}
    thread_caches::Vector{TC}
end

mutable struct NonlinearTimingThreadStats
    flux_ns::UInt64
    flux_calls::UInt64
    bgk_ns::UInt64
    bgk_calls::UInt64
    boundary_ns::UInt64
    boundary_calls::UInt64
    diffuse_ns::UInt64
    diffuse_calls::UInt64
    speed_ns::UInt64
    speed_calls::UInt64
    gradient_ns::UInt64
    gradient_calls::UInt64
end

NonlinearTimingThreadStats() = NonlinearTimingThreadStats(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

const NONLINEAR_TIMING_ENABLED = Ref(false)
const NONLINEAR_TIMING_STATS = Ref(Vector{NonlinearTimingThreadStats}())

@inline function nonlinear_timing_stats()
    stats = NONLINEAR_TIMING_STATS[]
    nthreads = Threads.nthreads()
    if length(stats) != nthreads
        stats = [NonlinearTimingThreadStats() for _ in 1:nthreads]
        NONLINEAR_TIMING_STATS[] = stats
    end
    return stats
end

@inline nonlinear_timing_enabled() = NONLINEAR_TIMING_ENABLED[]

function enable_nonlinear_timing!()
    nonlinear_timing_stats()
    NONLINEAR_TIMING_ENABLED[] = true
    return nothing
end

function disable_nonlinear_timing!()
    NONLINEAR_TIMING_ENABLED[] = false
    return nothing
end

function reset_nonlinear_timing!()
    stats = nonlinear_timing_stats()
    for stat in stats
        stat.flux_ns = 0
        stat.flux_calls = 0
        stat.bgk_ns = 0
        stat.bgk_calls = 0
        stat.boundary_ns = 0
        stat.boundary_calls = 0
        stat.diffuse_ns = 0
        stat.diffuse_calls = 0
        stat.speed_ns = 0
        stat.speed_calls = 0
        stat.gradient_ns = 0
        stat.gradient_calls = 0
    end
    return nothing
end

@inline function record_nonlinear_timing!(category::Symbol, elapsed_ns::UInt64)
    stat = nonlinear_timing_stats()[Threads.threadid()]
    if category === :flux
        stat.flux_ns += elapsed_ns
        stat.flux_calls += 1
    elseif category === :bgk
        stat.bgk_ns += elapsed_ns
        stat.bgk_calls += 1
    elseif category === :boundary
        stat.boundary_ns += elapsed_ns
        stat.boundary_calls += 1
    elseif category === :diffuse
        stat.diffuse_ns += elapsed_ns
        stat.diffuse_calls += 1
    elseif category === :speed
        stat.speed_ns += elapsed_ns
        stat.speed_calls += 1
    elseif category === :gradient
        stat.gradient_ns += elapsed_ns
        stat.gradient_calls += 1
    else
        error("unknown nonlinear timing category: $category")
    end
    return nothing
end

function nonlinear_timing_snapshot()
    totals = Dict(
        :flux => (ns=UInt64(0), calls=UInt64(0)),
        :bgk => (ns=UInt64(0), calls=UInt64(0)),
        :boundary => (ns=UInt64(0), calls=UInt64(0)),
        :diffuse => (ns=UInt64(0), calls=UInt64(0)),
        :speed => (ns=UInt64(0), calls=UInt64(0)),
        :gradient => (ns=UInt64(0), calls=UInt64(0)),
    )
    for stat in nonlinear_timing_stats()
        totals[:flux] = (ns=totals[:flux].ns + stat.flux_ns, calls=totals[:flux].calls + stat.flux_calls)
        totals[:bgk] = (ns=totals[:bgk].ns + stat.bgk_ns, calls=totals[:bgk].calls + stat.bgk_calls)
        totals[:boundary] = (ns=totals[:boundary].ns + stat.boundary_ns, calls=totals[:boundary].calls + stat.boundary_calls)
        totals[:diffuse] = (ns=totals[:diffuse].ns + stat.diffuse_ns, calls=totals[:diffuse].calls + stat.diffuse_calls)
        totals[:speed] = (ns=totals[:speed].ns + stat.speed_ns, calls=totals[:speed].calls + stat.speed_calls)
        totals[:gradient] = (ns=totals[:gradient].ns + stat.gradient_ns, calls=totals[:gradient].calls + stat.gradient_calls)
    end
    return totals
end

function print_nonlinear_timing_summary(io::IO=stdout)
    snapshot = nonlinear_timing_snapshot()
    total_ns = zero(UInt64)
    for key in (:flux, :bgk, :boundary, :diffuse, :speed, :gradient)
        entry = snapshot[key]
        total_ns += entry.ns
        mean_ns = entry.calls == 0 ? 0.0 : Float64(entry.ns) / Float64(entry.calls)
        println(io, rpad(string(key), 10), ": calls=", entry.calls,
                ", total=", round(Float64(entry.ns) * 1.0e-9; digits=4), " s",
                ", mean=", round(mean_ns * 1.0e-6; digits=4), " ms")
    end
    println(io, "total instrumented time: ", round(Float64(total_ns) * 1.0e-9; digits=4), " s")
    return nothing
end

mutable struct BCProjectorCache
    state_buffers::Vector{Vector{Float64}}
    target_buffers::Vector{Vector{Float64}}
    out_buffers::Vector{Vector{Float64}}
    projectors::Dict{Int, SparseMatrixCSC{Float64, Int}}
    nonlinear_faces::Dict{Int, Any}
    initialized::Bool
    nvars::Int
    signature::Tuple{Symbol, Int, Int}
end

BCProjectorCache() = BCProjectorCache(
    [Float64[] for _ in 1:Threads.nthreads()],
    [Float64[] for _ in 1:Threads.nthreads()],
    [Float64[] for _ in 1:Threads.nthreads()],
    Dict{Int, SparseMatrixCSC{Float64, Int}}(),
    Dict{Int, Any}(),
    false,
    0,
    (:unset, 0, 0),
)

struct NonlinearBoundaryFaceData
    unit_normal::SVector{2, Float64}
    incoming_mask::BitVector
    stencil_indices::Matrix{Int}
    stencil_weights::Matrix{Float64}
    projections::Vector{Float64}
    incoming_weight::Float64
    sample_to_harmonics::Union{Nothing, Matrix{Float64}}
end

mutable struct MaxwellWallBC
    p_scatter::Float64
    tol::Float64
    cache::BCProjectorCache
end

MaxwellWallBC(p_scatter::Real; tol::Real=0.0) =
    MaxwellWallBC(Float64(p_scatter), Float64(tol), BCProjectorCache())

mutable struct OhmicContactBC
    p_ohmic_absorb::Float64
    bias::Float64
    tol::Float64
    cache::BCProjectorCache
end

OhmicContactBC(bias::Real; p_ohmic_absorb::Real=1.0, tol::Real=0.0) =
    OhmicContactBC(Float64(p_ohmic_absorb), Float64(bias), Float64(tol), BCProjectorCache())

mutable struct CurrentContactBC
    p_ohmic_absorb::Float64
    target_outward_flux::Float64
    tol::Float64
    cache::BCProjectorCache
end

CurrentContactBC(target_outward_flux::Real; p_ohmic_absorb::Real=1.0, tol::Real=0.0) =
    CurrentContactBC(Float64(p_ohmic_absorb), Float64(target_outward_flux), Float64(tol), BCProjectorCache())

boundary_condition_name(bc) = nameof(typeof(bc))
boundary_condition_name(::MaxwellWallBC) = :maxwell_wall
boundary_condition_name(::OhmicContactBC) = :ohmic_contact
boundary_condition_name(::CurrentContactBC) = :current_contact

@inline transport_symbol(::LinearBGKCollision) = :linear
@inline transport_symbol(::QuadraticBGKCollision) = :parabolic_nonlinear
@inline transport_symbol(::ExactAngleBGKCollision) = :parabolic_nonlinear
@inline transport_symbol(::TwoRateAngleBGKCollision) = :parabolic_nonlinear
@inline transport_symbol(::AngleRateBGKCollision) = :parabolic_nonlinear
@inline transport_symbol(model::KineticModel2D) = transport_symbol(model.collision)

@inline collision_symbol(::LinearBGKCollision) = :linear_mrt
@inline collision_symbol(::QuadraticBGKCollision) = :quadratic_bgk
@inline collision_symbol(::ExactAngleBGKCollision) = :exact_bgk
@inline collision_symbol(::TwoRateAngleBGKCollision) = :two_rate_bgk
@inline collision_symbol(::AngleRateBGKCollision) = :angle_rate_bgk

@inline harmonic_state_nvars(max_harmonic::Integer) = 1 + 2 * Int(max_harmonic)
@inline band_momentum_weight(band::Band) =
    surface_density_of_states(band.surface) * surface_mass(band.surface) * surface_vF(band.surface)
@inline mode_profile(collision::Union{LinearBGKCollision, QuadraticBGKCollision}) = collision.profile
@inline mode_profile(collision::AngleRateBGKCollision) = collision.profile
@inline collision_gamma_mr(collision::Union{LinearBGKCollision, QuadraticBGKCollision, ExactAngleBGKCollision, TwoRateAngleBGKCollision, AngleRateBGKCollision}) = collision.gamma_mr
@inline collision_gamma_mc(collision::Union{ExactAngleBGKCollision, TwoRateAngleBGKCollision}) = collision.gamma_mc
@inline collision_gamma_mc(collision::AngleRateBGKCollision) = profile_reference_rate(collision.profile)
@inline collision_mu0(collision::Union{QuadraticBGKCollision, ExactAngleBGKCollision, TwoRateAngleBGKCollision, AngleRateBGKCollision}) = collision.mu0
@inline collision_mass(collision::Union{QuadraticBGKCollision, ExactAngleBGKCollision, TwoRateAngleBGKCollision, AngleRateBGKCollision}) = collision.mass
@inline collision_electrostatic_coupling(collision::Union{QuadraticBGKCollision, ExactAngleBGKCollision, TwoRateAngleBGKCollision, AngleRateBGKCollision}) = collision.electrostatic_coupling
@inline collision_theta_oversample(collision::QuadraticBGKCollision) = collision.theta_oversample
@inline profile_reference_rate(profile::AbstractModeRateProfile) = mode_rate(profile, 2)

@inline zero_state_speed(mu0::Real, mass::Real) = sqrt(2.0 * Float64(mu0) / Float64(mass))
@inline nonlinear_timestep_speed(vF::Real, chi::Real) = Float64(vF) * (1.0 + abs(Float64(chi)))

function create_angle_transport_data(theta_count::Int)
    ntheta = Int(theta_count)
    ntheta >= 8 || throw(ArgumentError("n_angles must be >= 8"))
    iseven(ntheta) || throw(ArgumentError("n_angles must be even"))
    dtheta = 2.0 * pi / ntheta
    theta = collect(range(0.0, step=dtheta, length=ntheta))
    cos_theta = cos.(theta)
    sin_theta = sin.(theta)
    spectrum = Vector{ComplexF64}(undef, ntheta)
    scratch_spectrum = Vector{ComplexF64}(undef, ntheta)
    real_buffer = Vector{Float64}(undef, ntheta)
    cache_template = AngleThreadCache(
        spectrum,
        scratch_spectrum,
        real_buffer,
        FFTW.plan_fft!(spectrum; flags=FFTW.ESTIMATE),
        FFTW.plan_ifft!(spectrum; flags=FFTW.ESTIMATE),
    )
    TC = typeof(cache_template)
    thread_caches = Vector{TC}(undef, Threads.nthreads())
    thread_caches[1] = cache_template
    for tid in 2:length(thread_caches)
        spectrum_tid = Vector{ComplexF64}(undef, ntheta)
        thread_caches[tid] = AngleThreadCache(
            spectrum_tid,
            Vector{ComplexF64}(undef, ntheta),
            Vector{Float64}(undef, ntheta),
            FFTW.plan_fft!(spectrum_tid; flags=FFTW.ESTIMATE),
            FFTW.plan_ifft!(spectrum_tid; flags=FFTW.ESTIMATE),
        )
    end
    return AngleTransportData(ntheta, theta, cos_theta, sin_theta, dtheta, thread_caches)
end

@inline function nonlinear_theta_count(max_harmonic::Int, theta_oversample::Int)
    base_count = 2 * max_harmonic + 1
    dealiased_count = ceil(Int, 3 * base_count / 2)
    return nextpow(2, max(theta_oversample * base_count, dealiased_count))
end

function create_nonlinear_transport_data(max_harmonic::Int, theta_oversample::Int)
    ntheta = nonlinear_theta_count(max_harmonic, theta_oversample)
    theta = collect(range(0.0, 2.0 * pi, length=ntheta + 1))[1:end-1]
    cos_theta = cos.(theta)
    sin_theta = sin.(theta)
    spectrum = Vector{ComplexF64}(undef, ntheta)
    samples = Vector{ComplexF64}(undef, ntheta)
    scratch_samples = Vector{ComplexF64}(undef, ntheta)
    work_samples = Vector{ComplexF64}(undef, ntheta)
    gradx_samples = Vector{ComplexF64}(undef, ntheta)
    grady_samples = Vector{ComplexF64}(undef, ntheta)
    theta_derivative_samples = Vector{ComplexF64}(undef, ntheta)
    real_work = zeros(Float64, 1 + 2 * max_harmonic)
    real_scratch = zeros(Float64, 1 + 2 * max_harmonic)
    fft_plan = FFTW.plan_fft!(scratch_samples; flags=FFTW.ESTIMATE)
    ifft_plan = FFTW.plan_ifft!(samples; flags=FFTW.ESTIMATE)
    cache_template = NonlinearThreadCache(
        spectrum,
        samples,
        scratch_samples,
        work_samples,
        gradx_samples,
        grady_samples,
        theta_derivative_samples,
        real_work,
        real_scratch,
        fft_plan,
        ifft_plan,
    )
    TC = typeof(cache_template)
    thread_caches = Vector{TC}(undef, Threads.nthreads())
    thread_caches[1] = cache_template
    for tid in 2:length(thread_caches)
        spectrum_tid = Vector{ComplexF64}(undef, ntheta)
        samples_tid = Vector{ComplexF64}(undef, ntheta)
        scratch_samples_tid = Vector{ComplexF64}(undef, ntheta)
        work_samples_tid = Vector{ComplexF64}(undef, ntheta)
        gradx_samples_tid = Vector{ComplexF64}(undef, ntheta)
        grady_samples_tid = Vector{ComplexF64}(undef, ntheta)
        theta_derivative_samples_tid = Vector{ComplexF64}(undef, ntheta)
        real_work_tid = zeros(Float64, 1 + 2 * max_harmonic)
        real_scratch_tid = zeros(Float64, 1 + 2 * max_harmonic)
        thread_caches[tid] = NonlinearThreadCache(
            spectrum_tid,
            samples_tid,
            scratch_samples_tid,
            work_samples_tid,
            gradx_samples_tid,
            grady_samples_tid,
            theta_derivative_samples_tid,
            real_work_tid,
            real_scratch_tid,
            FFTW.plan_fft!(scratch_samples_tid; flags=FFTW.ESTIMATE),
            FFTW.plan_ifft!(samples_tid; flags=FFTW.ESTIMATE),
        )
    end
    return NonlinearTransportData(ntheta, theta, cos_theta, sin_theta, thread_caches)
end

function streaming_matrices(M::Int, vF::Float64=1.0)
    n = 1 + 2 * M
    Ax = zeros(Float64, n, n)
    Ay = zeros(Float64, n, n)

    M >= 1 && begin
        Ax[cosine_index(0), cosine_index(1)] = vF
        Ay[cosine_index(0), sine_index(1)] = vF
    end

    @inbounds @simd for m in 1:M
        Ax[cosine_index(m), cosine_index(m - 1)] = 0.5 * vF
        m + 1 <= M && (Ax[cosine_index(m), cosine_index(m + 1)] = 0.5 * vF)

        m - 1 >= 1 && (Ay[cosine_index(m), sine_index(m - 1)] = -0.5 * vF)
        m + 1 <= M && (Ay[cosine_index(m), sine_index(m + 1)] = 0.5 * vF)

        m - 1 >= 1 && (Ax[sine_index(m), sine_index(m - 1)] = 0.5 * vF)
        m + 1 <= M && (Ax[sine_index(m), sine_index(m + 1)] = 0.5 * vF)

        Ay[sine_index(m), cosine_index(m - 1)] = 0.5 * vF
        m + 1 <= M && (Ay[sine_index(m), cosine_index(m + 1)] = -0.5 * vF)
    end

    return Ax, Ay
end

# ------------------------------------------------------------------
# Surface-dispatch overloads for streaming_matrices
# ------------------------------------------------------------------

# Isotropic: delegate to the analytic tridiagonal formula above
@inline streaming_matrices(M::Int, s::Isotropic2DFermiSurface) = streaming_matrices(M, s.vF)

# Elliptic / General: numerical quadrature via surface_vF_angle
streaming_matrices(M::Int, s::EllipticFermiSurface2D) =
    _anisotropic_streaming_matrices(M, θ -> surface_vF_angle(s, θ))

streaming_matrices(M::Int, s::GeneralFermiSurface2D) =
    _anisotropic_streaming_matrices(M, θ -> surface_vF_angle(s, θ))

# Generic fallback for any future AbstractFermiSurface2D subtype implementing surface_vF_angle
streaming_matrices(M::Int, s::AbstractFermiSurface2D) =
    _anisotropic_streaming_matrices(M, θ -> surface_vF_angle(s, θ))

"""
    _anisotropic_streaming_matrices(M, vF_func)

Compute the harmonic-basis streaming matrices Ax and Ay for an anisotropic Fermi surface
by numerical quadrature of `vF_func(θ)`.

Normalization convention matches `streaming_matrices(M, vF::Float64)`:
  - Basis: φ₁=1, φ_{cosine_index(m)}=cos(mθ), φ_{sine_index(m)}=sin(mθ)  (unnormalized)
  - Weight: 1/π
For constant vF, the result is identical to the analytic tridiagonal formula.
"""
function _anisotropic_streaming_matrices(M::Int, vF_func)
    n = 1 + 2M
    N_quad = max(4n, 128)
    θ    = range(0.0, 2π; length = N_quad + 1)[1:N_quad]
    dθ   = 2π / N_quad
    vF_v = vF_func.(θ)

    # Projection weights (row basis): χ[1]=1, χ[j≥2]=cos(mθ) or sin(mθ)
    # Physical column basis (encoding f = a₀ + 2Σ(aₘcos+bₘsin)):
    #   Φ_phys[1]=1, Φ_phys[j≥2]=2cos(mθ) or 2sin(mθ)
    # Streaming matrix: Ax[i,j] = (1/2π) ∫ vF(θ) cos(θ) χᵢ(θ) Φ_phys,j(θ) dθ
    # This gives the asymmetric (monopole, dipole) coupling Ax[1,2]=vF, Ax[2,1]=vF/2
    # matching streaming_matrices(M, vF::Float64) exactly.
    χ     = zeros(N_quad, n)   # row (test) basis: 1, cos, sin, cos2, sin2, ...
    Φphys = zeros(N_quad, n)   # column (physical) basis: 1, 2cos, 2sin, 2cos2, ...
    χ[:, 1] .= 1.0
    Φphys[:, 1] .= 1.0
    for m in 1:M
        χ[:, cosine_index(m)]     .=  cos.(m .* θ)
        χ[:, sine_index(m)]       .=  sin.(m .* θ)
        Φphys[:, cosine_index(m)] .= 2 .* cos.(m .* θ)
        Φphys[:, sine_index(m)]   .= 2 .* sin.(m .* θ)
    end
    cos_θ = cos.(θ)
    sin_θ = sin.(θ)
    Ax = zeros(n, n)
    Ay = zeros(n, n)
    @inbounds for j in 1:n, i in 1:n
        ax = 0.0
        ay = 0.0
        for k in 1:N_quad
            c = vF_v[k] * χ[k, i] * Φphys[k, j]
            ax += c * cos_θ[k]
            ay += c * sin_θ[k]
        end
        Ax[i, j] = ax * dθ / (2π)
        Ay[i, j] = ay * dθ / (2π)
    end
    return Ax, Ay
end

@inline function harmonics_flux!(
    out::AbstractVector{Float64},
    state::AbstractVector{Float64},
    normal::SVector{2, Float64},
    vF::Float64=1.0,
)
    n_vars = length(state)
    max_harmonic_local = (n_vars - 1) ÷ 2
    normal_x, normal_y = normal
    @inbounds begin
        out[cosine_index(0)] = (max_harmonic_local >= 1) ?
            (normal_x * vF * state[cosine_index(1)] +
             normal_y * vF * state[sine_index(1)]) : 0.0
        for m in 1:max_harmonic_local
            out[cosine_index(m)] =
                normal_x * (0.5 * vF) * state[cosine_index(m - 1)] +
                (m + 1 <= max_harmonic_local ? normal_x * (0.5 * vF) * state[cosine_index(m + 1)] : 0.0) +
                (m - 1 >= 1 ? normal_y * (-0.5 * vF) * state[sine_index(m - 1)] : 0.0) +
                (m + 1 <= max_harmonic_local ? normal_y * (0.5 * vF) * state[sine_index(m + 1)] : 0.0)
            out[sine_index(m)] =
                (m - 1 >= 1 ? normal_x * (0.5 * vF) * state[sine_index(m - 1)] : 0.0) +
                (m + 1 <= max_harmonic_local ? normal_x * (0.5 * vF) * state[sine_index(m + 1)] : 0.0) +
                normal_y * (0.5 * vF) * state[cosine_index(m - 1)] +
                (m + 1 <= max_harmonic_local ? normal_y * (-0.5 * vF) * state[cosine_index(m + 1)] : 0.0)
        end
    end
    return out
end

@inline cosine_index(m::Int) = (m == 0) ? 1 : 2m
@inline sine_index(m::Int) = 2m + 1

@inline nonlinear_bias_scale(model::KineticModel2D) =
    collision_mu0(model.collision) * (1.0 + abs(collision_electrostatic_coupling(model.collision)))

function collision_sources!(args...)
    throw(MethodError(collision_sources!, args))
end

function solve(args...; kwargs...)
    throw(ArgumentError("No backend solve method is available. Load a backend package such as Trixi to activate extensions."))
end

function solve_status(args...; kwargs...)
    throw(ArgumentError("solve_status is provided by backend extensions. Load Trixi to use the current backend implementation."))
end

function save_solution_custom(args...; kwargs...)
    throw(ArgumentError("save_solution_custom is provided by backend extensions. Load Trixi to use the current backend implementation."))
end

function save_for_analysis(args...; kwargs...)
    throw(ArgumentError("save_for_analysis is provided by backend extensions. Load Trixi to use the current backend implementation."))
end

function save_mesh_native_analysis(args...; kwargs...)
    throw(ArgumentError("save_mesh_native_analysis is provided by backend extensions. Load Trixi to use the current backend implementation."))
end

function compute_analysis_grids(args...; kwargs...)
    throw(ArgumentError("compute_analysis_grids is provided by backend extensions. Load Trixi to use the current backend implementation."))
end

function current_norm_variables(args...; kwargs...)
    throw(ArgumentError("current_norm_variables is provided by backend extensions. Load Trixi to use the current backend implementation."))
end

function visualization_callback(args...; kwargs...)
    throw(ArgumentError("visualization_callback is provided by backend extensions. Load Trixi to use the current backend implementation."))
end

function evaluate_solution(args...; kwargs...)
    throw(ArgumentError("evaluate_solution is provided by backend extensions. Load Trixi to use the current backend implementation."))
end

function evaluate_observables(args...; kwargs...)
    throw(ArgumentError("evaluate_observables is provided by backend extensions. Load Trixi to use the current backend implementation."))
end

function create_live_dashboard(args...; kwargs...)
    throw(ArgumentError("Live visualization requires loading GLMakie to activate the Makie extension."))
end

function update_live_dashboard!(args...; kwargs...)
    throw(ArgumentError("Live visualization requires loading GLMakie to activate the Makie extension."))
end

function finalize_live_dashboard!(args...; kwargs...)
    return nothing
end

function live_dashboard_is_open(args...; kwargs...)
    return true
end
