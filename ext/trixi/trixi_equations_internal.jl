abstract type AbstractFermiTransportEquations2D{NVARS} <: Trixi.AbstractEquations{2, NVARS} end

struct FermiHarmonics2D{NVARS, TTransport<:AbstractTransportMode, TData, TModel} <: AbstractFermiTransportEquations2D{NVARS}
    gamma_mr::Float64
    gamma_mc::Float64
    gamma3::Float64
    max_speed::Float64
    timestep_speed::Float64
    Ax::Matrix{Float64}
    Ay::Matrix{Float64}
    mu0::Float64
    mass::Float64
    electrostatic_coupling::Float64
    theta_oversample::Int
    nonlinear_data::TData
    model::TModel
end

const LinearFermiHarmonics2D{NVARS, TModel} =
    FermiHarmonics2D{NVARS, LinearTransport, Nothing, TModel}
const NonlinearFermiHarmonics2D{NVARS, TData, TModel} =
    FermiHarmonics2D{NVARS, NonlinearParabolicTransport, TData, TModel}

struct MultiBandFermiHarmonics2D{NVARS, TModel} <: AbstractFermiTransportEquations2D{NVARS}
    bands::Vector{Band}
    gamma_drag::Float64
    max_harmonic::Int
    max_speed::Float64
    timestep_speed::Float64
    Ax::Matrix{Float64}
    Ay::Matrix{Float64}
    model::TModel
end

struct FermiAngles2D{NVARS, TData, TModel} <: AbstractFermiTransportEquations2D{NVARS}
    gamma_mr::Float64
    gamma_mc::Float64
    max_speed::Float64
    timestep_speed::Float64
    mu0::Float64
    mass::Float64
    electrostatic_coupling::Float64
    nonlinear_data::TData
    model::TModel
end

@inline band_count(equations::MultiBandFermiHarmonics2D) = length(equations.bands)
@inline band_nvars(equations::MultiBandFermiHarmonics2D) = harmonic_state_nvars(equations.max_harmonic)
@inline band_offset(equations::MultiBandFermiHarmonics2D, band_index::Integer) =
    (Int(band_index) - 1) * band_nvars(equations)
@inline band_global_cosine_index(equations::MultiBandFermiHarmonics2D, band_index::Integer, m::Int) =
    band_offset(equations, band_index) + cosine_index(m)
@inline band_global_sine_index(equations::MultiBandFermiHarmonics2D, band_index::Integer, m::Int) =
    band_offset(equations, band_index) + sine_index(m)

function block_streaming_matrices(bands::AbstractVector{<:Band}, max_harmonic::Int)
    band_size = harmonic_state_nvars(max_harmonic)
    total_nvars = length(bands) * band_size
    Ax = zeros(Float64, total_nvars, total_nvars)
    Ay = zeros(Float64, total_nvars, total_nvars)

    @inbounds for (band_index, band) in enumerate(bands)
        local_Ax, local_Ay = streaming_matrices(max_harmonic, band.surface)
        offset = (band_index - 1) * band_size
        Ax[(offset + 1):(offset + band_size), (offset + 1):(offset + band_size)] .= local_Ax
        Ay[(offset + 1):(offset + band_size), (offset + 1):(offset + band_size)] .= local_Ay
    end

    return Ax, Ay
end

function harmonic_equations(model::KineticModel2D, max_harmonic::Int)
    collision = model.collision
    surface = model.surface
    collision isa Union{LinearBGKCollision, QuadraticBGKCollision} ||
        throw(ArgumentError("harmonic_equations requires a harmonic BGK collision model"))

    nvars = harmonic_state_nvars(max_harmonic)
    if !isempty(model.bands)
        Ax, Ay = block_streaming_matrices(model.bands, max_harmonic)
        max_speed = maximum(surface_max_speed(band.surface) for band in model.bands)
        return MultiBandFermiHarmonics2D{length(model.bands) * nvars, typeof(model)}(
            model.bands,
            model.gamma_drag,
            max_harmonic,
            max_speed,
            max_speed,
            Ax,
            Ay,
            model,
        )
    end

    vF = collision isa QuadraticBGKCollision ? zero_state_speed(collision.mu0, collision.mass) : surface_vF(surface)
    Ax, Ay = collision isa QuadraticBGKCollision ? streaming_matrices(max_harmonic, vF) : streaming_matrices(max_harmonic, surface)
    gamma_mr = collision_gamma_mr(collision)
    gamma_mc = profile_reference_rate(mode_profile(collision))
    gamma3 = mode_profile(collision) isa OddQuarticRateProfile ? mode_profile(collision).gamma3 : gamma_mc
    mu0_value = collision isa QuadraticBGKCollision ? collision.mu0 : NaN
    mass_value = collision isa QuadraticBGKCollision ? collision.mass : NaN
    chi = collision isa QuadraticBGKCollision ? collision.electrostatic_coupling : 0.0
    theta_oversample = collision isa QuadraticBGKCollision ? collision.theta_oversample : 1
    nonlinear_data = collision isa QuadraticBGKCollision ? create_nonlinear_transport_data(max_harmonic, theta_oversample) : nothing
    transport_type = collision isa QuadraticBGKCollision ? NonlinearParabolicTransport : LinearTransport
    timestep_speed = collision isa QuadraticBGKCollision ? nonlinear_timestep_speed(vF, chi) : vF

    return FermiHarmonics2D{nvars, transport_type, typeof(nonlinear_data), typeof(model)}(
        gamma_mr,
        gamma_mc,
        gamma3,
        vF,
        timestep_speed,
        Ax,
        Ay,
        mu0_value,
        mass_value,
        chi,
        theta_oversample,
        nonlinear_data,
        model,
    )
end

function angle_equations(model::KineticModel2D)
    discretization = model.discretization
    collision = model.collision
    discretization isa AngleGrid || throw(ArgumentError("angle_equations requires an AngleGrid discretization"))
    collision isa Union{ExactAngleBGKCollision, TwoRateAngleBGKCollision} ||
        throw(ArgumentError("angle_equations requires an angle BGK collision model"))

    ntheta = discretization.theta_count
    data = create_angle_transport_data(ntheta)
    vF = zero_state_speed(collision.mu0, collision.mass)
    return FermiAngles2D{ntheta, typeof(data), typeof(model)}(
        collision.gamma_mr,
        collision.gamma_mc,
        vF,
        nonlinear_timestep_speed(vF, collision.electrostatic_coupling),
        collision.mu0,
        collision.mass,
        collision.electrostatic_coupling,
        data,
        model,
    )
end

function build_equations(model::KineticModel2D, config::SolverConfig)
    validate(config)
    if model.discretization isa HarmonicBasis
        if !isempty(model.bands)
            max_harmonic_value = model.discretization.max_harmonic
            if max_harmonic_value isa Integer
                return harmonic_equations(model, max_harmonic_value), :manual
            end
            max_harmonic_resolved = maximum(
                estimate_max_harmonic(
                    band.gamma_mr,
                    band.gamma_mc;
                    min_harmonic=config.min_harmonic,
                    max_harmonic=config.max_harmonic_auto,
                ) for band in model.bands
            )
            return harmonic_equations(model, max_harmonic_resolved), :auto
        end

        collision = model.collision
        gamma_mr = collision_gamma_mr(collision)
        gamma_mc = profile_reference_rate(mode_profile(collision))
        max_harmonic_resolved, harmonic_mode = resolve_max_harmonic(
            model.discretization.max_harmonic,
            config,
            gamma_mr,
            gamma_mc,
        )
        return harmonic_equations(model, max_harmonic_resolved), harmonic_mode
    end

    return angle_equations(model), :angles
end
