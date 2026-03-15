"""
Utilities for the opt-in nonlinear parabolic-band transport mode.
"""

mutable struct NonlinearThreadCache{PF, PI}
    samples::Vector{ComplexF64}
    scratch_samples::Vector{ComplexF64}
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

@inline transport_is_nonlinear(equations::FermiHarmonics2D) = equations.transport === :parabolic_nonlinear
@inline nonlinear_data(equations::FermiHarmonics2D) = something(equations.nonlinear_data)

function validate_transport_mode(
    transport::Symbol,
    mu0::Union{Nothing, Real},
    mass::Union{Nothing, Real},
    theta_oversample::Integer,
)
    transport in (:linear, :parabolic_nonlinear) ||
        throw(ArgumentError("transport must be :linear or :parabolic_nonlinear"))
    theta_oversample >= 1 || throw(ArgumentError("theta_oversample must be >= 1"))

    if transport === :linear
        return nothing
    end

    isnothing(mu0) && throw(ArgumentError("mu0 is required for :parabolic_nonlinear transport"))
    isnothing(mass) && throw(ArgumentError("mass is required for :parabolic_nonlinear transport"))
    Float64(mu0) > 0 || throw(ArgumentError("mu0 must be > 0 for :parabolic_nonlinear transport"))
    Float64(mass) > 0 || throw(ArgumentError("mass must be > 0 for :parabolic_nonlinear transport"))
    return nothing
end

@inline function zero_state_speed(mu0::Real, mass::Real)
    return sqrt(2.0 * Float64(mu0) / Float64(mass))
end

@inline function nonlinear_theta_count(max_harmonic::Int, theta_oversample::Int)
    return nextpow(2, theta_oversample * (2 * max_harmonic + 1))
end

function create_nonlinear_transport_data(max_harmonic::Int, theta_oversample::Int)
    ntheta = nonlinear_theta_count(max_harmonic, theta_oversample)
    theta = collect(range(0.0, 2.0 * pi, length=ntheta + 1))[1:end-1]
    cos_theta = cos.(theta)
    sin_theta = sin.(theta)
    samples = Vector{ComplexF64}(undef, ntheta)
    scratch_samples = Vector{ComplexF64}(undef, ntheta)
    fft_plan = FFTW.plan_fft!(scratch_samples; flags=FFTW.ESTIMATE)
    ifft_plan = FFTW.plan_ifft!(samples; flags=FFTW.ESTIMATE)
    cache_template = NonlinearThreadCache(samples, scratch_samples, fft_plan, ifft_plan)
    TC = typeof(cache_template)
    thread_caches = Vector{TC}(undef, Threads.nthreads())
    thread_caches[1] = cache_template
    for tid in 2:length(thread_caches)
        samples_tid = Vector{ComplexF64}(undef, ntheta)
        scratch_samples_tid = Vector{ComplexF64}(undef, ntheta)
        thread_caches[tid] = NonlinearThreadCache(
            samples_tid,
            scratch_samples_tid,
            FFTW.plan_fft!(scratch_samples_tid; flags=FFTW.ESTIMATE),
            FFTW.plan_ifft!(samples_tid; flags=FFTW.ESTIMATE),
        )
    end
    return NonlinearTransportData(ntheta, theta, cos_theta, sin_theta, thread_caches)
end

@inline function get_nonlinear_cache(equations::FermiHarmonics2D)
    data = nonlinear_data(equations)
    return data.thread_caches[Threads.threadid()]
end

@inline function parabolic_argument(phi::Real, equations::FermiHarmonics2D)
    arg = equations.mu0 + Float64(phi)
    if !(arg > 0.0)
        throw(DomainError(
            Float64(phi),
            "parabolic-band transport requires mu0 + phi > 0 everywhere; got mu0 + phi = $arg",
        ))
    end
    return arg
end

@inline function parabolic_shifted_flux(phi::Real, equations::FermiHarmonics2D)
    arg = parabolic_argument(phi, equations)
    prefactor = (2.0 / 3.0) * sqrt(2.0 / equations.mass)
    arg32 = arg * sqrt(arg)
    mu032 = equations.mu0 * sqrt(equations.mu0)
    return prefactor * (arg32 - mu032)
end

@inline function parabolic_shifted_flux_inverse(flux_value::Real, equations::FermiHarmonics2D)
    prefactor = (2.0 / 3.0) * sqrt(2.0 / equations.mass)
    mu032 = equations.mu0 * sqrt(equations.mu0)
    base = mu032 + Float64(flux_value) / prefactor
    if !(base > 0.0)
        throw(DomainError(
            Float64(flux_value),
            "shifted parabolic flux inverse requires a positive branch argument; got $base",
        ))
    end
    return base^(2 / 3) - equations.mu0
end

@inline function parabolic_speed(phi::Real, equations::FermiHarmonics2D)
    arg = parabolic_argument(phi, equations)
    return sqrt(2.0 * arg / equations.mass)
end

function harmonic_state_to_samples!(
    samples::Vector{ComplexF64},
    state::AbstractVector{<:Real},
    equations::FermiHarmonics2D,
)
    data = nonlinear_data(equations)
    ntheta = data.theta_count
    max_harmonic = (length(state) - 1) ÷ 2
    fill!(samples, 0.0 + 0.0im)
    samples[1] = ComplexF64(0.5 * ntheta * Float64(state[1]), 0.0)
    @inbounds for m in 1:max_harmonic
        coeff = 0.5 * ComplexF64(Float64(state[cosine_index(m)]), -Float64(state[sine_index(m)]))
        scaled = ntheta * coeff
        samples[m + 1] = scaled
        samples[ntheta - m + 1] = conj(scaled)
    end
    mul!(samples, get_nonlinear_cache(equations).ifft_plan, samples)
    return samples
end

function samples_to_harmonics!(
    out::AbstractVector{Float64},
    samples::Vector{ComplexF64},
    equations::FermiHarmonics2D,
)
    ntheta = nonlinear_data(equations).theta_count
    max_harmonic = (length(out) - 1) ÷ 2
    mul!(samples, get_nonlinear_cache(equations).fft_plan, samples)
    inv_ntheta = 1.0 / ntheta
    @inbounds begin
        out[1] = 2.0 * real(samples[1]) * inv_ntheta
        for m in 1:max_harmonic
            coeff = samples[m + 1] * inv_ntheta
            out[cosine_index(m)] = 2.0 * real(coeff)
            out[sine_index(m)] = -2.0 * imag(coeff)
        end
    end
    return out
end

function nonlinear_flux!(
    out::AbstractVector{Float64},
    state::AbstractVector{<:Real},
    normal::SVector{2, Float64},
    equations::FermiHarmonics2D,
)
    data = nonlinear_data(equations)
    cache = get_nonlinear_cache(equations)
    harmonic_state_to_samples!(cache.samples, state, equations)

    nx, ny = normal
    @inbounds for j in eachindex(cache.scratch_samples)
        phi = real(cache.samples[j])
        directional_cosine = nx * data.cos_theta[j] + ny * data.sin_theta[j]
        cache.scratch_samples[j] = ComplexF64(
            directional_cosine * parabolic_shifted_flux(phi, equations),
            0.0,
        )
    end

    return samples_to_harmonics!(out, cache.scratch_samples, equations)
end

function nonlinear_max_abs_speed(
    state::AbstractVector{<:Real},
    normal::SVector{2, Float64},
    equations::FermiHarmonics2D,
)
    data = nonlinear_data(equations)
    cache = get_nonlinear_cache(equations)
    harmonic_state_to_samples!(cache.samples, state, equations)

    max_speed = 0.0
    nx, ny = normal
    @inbounds for j in eachindex(cache.samples)
        phi = real(cache.samples[j])
        speed = abs(nx * data.cos_theta[j] + ny * data.sin_theta[j]) *
                parabolic_speed(phi, equations)
        max_speed = max(max_speed, speed)
    end
    return max_speed
end

function nonlinear_max_abs_speeds(
    state::AbstractVector{<:Real},
    equations::FermiHarmonics2D,
)
    cache = get_nonlinear_cache(equations)
    harmonic_state_to_samples!(cache.samples, state, equations)

    vmax = 0.0
    @inbounds for j in eachindex(cache.samples)
        vmax = max(vmax, parabolic_speed(real(cache.samples[j]), equations))
    end
    return (vmax, vmax)
end

function nonlinear_current_components(
    state::AbstractVector{<:Real},
    equations::FermiHarmonics2D,
)
    length(state) >= 3 || return (0.0, 0.0)
    cache = get_nonlinear_cache(equations)
    harmonic_state_to_samples!(cache.samples, state, equations)

    @inbounds for j in eachindex(cache.scratch_samples)
        cache.scratch_samples[j] = ComplexF64(parabolic_shifted_flux(real(cache.samples[j]), equations), 0.0)
    end
    mul!(cache.scratch_samples, cache.fft_plan, cache.scratch_samples)
    coeff = cache.scratch_samples[2] / nonlinear_data(equations).theta_count
    return (2.0 * real(coeff), -2.0 * imag(coeff))
end

@inline function analysis_variables(u, equations::FermiHarmonics2D)
    if transport_is_nonlinear(equations)
        jx, jy = nonlinear_current_components(u, equations)
        return SVector(u[1], jx, jy)
    end

    a1 = length(u) >= 2 ? u[2] : 0.0
    b1 = length(u) >= 3 ? u[3] : 0.0
    return SVector(u[1], a1, b1)
end
