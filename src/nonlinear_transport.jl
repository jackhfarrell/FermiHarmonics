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
@inline transport_is_nonlinear(::FermiAngles2D) = true
@inline nonlinear_data(equations::FermiHarmonics2D) = something(equations.nonlinear_data)
@inline nonlinear_data(equations::FermiAngles2D) = equations.nonlinear_data
@inline nonlinear_collision_is_exact_bgk(equations::FermiHarmonics2D) =
    false
@inline nonlinear_collision_is_exact_bgk(equations::FermiAngles2D) = true
@inline nonlinear_collision_is_quadratic_bgk(equations::FermiHarmonics2D) =
    transport_is_nonlinear(equations) && equations.collision_model === :quadratic_bgk
@inline nonlinear_collision_is_quadratic_bgk(::FermiAngles2D) = false
@inline nonlinear_has_electrostatic_force(equations::FermiHarmonics2D) =
    transport_is_nonlinear(equations) && equations.electrostatic_coupling != 0.0
@inline nonlinear_has_electrostatic_force(equations::FermiAngles2D) =
    equations.electrostatic_coupling != 0.0

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

function validate_collision_model(transport::Symbol, collision_model::Union{Nothing, Symbol})
    default_model = transport === :parabolic_nonlinear ? :quadratic_bgk : :linear_mrt
    model = something(collision_model, default_model)
    if transport === :linear
        model === :linear_mrt ||
            throw(ArgumentError("collision_model must be :linear_mrt for :linear transport"))
        return model
    end

    model in (:quadratic_bgk, :exact_bgk) ||
        throw(ArgumentError("collision_model must be :quadratic_bgk or :exact_bgk for :parabolic_nonlinear transport"))
    return model
end

@inline function zero_state_speed(mu0::Real, mass::Real)
    return sqrt(2.0 * Float64(mu0) / Float64(mass))
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

@inline function parabolic_argument(phi::Real, equations::AbstractFermiTransportEquations2D)
    arg = equations.mu0 + Float64(phi)
    if !(arg > 0.0)
        throw(DomainError(
            Float64(phi),
            "parabolic-band transport requires mu0 + phi > 0 everywhere; got mu0 + phi = $arg",
        ))
    end
    return arg
end

@inline function parabolic_shifted_flux(phi::Real, equations::AbstractFermiTransportEquations2D)
    arg = parabolic_argument(phi, equations)
    prefactor = (2.0 / 3.0) * sqrt(2.0 / equations.mass)
    arg32 = arg * sqrt(arg)
    mu032 = equations.mu0 * sqrt(equations.mu0)
    return prefactor * (arg32 - mu032)
end

@inline function parabolic_shifted_flux_inverse(flux_value::Real, equations::AbstractFermiTransportEquations2D)
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

@inline function parabolic_speed(phi::Real, equations::AbstractFermiTransportEquations2D)
    arg = parabolic_argument(phi, equations)
    return sqrt(2.0 * arg / equations.mass)
end

@inline function parabolic_momentum(phi::Real, equations::AbstractFermiTransportEquations2D)
    arg = parabolic_argument(phi, equations)
    return sqrt(2.0 * equations.mass * arg)
end

@inline quadratic_flux_linear_speed(equations::FermiHarmonics2D) = equations.max_speed

@inline function quadratic_shifted_flux(phi::Real, equations::FermiHarmonics2D)
    vF = quadratic_flux_linear_speed(equations)
    return vF * Float64(phi) + Float64(phi)^2 / (2.0 * equations.mass * vF)
end

@inline function quadratic_speed(phi::Real, equations::FermiHarmonics2D)
    vF = quadratic_flux_linear_speed(equations)
    phi_value = Float64(phi)
    return vF +
           phi_value / (equations.mass * vF) -
           phi_value^2 / (2.0 * equations.mass^2 * vF^3)
end

@inline function quadratic_inverse_momentum(phi::Real, equations::FermiHarmonics2D)
    vF = quadratic_flux_linear_speed(equations)
    phi_value = Float64(phi)
    return 1.0 / (equations.mass * vF) -
           phi_value / (equations.mass^2 * vF^3) +
           1.5 * phi_value^2 / (equations.mass^3 * vF^5)
end

@inline function harmonic_mean_phi(state::AbstractVector{<:Real})
    return 0.5 * Float64(state[1])
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

function harmonic_theta_derivative_to_samples!(
    samples::Vector{ComplexF64},
    state::AbstractVector{<:Real},
    equations::FermiHarmonics2D,
)
    ntheta = nonlinear_data(equations).theta_count
    max_harmonic = (length(state) - 1) ÷ 2
    fill!(samples, 0.0 + 0.0im)
    @inbounds for m in 1:max_harmonic
        coeff = 0.5 * ComplexF64(Float64(state[cosine_index(m)]), -Float64(state[sine_index(m)]))
        derivative_coeff = ComplexF64(-imag(coeff) * m, real(coeff) * m)
        scaled = ntheta * derivative_coeff
        samples[m + 1] = scaled
        samples[ntheta - m + 1] = conj(scaled)
    end
    mul!(samples, get_nonlinear_cache(equations).ifft_plan, samples)
    return samples
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
            directional_cosine * quadratic_shifted_flux(phi, equations),
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
                abs(quadratic_speed(phi, equations))
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
        vmax = max(vmax, abs(quadratic_speed(real(cache.samples[j]), equations)))
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
    return nonlinear_current_from_samples!(cache.scratch_samples, cache.samples, equations)
end

function nonlinear_current_from_samples!(
    scratch_samples::Vector{ComplexF64},
    samples::Vector{ComplexF64},
    equations::FermiHarmonics2D,
)
    @inbounds for j in eachindex(scratch_samples)
        scratch_samples[j] = ComplexF64(quadratic_shifted_flux(real(samples[j]), equations), 0.0)
    end
    mul!(scratch_samples, get_nonlinear_cache(equations).fft_plan, scratch_samples)
    coeff = scratch_samples[2] / nonlinear_data(equations).theta_count
    return (2.0 * real(coeff), -2.0 * imag(coeff))
end

@inline nonlinear_density(state::AbstractVector{<:Real}, equations::FermiHarmonics2D) =
    equations.mass * (equations.mu0 + harmonic_mean_phi(state)) / (2.0 * pi)

@inline nonlinear_current(state::AbstractVector{<:Real}, equations::FermiHarmonics2D) =
    nonlinear_current_components(state, equations)

@inline function nonlinear_density_from_samples(samples::Vector{ComplexF64}, equations::FermiHarmonics2D)
    phi_mean = sum(real, samples) / length(samples)
    return equations.mass * (equations.mu0 + phi_mean) / (2.0 * pi)
end

@inline function quadratic_admissible_drift_speed(mu::Real, equations::FermiHarmonics2D)
    return sqrt(2.0 * Float64(mu) / equations.mass)
end

@inline function quadratic_project_velocity(
    velocity::SVector{2, Float64},
    mu::Real,
    equations::FermiHarmonics2D;
    safety_factor::Float64 = 1.0 - 1.0e-12,
)
    vmax = quadratic_admissible_drift_speed(mu, equations)
    speed = norm(velocity)
    if !(speed > vmax)
        return velocity
    end
    speed > 0.0 || return velocity
    return velocity * (safety_factor * vmax / speed)
end

function recover_mu_u_closed_form(state::AbstractVector{<:Real}, equations::FermiHarmonics2D)
    density = nonlinear_density(state, equations)
    density > 0.0 || throw(DomainError(density, "recover_mu_u requires positive density"))
    mu = 2.0 * pi * density / equations.mass
    delta_mu = mu - equations.mu0
    velocity_scale = equations.mass * equations.max_speed + delta_mu / equations.max_speed
    velocity_scale > 0.0 ||
        throw(DomainError(velocity_scale, "recover_mu_u requires positive drift scale"))
    ux = length(state) >= 2 ? Float64(state[2]) / velocity_scale : 0.0
    uy = length(state) >= 3 ? Float64(state[3]) / velocity_scale : 0.0
    return mu, SVector(ux, uy)
end

function match_quadratic_equilibrium_moments(
    state::AbstractVector{<:Real},
    equations::FermiHarmonics2D,
)
    target_density = nonlinear_density(state, equations)
    target_current = nonlinear_current(state, equations)
    target_density > 0.0 ||
        throw(DomainError(target_density, "quadratic recovery requires positive density"))

    mu_guess, velocity_guess = recover_mu_u_closed_form(state, equations)
    velocity_guess = quadratic_project_velocity(velocity_guess, mu_guess, equations)
    cache = get_nonlinear_cache(equations)
    equilibrium_state = Vector{Float64}(undef, length(state))

    # The quadratic harmonic BGK model recovers macroscopic parameters by matching
    # density and current only; higher harmonics remain nonequilibrium content.
    function residual!(residual, x)
        mu_trial = max(Float64(x[1]), 1.0e-12)
        velocity_trial = quadratic_project_velocity(
            SVector(Float64(x[2]), Float64(x[3])),
            mu_trial,
            equations,
        )
        local_equilibrium_state!(equilibrium_state, mu_trial, velocity_trial, equations)
        density_trial = nonlinear_density(equilibrium_state, equations)
        current_trial = nonlinear_current(equilibrium_state, equations)
        residual[1] = density_trial - target_density
        residual[2] = current_trial[1] - target_current[1]
        residual[3] = current_trial[2] - target_current[2]
        return residual
    end

    initial_guess = [mu_guess, velocity_guess[1], velocity_guess[2]]
    result = nlsolve(
        residual!,
        initial_guess;
        method=:newton,
        ftol=1.0e-12,
        xtol=1.0e-12,
        iterations=50,
    )

    if converged(result)
        mu = max(Float64(result.zero[1]), 1.0e-12)
        velocity = quadratic_project_velocity(
            SVector(Float64(result.zero[2]), Float64(result.zero[3])),
            mu,
            equations,
        )
        return mu, velocity
    end

    return mu_guess, velocity_guess
end

function recover_mu_u(state::AbstractVector{<:Real}, equations::FermiHarmonics2D)
    return match_quadratic_equilibrium_moments(state, equations)
end

function isotropic_equilibrium_state!(
    out::AbstractVector{Float64},
    mu::Real,
    equations::FermiHarmonics2D,
)
    fill!(out, 0.0)
    out[1] = 2.0 * (Float64(mu) - equations.mu0)
    return out
end

function local_equilibrium_samples!(
    samples::Vector{ComplexF64},
    mu::Real,
    velocity::SVector{2, Float64},
    equations::FermiHarmonics2D,
)
    data = nonlinear_data(equations)
    mu_value = Float64(mu)
    delta_mu = mu_value - equations.mu0
    ux, uy = velocity
    u_sq = ux * ux + uy * uy
    vF = equations.max_speed
    linear_drift_scale = equations.mass * vF + delta_mu / vF
    @inbounds for j in eachindex(samples)
        alpha = ux * data.cos_theta[j] + uy * data.sin_theta[j]
        phi_eq = delta_mu + linear_drift_scale * alpha +
                 equations.mass * (alpha^2 - 0.5 * u_sq)
        samples[j] = ComplexF64(phi_eq, 0.0)
    end
    return samples
end

function local_equilibrium_state!(
    out::AbstractVector{Float64},
    mu::Real,
    velocity::SVector{2, Float64},
    equations::FermiHarmonics2D,
)
    cache = get_nonlinear_cache(equations)
    local_equilibrium_samples!(cache.samples, mu, velocity, equations)
    return samples_to_harmonics!(out, cache.samples, equations)
end

function local_equilibrium_state!(
    out::AbstractVector{Float64},
    state::AbstractVector{<:Real},
    equations::FermiHarmonics2D,
)
    mu, velocity = recover_mu_u(state, equations)
    return local_equilibrium_state!(out, mu, velocity, equations)
end

struct ElectrostaticGradientEquation2D{E, N} <: Trixi.AbstractLaplaceDiffusion{2, N}
    diffusivity::Float64
    equations_hyperbolic::E
end

function ElectrostaticGradientEquation2D(equations_hyperbolic::FermiHarmonics2D)
    return ElectrostaticGradientEquation2D{typeof(equations_hyperbolic),
                                           Trixi.nvariables(equations_hyperbolic)}(
        0.0,
        equations_hyperbolic,
    )
end

Trixi.varnames(variable_mapping, equations_parabolic::ElectrostaticGradientEquation2D) =
    Trixi.varnames(variable_mapping, equations_parabolic.equations_hyperbolic)

Trixi.gradient_variable_transformation(::ElectrostaticGradientEquation2D) = Trixi.cons2cons

@inline Trixi.have_constant_diffusivity(::ElectrostaticGradientEquation2D) = Trixi.True()
@inline Trixi.max_diffusivity(::ElectrostaticGradientEquation2D) = 0.0
@inline Trixi.max_diffusivity(u, ::ElectrostaticGradientEquation2D) = 0.0

@inline function Trixi.flux(
    u,
    gradients,
    orientation::Integer,
    equations_parabolic::ElectrostaticGradientEquation2D,
)
    return 0.0 * u
end

@inline function Trixi.penalty(
    u_outer,
    u_inner,
    inv_h,
    equations_parabolic::ElectrostaticGradientEquation2D,
    dg,
)
    return 0.0 * u_inner
end

function electrostatic_force_sources!(
    out::AbstractVector{Float64},
    state::AbstractVector{<:Real},
    gradients,
    equations::FermiHarmonics2D,
)
    if !nonlinear_has_electrostatic_force(equations)
        fill!(out, 0.0)
        return out
    end

    grad_phi0_x = 0.5 * Float64(gradients[1][1])
    grad_phi0_y = 0.5 * Float64(gradients[2][1])
    force_x = -equations.electrostatic_coupling * grad_phi0_x
    force_y = -equations.electrostatic_coupling * grad_phi0_y
    if force_x == 0.0 && force_y == 0.0
        fill!(out, 0.0)
        return out
    end

    cache = get_nonlinear_cache(equations)
    harmonic_state_to_samples!(cache.samples, state, equations)
    harmonic_theta_derivative_to_samples!(cache.scratch_samples, state, equations)
    data = nonlinear_data(equations)

    @inbounds for j in eachindex(cache.samples)
        phi = real(cache.samples[j])
        dphi_dtheta = real(cache.scratch_samples[j])
        cos_theta = data.cos_theta[j]
        sin_theta = data.sin_theta[j]
        force_dot_hat = force_x * cos_theta + force_y * sin_theta
        force_dot_theta = -force_x * sin_theta + force_y * cos_theta
        source = -quadratic_inverse_momentum(phi, equations) * force_dot_theta * dphi_dtheta -
                 quadratic_speed(phi, equations) * force_dot_hat
        cache.scratch_samples[j] = ComplexF64(source, 0.0)
    end

    return samples_to_harmonics!(out, cache.scratch_samples, equations)
end

@inline function analysis_variables(u, equations::FermiHarmonics2D)
    if transport_is_nonlinear(equations)
        density = nonlinear_density(u, equations)
        jx, jy = nonlinear_current(u, equations)
        return SVector(density, jx, jy)
    end

    a1 = length(u) >= 2 ? u[2] : 0.0
    b1 = length(u) >= 3 ? u[3] : 0.0
    return SVector(u[1], a1, b1)
end

function derived_harmonics(
    state::AbstractVector{<:Real},
    equations::FermiHarmonics2D,
)
    a0 = length(state) >= 1 ? Float64(state[1]) : 0.0
    a1 = length(state) >= 2 ? Float64(state[2]) : 0.0
    b1 = length(state) >= 3 ? Float64(state[3]) : 0.0
    return (a0, a1, b1)
end

@inline function get_nonlinear_cache(equations::FermiAngles2D)
    return nonlinear_data(equations).thread_caches[Threads.threadid()]
end

@inline function angle_weighted_mean(
    values::AbstractVector{<:Real},
    equations::FermiAngles2D,
)
    return sum(values) / nonlinear_data(equations).theta_count
end

function derived_harmonics(
    state::AbstractVector{<:Real},
    equations::FermiAngles2D,
)
    data = nonlinear_data(equations)
    inv_ntheta = 1.0 / data.theta_count
    a0 = 0.0
    a1 = 0.0
    b1 = 0.0
    @inbounds for j in eachindex(state)
        phi = Float64(state[j])
        a0 += phi
        a1 += phi * data.cos_theta[j]
        b1 += phi * data.sin_theta[j]
    end
    return (2.0 * a0 * inv_ntheta, 2.0 * a1 * inv_ntheta, 2.0 * b1 * inv_ntheta)
end

function periodic_theta_derivative!(
    out::AbstractVector{Float64},
    state::AbstractVector{<:Real},
    equations::FermiAngles2D,
)
    cache = get_nonlinear_cache(equations)
    ntheta = nonlinear_data(equations).theta_count
    @inbounds for j in eachindex(state)
        cache.spectrum[j] = ComplexF64(Float64(state[j]), 0.0)
    end
    mul!(cache.spectrum, cache.fft_plan, cache.spectrum)
    half_ntheta = ntheta ÷ 2
    @inbounds for mode in 0:(ntheta - 1)
        wave_number = mode <= half_ntheta ? mode : mode - ntheta
        if abs(wave_number) == half_ntheta
            cache.spectrum[mode + 1] = 0.0 + 0.0im
        else
            cache.spectrum[mode + 1] *= ComplexF64(0.0, wave_number)
        end
    end
    mul!(cache.spectrum, cache.ifft_plan, cache.spectrum)
    @inbounds for j in eachindex(out)
        out[j] = real(cache.spectrum[j])
    end
    return out
end

function nonlinear_flux!(
    out::AbstractVector{Float64},
    state::AbstractVector{<:Real},
    normal::SVector{2, Float64},
    equations::FermiAngles2D,
)
    data = nonlinear_data(equations)
    nx, ny = normal
    @inbounds for j in eachindex(out)
        out[j] = (nx * data.cos_theta[j] + ny * data.sin_theta[j]) *
                 parabolic_shifted_flux(state[j], equations)
    end
    return out
end

function nonlinear_max_abs_speed(
    state::AbstractVector{<:Real},
    normal::SVector{2, Float64},
    equations::FermiAngles2D,
)
    data = nonlinear_data(equations)
    nx, ny = normal
    vmax = 0.0
    @inbounds for j in eachindex(state)
        vmax = max(vmax, abs(nx * data.cos_theta[j] + ny * data.sin_theta[j]) *
                         parabolic_speed(state[j], equations))
    end
    return vmax
end

function nonlinear_max_abs_speeds(
    state::AbstractVector{<:Real},
    equations::FermiAngles2D,
)
    vmax = 0.0
    @inbounds for j in eachindex(state)
        vmax = max(vmax, parabolic_speed(state[j], equations))
    end
    return (vmax, vmax)
end

@inline function nonlinear_density(state::AbstractVector{<:Real}, equations::FermiAngles2D)
    return equations.mass * (equations.mu0 + angle_weighted_mean(state, equations)) / (2.0 * pi)
end

function nonlinear_current(state::AbstractVector{<:Real}, equations::FermiAngles2D)
    data = nonlinear_data(equations)
    inv_ntheta = 1.0 / data.theta_count
    jx = 0.0
    jy = 0.0
    @inbounds for j in eachindex(state)
        flux_value = parabolic_shifted_flux(state[j], equations)
        jx += flux_value * data.cos_theta[j]
        jy += flux_value * data.sin_theta[j]
    end
    return (2.0 * jx * inv_ntheta, 2.0 * jy * inv_ntheta)
end

@inline function admissible_drift_speed(mu::Real, equations::FermiAngles2D)
    return sqrt(2.0 * Float64(mu) / equations.mass)
end

@inline function project_admissible_velocity(
    velocity::SVector{2, Float64},
    mu::Real,
    equations::FermiAngles2D;
    safety_factor::Float64 = 1.0 - 1.0e-12,
)
    vmax = admissible_drift_speed(mu, equations)
    speed = norm(velocity)
    if !(speed > vmax)
        return velocity
    end
    speed > 0.0 || return velocity
    return velocity * (safety_factor * vmax / speed)
end

function match_local_equilibrium_moments(
    state::AbstractVector{<:Real},
    equations::FermiAngles2D,
)
    target_density = nonlinear_density(state, equations)
    target_current = nonlinear_current(state, equations)

    target_density > 0.0 || throw(DomainError(target_density, "moment matching requires positive density"))
    mu_guess = 2.0 * pi * target_density / equations.mass
    current_scale = 0.25 * equations.mass * mu_guess
    velocity_guess = current_scale > 0.0 ?
        SVector(target_current[1] / current_scale, target_current[2] / current_scale) :
        SVector(0.0, 0.0)
    velocity_guess = project_admissible_velocity(velocity_guess, mu_guess, equations)
    cache = get_nonlinear_cache(equations)

    function residual!(residual, x)
        mu_trial = max(Float64(x[1]), 1.0e-12)
        velocity_trial = project_admissible_velocity(
            SVector(Float64(x[2]), Float64(x[3])),
            mu_trial,
            equations,
        )
        try
            local_equilibrium_state!(cache.real_buffer, mu_trial, velocity_trial, equations)
            density_trial = nonlinear_density(cache.real_buffer, equations)
            current_trial = nonlinear_current(cache.real_buffer, equations)
            residual[1] = density_trial - target_density
            residual[2] = current_trial[1] - target_current[1]
            residual[3] = current_trial[2] - target_current[2]
        catch err
            if err isa DomainError
                residual[1] = 1.0
                residual[2] = 1.0
                residual[3] = 1.0
            else
                rethrow(err)
            end
        end
        return residual
    end

    initial_guess = [mu_guess, velocity_guess[1], velocity_guess[2]]
    result = nlsolve(residual!, initial_guess; method=:newton, ftol=1.0e-12, xtol=1.0e-12, iterations=50)
    if converged(result)
        mu = max(Float64(result.zero[1]), 1.0e-12)
        velocity = project_admissible_velocity(
            SVector(Float64(result.zero[2]), Float64(result.zero[3])),
            mu,
            equations,
        )
        return mu, velocity
    end

    return mu_guess, velocity_guess
end

function recover_mu_u_closed_form(state::AbstractVector{<:Real}, equations::FermiAngles2D)
    density = nonlinear_density(state, equations)
    density > 0.0 || throw(DomainError(density, "recover_mu_u requires positive density"))
    mu = 2.0 * pi * density / equations.mass
    jx, jy = nonlinear_current(state, equations)
    current_scale = 0.25 * equations.mass * mu
    current_scale > 0.0 || throw(DomainError(current_scale, "recover_mu_u requires positive current scale"))
    return mu, SVector(jx / current_scale, jy / current_scale)
end

function recover_mu_u(state::AbstractVector{<:Real}, equations::FermiAngles2D)
    mu, velocity = recover_mu_u_closed_form(state, equations)
    if norm(velocity) <= admissible_drift_speed(mu, equations)
        return mu, velocity
    end
    return match_local_equilibrium_moments(state, equations)
end

function isotropic_equilibrium_state!(
    out::AbstractVector{Float64},
    mu::Real,
    equations::FermiAngles2D,
)
    fill!(out, Float64(mu) - equations.mu0)
    return out
end

function local_equilibrium_state!(
    out::AbstractVector{Float64},
    mu::Real,
    velocity::SVector{2, Float64},
    equations::FermiAngles2D,
)
    data = nonlinear_data(equations)
    mu_value = Float64(mu)
    ux, uy = velocity
    u_sq = ux * ux + uy * uy
    @inbounds for j in eachindex(out)
        u_dot_hat = ux * data.cos_theta[j] + uy * data.sin_theta[j]
        radicand = 2.0 * equations.mass * mu_value -
                   equations.mass^2 * (u_sq - u_dot_hat^2)
        radicand >= 0.0 || throw(DomainError(
            radicand,
            "drifting local equilibrium is undefined because the Fermi-disk radicand became negative",
        ))
        p_eq = equations.mass * u_dot_hat + sqrt(radicand)
        out[j] = p_eq^2 / (2.0 * equations.mass) - equations.mu0
    end
    return out
end

function local_equilibrium_state!(
    out::AbstractVector{Float64},
    state::AbstractVector{<:Real},
    equations::FermiAngles2D,
)
    mu, velocity = recover_mu_u(state, equations)
    return local_equilibrium_state!(out, mu, velocity, equations)
end

function ElectrostaticGradientEquation2D(equations_hyperbolic::FermiAngles2D)
    return ElectrostaticGradientEquation2D{typeof(equations_hyperbolic),
                                           Trixi.nvariables(equations_hyperbolic)}(
        0.0,
        equations_hyperbolic,
    )
end

function electrostatic_force_sources!(
    out::AbstractVector{Float64},
    state::AbstractVector{<:Real},
    gradients,
    equations::FermiAngles2D,
)
    if !nonlinear_has_electrostatic_force(equations)
        fill!(out, 0.0)
        return out
    end

    grad_phi0_x = 0.5 * derived_harmonics(gradients[1], equations)[1]
    grad_phi0_y = 0.5 * derived_harmonics(gradients[2], equations)[1]
    force_x = -equations.electrostatic_coupling * grad_phi0_x
    force_y = -equations.electrostatic_coupling * grad_phi0_y
    if force_x == 0.0 && force_y == 0.0
        fill!(out, 0.0)
        return out
    end

    cache = get_nonlinear_cache(equations)
    periodic_theta_derivative!(cache.real_buffer, state, equations)
    data = nonlinear_data(equations)

    @inbounds for j in eachindex(out)
        phi = state[j]
        dphi_dtheta = cache.real_buffer[j]
        cos_theta = data.cos_theta[j]
        sin_theta = data.sin_theta[j]
        force_dot_hat = force_x * cos_theta + force_y * sin_theta
        force_dot_theta = -force_x * sin_theta + force_y * cos_theta
        out[j] = -(force_dot_theta / parabolic_momentum(phi, equations)) * dphi_dtheta -
                 parabolic_speed(phi, equations) * force_dot_hat
    end
    return out
end

@inline function analysis_variables(u, equations::FermiAngles2D)
    density = nonlinear_density(u, equations)
    jx, jy = nonlinear_current(u, equations)
    return SVector(density, jx, jy)
end
