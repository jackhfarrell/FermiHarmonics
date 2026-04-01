mutable struct MaxwellWallBC <: AbstractWallBC
    p_scatter::Float64
    tol::Float64
    cache::BCProjectorCache
end

MaxwellWallBC(p_scatter::Real; tol::Real=0.0) =
    MaxwellWallBC(Float64(p_scatter), Float64(tol), BCProjectorCache())

mutable struct OhmicContactBC <: AbstractContactBC
    p_ohmic_absorb::Float64
    bias::Float64
    tol::Float64
    cache::BCProjectorCache
end

OhmicContactBC(bias::Real; p_ohmic_absorb::Real=1.0, tol::Real=0.0) =
    OhmicContactBC(Float64(p_ohmic_absorb), Float64(bias), Float64(tol), BCProjectorCache())

mutable struct CurrentContactBC <: AbstractContactBC
    p_ohmic_absorb::Float64
    target_outward_flux::Float64
    tol::Float64
    cache::BCProjectorCache
end

CurrentContactBC(target_outward_flux::Real; p_ohmic_absorb::Real=1.0, tol::Real=0.0) =
    CurrentContactBC(Float64(p_ohmic_absorb), Float64(target_outward_flux), Float64(tol), BCProjectorCache())

# Boundary condition name dispatch (safe accessors)
function boundary_condition_name(bc::AbstractBoundaryCondition)
    error("$(typeof(bc)) does not implement boundary_condition_name.")
end
boundary_condition_name(::MaxwellWallBC) = :maxwell_wall
boundary_condition_name(::OhmicContactBC) = :ohmic_contact
boundary_condition_name(::CurrentContactBC) = :current_contact

# Scatter probability accessor (wall-type BCs only)
function bc_scatter_probability(bc::AbstractBoundaryCondition)
    error("$(typeof(bc)) does not have scatter probability. Only AbstractWallBC types do.")
end
@inline bc_scatter_probability(bc::AbstractWallBC) = bc.p_scatter

# Absorption probability accessor (contact-type BCs only)
function bc_absorption_probability(bc::AbstractBoundaryCondition)
    error("$(typeof(bc)) does not have absorption probability. Only AbstractContactBC types do.")
end
@inline bc_absorption_probability(bc::AbstractContactBC) = bc.p_ohmic_absorb

# Bias voltage accessor (OhmicContactBC only)
function bc_bias(bc::AbstractBoundaryCondition)
    error("$(typeof(bc)) does not have a bias voltage. Only OhmicContactBC does.")
end
@inline bc_bias(bc::OhmicContactBC) = bc.bias

# Target flux accessor (CurrentContactBC only)
function bc_target_flux(bc::AbstractBoundaryCondition)
    error("$(typeof(bc)) does not have target flux. Only CurrentContactBC does.")
end
@inline bc_target_flux(bc::CurrentContactBC) = bc.target_outward_flux

@inline transport_symbol(::LinearBGKCollision) = :linear
@inline transport_symbol(::LinearCollisionMatrix) = :linear
@inline transport_symbol(::QuadraticBGKCollision) = :parabolic_nonlinear
@inline transport_symbol(::ExactAngleBGKCollision) = :parabolic_nonlinear
@inline transport_symbol(::TwoRateAngleBGKCollision) = :parabolic_nonlinear
@inline transport_symbol(::AngleRateBGKCollision) = :parabolic_nonlinear
@inline transport_symbol(model::KineticModel2D) = transport_symbol(model.collision)

@inline collision_symbol(::LinearBGKCollision) = :linear_mrt
@inline collision_symbol(::LinearCollisionMatrix) = :linear_matrix
@inline collision_symbol(::QuadraticBGKCollision) = :quadratic_bgk
@inline collision_symbol(::ExactAngleBGKCollision) = :exact_bgk
@inline collision_symbol(::TwoRateAngleBGKCollision) = :two_rate_bgk
@inline collision_symbol(::AngleRateBGKCollision) = :angle_rate_bgk

@inline harmonic_state_nvars(max_harmonic::Integer) = 1 + 2 * Int(max_harmonic)
@inline band_momentum_weight(band::Band) =
    surface_density_of_states(band.surface) * surface_mass(band.surface) * surface_vF(band.surface)
# Mode rate profile accessor (works with LinearBGKCollision, QuadraticBGKCollision, AngleRateBGKCollision)
@inline function mode_profile(collision::AbstractCollisionModel2D)
    error("$(typeof(collision)) does not have a mode_profile. Only collisions with AbstractModeRateProfile have this.")
end
@inline mode_profile(c::LinearBGKCollision) = c.profile
@inline mode_profile(c::QuadraticBGKCollision) = c.profile
@inline mode_profile(c::AngleRateBGKCollision) = c.profile

# Gamma_mr accessor — available for all collision types
@inline function collision_gamma_mr(collision::AbstractCollisionModel2D)
    error("$(typeof(collision)) does not implement collision_gamma_mr.")
end
@inline collision_gamma_mr(c::AbstractLinearCollision) = c.gamma_mr
@inline collision_gamma_mr(c::AbstractNonlinearAngleCollision) = c.gamma_mr

# Gamma_ee accessor — varies by collision type
@inline function collision_gamma_ee(collision::AbstractCollisionModel2D)
    error("gamma_ee query is complex and type-dependent. Use mode_rate(c.profile, m) for profiles.")
end
@inline collision_gamma_ee(c::ExactAngleBGKCollision) = c.gamma_ee
@inline collision_gamma_ee(c::TwoRateAngleBGKCollision) = c.gamma_ee
@inline collision_gamma_ee(c::LinearCollisionMatrix) = c.gamma_ee
@inline collision_gamma_ee(c::AngleRateBGKCollision) = profile_reference_rate(c.profile)

# mu0 accessor — for nonlinear angle collisions only
@inline function collision_mu0(collision::AbstractCollisionModel2D)
    error("$(typeof(collision)) does not have mu0. Only AbstractNonlinearAngleCollision models do.")
end
@inline collision_mu0(c::AbstractNonlinearAngleCollision) = c.mu0

# mass accessor — for nonlinear angle collisions only
@inline function collision_mass(collision::AbstractCollisionModel2D)
    error("$(typeof(collision)) does not have mass. Only AbstractNonlinearAngleCollision models do.")
end
@inline collision_mass(c::AbstractNonlinearAngleCollision) = c.mass

# electrostatic_coupling accessor — for nonlinear angle collisions only
@inline function collision_electrostatic_coupling(collision::AbstractCollisionModel2D)
    error("$(typeof(collision)) does not have electrostatic_coupling. Only AbstractNonlinearAngleCollision models do.")
end
@inline collision_electrostatic_coupling(c::AbstractNonlinearAngleCollision) = c.electrostatic_coupling

# theta_oversample accessor — QuadraticBGKCollision only
@inline function collision_theta_oversample(collision::AbstractCollisionModel2D)
    error("$(typeof(collision)) does not have theta_oversample. Only QuadraticBGKCollision has this.")
end
@inline collision_theta_oversample(c::QuadraticBGKCollision) = c.theta_oversample
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

# Optimized path for analytic surfaces: use direct formula (zero quadrature overhead)
@inline streaming_matrices(M::Int, s::Isotropic2DFermiSurface) = streaming_matrices(M, s.vF)

function streaming_matrices(M::Int, s::AbstractAnalyticSurface)
    # All analytic surfaces can use the numerical quadrature path that's optimized
    # for smooth vF_angle functions. For isotropic, the dispatch above will be used.
    _anisotropic_streaming_matrices(M, θ -> surface_vF_angle(s, θ))
end

# General path for user-defined surfaces: use numerical quadrature
function streaming_matrices(M::Int, s::AbstractUserDefinedSurface)
    _anisotropic_streaming_matrices(M, θ -> surface_vF_angle(s, θ))
end

# Fallback with helpful error for unknown surface types
function streaming_matrices(M::Int, s::AbstractFermiSurface2D)
    error("streaming_matrices: unsupported surface type $(typeof(s)). " *
          "Implement surface_vF_angle(s::$(typeof(s)), θ) to support this surface type.")
end

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
    n = 1 + 2 * M
    N_quad = max(4 * n, 128)
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

