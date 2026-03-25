# Thread-local cache structs and factory functions for angle-grid and nonlinear
# harmonic transport. Requires FFTW, so this file lives in the Trixi extension.

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
