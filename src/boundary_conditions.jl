# custom boundary conditions for 2D FermiHarmonics equations

@inline function get_bc_thread_buffer!(buffers::Vector{Vector{Float64}}, nvars::Int)
    tid = Threads.threadid()
    @inbounds buf = buffers[tid]
    if length(buf) != nvars
        buf = Vector{Float64}(undef, nvars)
        buffers[tid] = buf
    end
    return buf
end

@inline function ensure_bc_state_vector(u_inner, cache::BCProjectorCache)
    u_inner isa AbstractVector{Float64} && return u_inner
    nvars = length(u_inner)
    buf = get_bc_thread_buffer!(cache.state_buffers, nvars)
    @inbounds for i in 1:nvars
        buf[i] = u_inner[i]
    end
    return buf
end

@inline nonlinear_bias_scale(equations::AbstractFermiTransportEquations2D) =
    equations.mu0 * (1.0 + abs(equations.electrostatic_coupling))

@inline nonlinear_electrochemical_bias(bias::Real, equations::AbstractFermiTransportEquations2D) =
    Float64(bias) * nonlinear_bias_scale(equations)


# ======================================================================================================================
# Maxwell and Ohmic BC implementations
# ======================================================================================================================

# these functions implement the types of boundary condition defined above

"""
    maxwell_wall!(out, state, unit_normal, P_in, p, target) -> out

Apply Maxwell blending boundary condition in-place: blend diffuse and specular targets, then project onto incoming characteristics.
"""
function maxwell_wall!(
    out::AbstractVector{Float64},
    state::AbstractVector{Float64},
    unit_normal::SVector{2, Float64},
    P_in::AbstractSparseMatrix,
    p::Real,
    target::AbstractVector{Float64},
    equations::Union{FermiHarmonics2D, MultiBandFermiHarmonics2D},
)
    diffuse_target!(target, state, unit_normal, equations)
    specular_target!(out, state, unit_normal, equations)
    p_scatter = Float64(p)
    one_minus = 1.0 - p_scatter
    N = length(state)
    @inbounds for i in 1:N
        target[i] = p_scatter * target[i] + one_minus * out[i]
    end
    apply_projector!(out, state, target, P_in)
    return out
end

function maxwell_wall!(
    out::AbstractVector{Float64},
    state::AbstractVector{Float64},
    unit_normal::SVector{2, Float64},
    P_in::AbstractSparseMatrix,
    p::Real,
    target::AbstractVector{Float64},
)
    diffuse_target!(target, state, unit_normal)
    specular_target!(out, state, unit_normal)
    p_scatter = Float64(p)
    one_minus = 1.0 - p_scatter
    N = length(state)
    @inbounds for i in 1:N
        target[i] = p_scatter * target[i] + one_minus * out[i]
    end
    apply_projector!(out, state, target, P_in)
    return out
end


"""
    ohmic_contact!(out, state, unit_normal, P_in, p_ohmic_absorb, bias, target) -> out 
Apply ohmic contact boundary condition in-place: blend diffuse (with clamped monopole) and specular targets, then project onto incoming characteristics.
"""
function ohmic_contact!(out::AbstractVector{Float64},
                        state::AbstractVector{Float64},
                        unit_normal::SVector{2, Float64},
                        P_in::AbstractSparseMatrix,
                        p_ohmic_absorb::Real,
                        bias::Real,
                        target::AbstractVector{Float64},
                        equations::Union{FermiHarmonics2D, MultiBandFermiHarmonics2D})
    diffuse_target!(target, state, unit_normal, equations)
    contact_bias_target!(target, Float64(bias), equations)
    specular_target!(out, state, unit_normal, equations)
    p_absorb = Float64(p_ohmic_absorb)
    one_minus = 1.0 - p_absorb
    N = length(state)
    @inbounds for i in 1:N
        target[i] = p_absorb * target[i] + one_minus * out[i]
    end
    apply_projector!(out, state, target, P_in)
    return out
end

function ohmic_contact!(out::AbstractVector{Float64},
                        state::AbstractVector{Float64},
                        unit_normal::SVector{2, Float64},
                        P_in::AbstractSparseMatrix,
                        p_ohmic_absorb::Real,
                        bias::Real,
                        target::AbstractVector{Float64})
    diffuse_target!(target, state, unit_normal)
    @inbounds target[1] = Float64(bias)
    specular_target!(out, state, unit_normal)
    p_absorb = Float64(p_ohmic_absorb)
    one_minus = 1.0 - p_absorb
    N = length(state)
    @inbounds for i in 1:N
        target[i] = p_absorb * target[i] + one_minus * out[i]
    end
    apply_projector!(out, state, target, P_in)
    return out
end

@inline function diffuse_target!(
    target::AbstractVector{Float64},
    state::AbstractVector{Float64},
    unit_normal::SVector{2, Float64},
    ::FermiHarmonics2D,
)::AbstractVector{Float64}
    return diffuse_target!(target, state, unit_normal)
end

@inline function diffuse_target!(
    target::AbstractVector{Float64},
    state::AbstractVector{Float64},
    unit_normal::SVector{2, Float64},
    equations::MultiBandFermiHarmonics2D,
)::AbstractVector{Float64}
    local_nvars = band_nvars(equations)
    @inbounds for band_index in 1:band_count(equations)
        offset = band_offset(equations, band_index)
        diffuse_target!(
            @view(target[(offset + 1):(offset + local_nvars)]),
            @view(state[(offset + 1):(offset + local_nvars)]),
            unit_normal,
        )
    end
    return target
end

@inline function specular_target!(
    target::AbstractVector{Float64},
    state::AbstractVector{Float64},
    unit_normal::SVector{2, Float64},
    ::FermiHarmonics2D,
)::AbstractVector{Float64}
    return specular_target!(target, state, unit_normal)
end

@inline function specular_target!(
    target::AbstractVector{Float64},
    state::AbstractVector{Float64},
    unit_normal::SVector{2, Float64},
    equations::MultiBandFermiHarmonics2D,
)::AbstractVector{Float64}
    local_nvars = band_nvars(equations)
    @inbounds for band_index in 1:band_count(equations)
        offset = band_offset(equations, band_index)
        specular_target!(
            @view(target[(offset + 1):(offset + local_nvars)]),
            @view(state[(offset + 1):(offset + local_nvars)]),
            unit_normal,
        )
    end
    return target
end

@inline function contact_bias_target!(
    target::AbstractVector{Float64},
    bias::Float64,
    ::FermiHarmonics2D,
)
    target[1] = bias
    return target
end

@inline function contact_bias_target!(
    target::AbstractVector{Float64},
    bias::Float64,
    equations::MultiBandFermiHarmonics2D,
)
    @inbounds for band_index in 1:band_count(equations)
        target[band_offset(equations, band_index) + 1] = bias
    end
    return target
end

function nonlinear_diffuse_incoming_value(
    state_samples::Vector{ComplexF64},
    unit_normal::SVector{2, Float64},
    equations::FermiHarmonics2D,
    tol::Float64,
)
    data = nonlinear_data(equations)
    nx, ny = unit_normal
    dtheta = 2.0 * pi / data.theta_count
    outgoing_flux = 0.0
    incoming_weight = 0.0

    @inbounds for j in eachindex(state_samples)
        projection = nx * data.cos_theta[j] + ny * data.sin_theta[j]
        if projection > tol
            outgoing_flux += projection * parabolic_shifted_flux(real(state_samples[j]), equations)
        elseif projection < -tol
            incoming_weight -= projection
        end
    end

    incoming_weight *= dtheta
    incoming_weight > 0.0 || return 0.0
    outgoing_flux *= dtheta
    return quadratic_shifted_flux_inverse(outgoing_flux / incoming_weight, equations)
end

function nonlinear_diffuse_incoming_value(
    state_samples::Vector{ComplexF64},
    face_data::NonlinearBoundaryFaceData,
    equations::FermiHarmonics2D,
    tol::Float64,
)
    t0 = nonlinear_timing_enabled() ? time_ns() : UInt64(0)
    outgoing_flux = 0.0
    @inbounds for j in eachindex(state_samples)
        projection = face_data.projections[j]
        if projection > tol
            outgoing_flux += projection * quadratic_shifted_flux(real(state_samples[j]), equations)
        end
    end

    face_data.incoming_weight > 0.0 || return 0.0
    outgoing_flux *= 2.0 * pi / nonlinear_data(equations).theta_count
    incoming_value = quadratic_shifted_flux_inverse(outgoing_flux / face_data.incoming_weight, equations)
    if nonlinear_timing_enabled()
        record_nonlinear_timing!(:diffuse, time_ns() - t0)
    end
    return incoming_value
end

function nonlinear_boundary_samples!(
    out::AbstractVector{Float64},
    state::AbstractVector{Float64},
    state_samples::Vector{ComplexF64},
    face_data::NonlinearBoundaryFaceData,
    incoming_value::Float64,
    specular_weight::Float64,
    equations::FermiHarmonics2D,
)
    cache = get_nonlinear_cache(equations)
    copy!(out, state)
    diffuse_weight = 1.0 - specular_weight
    @inbounds for j in eachindex(state_samples)
        if face_data.incoming_mask[j]
            specular_value = specular_weight > 0.0 ? real(apply_specular_stencil(state_samples, face_data, j)) : 0.0
            incoming_sample = diffuse_weight * incoming_value + specular_weight * specular_value
            cache.work_samples[j] = ComplexF64(incoming_sample - real(state_samples[j]), 0.0)
        else
            cache.work_samples[j] = 0.0 + 0.0im
        end
    end

    apply_sample_to_harmonics_transform!(cache.real_scratch, cache.work_samples, face_data)
    @inbounds for i in eachindex(out)
        out[i] += cache.real_scratch[i]
    end
    return out
end

function nonlinear_boundary_flux!(
    out_flux::AbstractVector{Float64},
    state_samples::Vector{ComplexF64},
    face_data::NonlinearBoundaryFaceData,
    normal::SVector{2, Float64},
    incoming_value::Float64,
    specular_weight::Float64,
    equations::FermiHarmonics2D,
)
    cache = get_nonlinear_cache(equations)
    normal_x, normal_y = normal
    scale = hypot(normal_x, normal_y)
    diffuse_weight = 1.0 - specular_weight
    @inbounds for j in eachindex(state_samples)
        phi_trace = real(state_samples[j])
        if face_data.incoming_mask[j]
            specular_value = specular_weight > 0.0 ? real(apply_specular_stencil(state_samples, face_data, j)) : 0.0
            phi_trace = diffuse_weight * incoming_value + specular_weight * specular_value
        end
        directional = scale * face_data.projections[j]
        cache.scratch_samples[j] = ComplexF64(
            directional * quadratic_shifted_flux(phi_trace, equations),
            0.0,
        )
    end

    return apply_sample_to_harmonics_transform!(out_flux, cache.scratch_samples, face_data)
end

function nonlinear_maxwell_wall!(
    out::AbstractVector{Float64},
    state::AbstractVector{Float64},
    unit_normal::SVector{2, Float64},
    p_scatter::Real,
    target::AbstractVector{Float64},
    equations::FermiHarmonics2D,
    tol::Float64,
    face_data::Union{Nothing, NonlinearBoundaryFaceData} = nothing,
)
    t0 = nonlinear_timing_enabled() ? time_ns() : UInt64(0)
    cache = get_nonlinear_cache(equations)
    harmonic_state_to_samples!(cache.samples, state, equations)
    local_face_data = isnothing(face_data) ? build_nonlinear_face_data(equations, unit_normal, tol) : face_data
    diffuse_value = nonlinear_diffuse_incoming_value(cache.samples, local_face_data, equations, tol)
    result = nonlinear_boundary_samples!(
        out,
        state,
        cache.samples,
        local_face_data,
        diffuse_value,
        1.0 - Float64(p_scatter),
        equations,
    )
    if nonlinear_timing_enabled()
        record_nonlinear_timing!(:boundary, time_ns() - t0)
    end
    return result
end

function nonlinear_maxwell_wall_flux!(
    out_flux::AbstractVector{Float64},
    state::AbstractVector{Float64},
    normal::SVector{2, Float64},
    unit_normal::SVector{2, Float64},
    p_scatter::Real,
    target::AbstractVector{Float64},
    equations::FermiHarmonics2D,
    tol::Float64,
    face_data::Union{Nothing, NonlinearBoundaryFaceData} = nothing,
)
    cache = get_nonlinear_cache(equations)
    harmonic_state_to_samples!(cache.samples, state, equations)
    local_face_data = isnothing(face_data) ? build_nonlinear_face_data(equations, unit_normal, tol) : face_data
    diffuse_value = nonlinear_diffuse_incoming_value(cache.samples, local_face_data, equations, tol)
    return nonlinear_boundary_flux!(
        out_flux,
        cache.samples,
        local_face_data,
        normal,
        diffuse_value,
        1.0 - Float64(p_scatter),
        equations,
    )
end

function nonlinear_ohmic_incoming_value(
    state_samples::Vector{ComplexF64},
    face_data::NonlinearBoundaryFaceData,
    p_ohmic_absorb::Real,
    bias::Real,
    equations::FermiHarmonics2D,
)
    electrochemical_bias = nonlinear_electrochemical_bias(bias, equations)
    if !nonlinear_has_electrostatic_force(equations)
        return electrochemical_bias
    end

    specular_weight = 1.0 - Float64(p_ohmic_absorb)
    diffuse_weight = 1.0 - specular_weight
    base_sum = 0.0
    phi0_coeff = 0.0
    inv_ntheta = 1.0 / nonlinear_data(equations).theta_count

    @inbounds for j in eachindex(state_samples)
        if face_data.incoming_mask[j]
            specular_value = specular_weight > 0.0 ? real(apply_specular_stencil(state_samples, face_data, j)) : 0.0
            base_sum += specular_weight * specular_value
            phi0_coeff += diffuse_weight * inv_ntheta
        else
            base_sum += real(state_samples[j])
        end
    end

    phi0_base = base_sum * inv_ntheta
    denominator = 1.0 + equations.electrostatic_coupling * phi0_coeff
    abs(denominator) > 1.0e-14 || throw(DomainError(
        denominator,
        "electrochemical contact solve became singular",
    ))
    return (electrochemical_bias - equations.electrostatic_coupling * phi0_base) / denominator
end

function nonlinear_ohmic_incoming_value(
    state::AbstractVector{Float64},
    unit_normal::SVector{2, Float64},
    p_ohmic_absorb::Real,
    bias::Real,
    target::AbstractVector{Float64},
    equations::FermiHarmonics2D,
    tol::Float64,
    face_data::Union{Nothing, NonlinearBoundaryFaceData} = nothing,
)
    cache = get_nonlinear_cache(equations)
    harmonic_state_to_samples!(cache.samples, state, equations)
    local_face_data = isnothing(face_data) ? build_nonlinear_face_data(equations, unit_normal, tol) : face_data
    return nonlinear_ohmic_incoming_value(
        cache.samples,
        local_face_data,
        p_ohmic_absorb,
        bias,
        equations,
    )
end

function nonlinear_ohmic_contact!(
    out::AbstractVector{Float64},
    state::AbstractVector{Float64},
    unit_normal::SVector{2, Float64},
    p_ohmic_absorb::Real,
    bias::Real,
    target::AbstractVector{Float64},
    equations::FermiHarmonics2D,
    tol::Float64,
    face_data::Union{Nothing, NonlinearBoundaryFaceData} = nothing,
)
    t0 = nonlinear_timing_enabled() ? time_ns() : UInt64(0)
    cache = get_nonlinear_cache(equations)
    harmonic_state_to_samples!(cache.samples, state, equations)
    local_face_data = isnothing(face_data) ? build_nonlinear_face_data(equations, unit_normal, tol) : face_data
    incoming_value = nonlinear_ohmic_incoming_value(
        cache.samples,
        local_face_data,
        p_ohmic_absorb,
        bias,
        equations,
    )
    result = nonlinear_boundary_samples!(
        out,
        state,
        cache.samples,
        local_face_data,
        incoming_value,
        1.0 - Float64(p_ohmic_absorb),
        equations,
    )
    if nonlinear_timing_enabled()
        record_nonlinear_timing!(:boundary, time_ns() - t0)
    end
    return result
end

function nonlinear_ohmic_contact_flux!(
    out_flux::AbstractVector{Float64},
    state::AbstractVector{Float64},
    normal::SVector{2, Float64},
    unit_normal::SVector{2, Float64},
    p_ohmic_absorb::Real,
    bias::Real,
    target::AbstractVector{Float64},
    equations::FermiHarmonics2D,
    tol::Float64,
    face_data::Union{Nothing, NonlinearBoundaryFaceData} = nothing,
)
    cache = get_nonlinear_cache(equations)
    harmonic_state_to_samples!(cache.samples, state, equations)
    local_face_data = isnothing(face_data) ? build_nonlinear_face_data(equations, unit_normal, tol) : face_data
    incoming_value = nonlinear_ohmic_incoming_value(
        cache.samples,
        local_face_data,
        p_ohmic_absorb,
        bias,
        equations,
    )
    return nonlinear_boundary_flux!(
        out_flux,
        cache.samples,
        local_face_data,
        normal,
        incoming_value,
        1.0 - Float64(p_ohmic_absorb),
        equations,
    )
end

@inline function wrap_periodic_index(index::Integer, count::Integer)
    return mod1(Int(index), Int(count))
end

function cubic_periodic_stencil(theta_ref::Float64, theta_count::Int, dtheta::Float64)
    position = theta_ref / dtheta
    base_index = floor(Int, position) + 1
    t = position - (base_index - 1)
    indices = (
        wrap_periodic_index(base_index - 1, theta_count),
        wrap_periodic_index(base_index, theta_count),
        wrap_periodic_index(base_index + 1, theta_count),
        wrap_periodic_index(base_index + 2, theta_count),
    )
    weights = (
        -t * (t - 1.0) * (t - 2.0) / 6.0,
        (t + 1.0) * (t - 1.0) * (t - 2.0) / 2.0,
        -(t + 1.0) * t * (t - 2.0) / 2.0,
        (t + 1.0) * t * (t - 1.0) / 6.0,
    )
    return indices, weights
end

function build_nonlinear_face_data(
    equations::FermiHarmonics2D,
    unit_normal::SVector{2, Float64},
    tol::Float64,
)
    data = nonlinear_data(equations)
    theta_count = data.theta_count
    incoming_mask = falses(theta_count)
    stencil_indices = Matrix{Int}(undef, 4, theta_count)
    stencil_weights = Matrix{Float64}(undef, 4, theta_count)
    projections = Vector{Float64}(undef, theta_count)
    alpha = atan(unit_normal[2], unit_normal[1])
    dtheta = 2.0 * pi / theta_count
    incoming_weight = 0.0
    nvars = length(get_nonlinear_cache(equations).real_work)
    max_harmonic = (nvars - 1) ÷ 2
    sample_to_harmonics = Matrix{Float64}(undef, nvars, theta_count)
    inv_theta_count = 1.0 / theta_count

    @inbounds for j in 1:theta_count
        projection = unit_normal[1] * data.cos_theta[j] + unit_normal[2] * data.sin_theta[j]
        projections[j] = projection
        incoming_mask[j] = projection < -tol
        if projection < -tol
            incoming_weight -= projection
        end
        theta_ref = mod(2.0 * alpha + pi - data.theta[j], 2.0 * pi)
        indices, weights = cubic_periodic_stencil(theta_ref, theta_count, dtheta)
        stencil_indices[:, j] .= indices
        stencil_weights[:, j] .= weights
        sample_to_harmonics[1, j] = 2.0 * inv_theta_count
        for m in 1:max_harmonic
            sample_to_harmonics[cosine_index(m), j] = 2.0 * inv_theta_count * cos(m * data.theta[j])
            sample_to_harmonics[sine_index(m), j] = 2.0 * inv_theta_count * sin(m * data.theta[j])
        end
    end

    return NonlinearBoundaryFaceData(
        unit_normal,
        incoming_mask,
        stencil_indices,
        stencil_weights,
        projections,
        incoming_weight * dtheta,
        sample_to_harmonics,
    )
end

function build_nonlinear_face_data(
    equations::FermiAngles2D,
    unit_normal::SVector{2, Float64},
    tol::Float64,
)
    data = nonlinear_data(equations)
    theta_count = data.theta_count
    incoming_mask = falses(theta_count)
    stencil_indices = Matrix{Int}(undef, 4, theta_count)
    stencil_weights = Matrix{Float64}(undef, 4, theta_count)
    projections = Vector{Float64}(undef, theta_count)
    alpha = atan(unit_normal[2], unit_normal[1])
    dtheta = data.weight
    incoming_weight = 0.0

    @inbounds for j in 1:theta_count
        projection = unit_normal[1] * data.cos_theta[j] + unit_normal[2] * data.sin_theta[j]
        projections[j] = projection
        incoming_mask[j] = projection < -tol
        if projection < -tol
            incoming_weight -= projection
        end
        theta_ref = mod(2.0 * alpha + pi - data.theta[j], 2.0 * pi)
        indices, weights = cubic_periodic_stencil(theta_ref, theta_count, dtheta)
        stencil_indices[:, j] .= indices
        stencil_weights[:, j] .= weights
    end

    return NonlinearBoundaryFaceData(
        unit_normal,
        incoming_mask,
        stencil_indices,
        stencil_weights,
        projections,
        incoming_weight * data.weight,
        nothing,
    )
end

@inline function apply_specular_stencil(
    state::AbstractVector,
    face_data::NonlinearBoundaryFaceData,
    angle_index::Int,
)
    return face_data.stencil_weights[1, angle_index] * state[face_data.stencil_indices[1, angle_index]] +
           face_data.stencil_weights[2, angle_index] * state[face_data.stencil_indices[2, angle_index]] +
           face_data.stencil_weights[3, angle_index] * state[face_data.stencil_indices[3, angle_index]] +
           face_data.stencil_weights[4, angle_index] * state[face_data.stencil_indices[4, angle_index]]
end

@inline function apply_sample_to_harmonics_transform!(
    out::AbstractVector{Float64},
    sample_values::AbstractVector,
    face_data::NonlinearBoundaryFaceData,
)
    transform = face_data.sample_to_harmonics
    isnothing(transform) && error("harmonic sample transform missing from face cache")
    nvars = size(transform, 1)
    theta_count = size(transform, 2)
    @inbounds for i in 1:nvars
        acc = 0.0
        for j in 1:theta_count
            acc += transform[i, j] * real(sample_values[j])
        end
        out[i] = acc
    end
    return out
end

function nonlinear_diffuse_incoming_value(
    state::AbstractVector{<:Real},
    face_data::NonlinearBoundaryFaceData,
    equations::FermiAngles2D,
    tol::Float64,
)
    data = nonlinear_data(equations)
    outgoing_flux = 0.0
    incoming_weight = 0.0
    nx, ny = face_data.unit_normal
    @inbounds for j in eachindex(state)
        projection = face_data.projections[j]
        if projection > tol
            outgoing_flux += projection * parabolic_shifted_flux(state[j], equations)
        end
    end
    incoming_weight = face_data.incoming_weight
    incoming_weight > 0.0 || return 0.0
    outgoing_flux *= data.weight
    return parabolic_shifted_flux_inverse(outgoing_flux / incoming_weight, equations)
end

function nonlinear_boundary_samples!(
    out::AbstractVector{Float64},
    state::AbstractVector{Float64},
    face_data::NonlinearBoundaryFaceData,
    incoming_value::Float64,
    specular_weight::Float64,
    equations::FermiAngles2D,
)
    copy!(out, state)
    diffuse_weight = 1.0 - specular_weight
    @inbounds for j in eachindex(out)
        if face_data.incoming_mask[j]
            specular_value = apply_specular_stencil(state, face_data, j)
            out[j] = diffuse_weight * incoming_value + specular_weight * specular_value
        end
    end
    return out
end

function nonlinear_maxwell_wall!(
    out::AbstractVector{Float64},
    state::AbstractVector{Float64},
    unit_normal::SVector{2, Float64},
    p_scatter::Real,
    target::AbstractVector{Float64},
    equations::FermiAngles2D,
    tol::Float64,
    face_data::Union{Nothing, NonlinearBoundaryFaceData} = nothing,
)
    local_face_data = isnothing(face_data) ? build_nonlinear_face_data(equations, unit_normal, tol) : face_data
    diffuse_value = nonlinear_diffuse_incoming_value(state, local_face_data, equations, tol)
    return nonlinear_boundary_samples!(
        out,
        state,
        local_face_data,
        diffuse_value,
        1.0 - Float64(p_scatter),
        equations,
    )
end

function nonlinear_ohmic_incoming_value(
    state::AbstractVector{Float64},
    face_data::NonlinearBoundaryFaceData,
    p_ohmic_absorb::Real,
    bias::Real,
    equations::FermiAngles2D,
)
    electrochemical_bias = nonlinear_electrochemical_bias(bias, equations)
    if !nonlinear_has_electrostatic_force(equations)
        return electrochemical_bias
    end

    specular_weight = 1.0 - Float64(p_ohmic_absorb)
    diffuse_weight = 1.0 - specular_weight
    inv_ntheta = 1.0 / nonlinear_data(equations).theta_count
    base_sum = 0.0
    phi0_coeff = 0.0

    @inbounds for j in eachindex(state)
        if face_data.incoming_mask[j]
            specular_value = specular_weight > 0.0 ? apply_specular_stencil(state, face_data, j) : 0.0
            base_sum += specular_weight * specular_value
            phi0_coeff += diffuse_weight * inv_ntheta
        else
            base_sum += state[j]
        end
    end

    phi0_base = base_sum * inv_ntheta
    denominator = 1.0 + equations.electrostatic_coupling * phi0_coeff
    abs(denominator) > 1.0e-14 || throw(DomainError(
        denominator,
        "electrochemical contact solve became singular",
    ))
    return (electrochemical_bias - equations.electrostatic_coupling * phi0_base) / denominator
end

function nonlinear_ohmic_contact!(
    out::AbstractVector{Float64},
    state::AbstractVector{Float64},
    unit_normal::SVector{2, Float64},
    p_ohmic_absorb::Real,
    bias::Real,
    target::AbstractVector{Float64},
    equations::FermiAngles2D,
    tol::Float64,
    face_data::Union{Nothing, NonlinearBoundaryFaceData} = nothing,
)
    local_face_data = isnothing(face_data) ? build_nonlinear_face_data(equations, unit_normal, tol) : face_data
    incoming_value = nonlinear_ohmic_incoming_value(
        state,
        local_face_data,
        p_ohmic_absorb,
        bias,
        equations,
    )
    return nonlinear_boundary_samples!(
        out,
        state,
        local_face_data,
        incoming_value,
        1.0 - Float64(p_ohmic_absorb),
        equations,
    )
end

@inline function (bc::MaxwellWallBC)(
    flux_inner,
    u_inner,
    normal_direction::AbstractVector,
    x,
    t,
    operator_type::Trixi.Gradient,
    equations_parabolic::ElectrostaticGradientEquation2D,
)
    equations = equations_parabolic.equations_hyperbolic
    state = ensure_bc_state_vector(u_inner, bc.cache)
    out = get_bc_thread_buffer!(bc.cache.out_buffers, length(state))
    target = get_bc_thread_buffer!(bc.cache.target_buffers, length(state))
    normal = SVector(Float64(normal_direction[1]), Float64(normal_direction[2]))
    unit_n = unit_normal(normal)
    nonlinear_maxwell_wall!(out, state, unit_n, bc.p_scatter, target, equations, max(bc.tol, 1.0e-12))
    return out
end

@inline function (bc::MaxwellWallBC)(
    flux_inner,
    u_inner,
    normal_direction::AbstractVector,
    x,
    t,
    operator_type::Trixi.Divergence,
    equations_parabolic::ElectrostaticGradientEquation2D,
)
    return flux_inner
end

@inline function (bc::MaxwellWallBC)(
    flux_inner,
    u_inner,
    normal_direction::AbstractVector,
    x,
    t,
    operator_type::Trixi.Gradient,
    equations_parabolic::MeanModeGradientEquation2D,
)
    return u_inner
end

@inline function (bc::MaxwellWallBC)(
    flux_inner,
    u_inner,
    normal_direction::AbstractVector,
    x,
    t,
    operator_type::Trixi.Divergence,
    equations_parabolic::MeanModeGradientEquation2D,
)
    return flux_inner
end

@inline function (bc::OhmicContactBC)(
    flux_inner,
    u_inner,
    normal_direction::AbstractVector,
    x,
    t,
    operator_type::Trixi.Gradient,
    equations_parabolic::ElectrostaticGradientEquation2D,
)
    equations = equations_parabolic.equations_hyperbolic
    state = ensure_bc_state_vector(u_inner, bc.cache)
    out = get_bc_thread_buffer!(bc.cache.out_buffers, length(state))
    target = get_bc_thread_buffer!(bc.cache.target_buffers, length(state))
    normal = SVector(Float64(normal_direction[1]), Float64(normal_direction[2]))
    unit_n = unit_normal(normal)
    nonlinear_ohmic_contact!(
        out,
        state,
        unit_n,
        bc.p_ohmic_absorb,
        bc.bias,
        target,
        equations,
        max(bc.tol, 1.0e-12),
    )
    return out
end

@inline function (bc::OhmicContactBC)(
    flux_inner,
    u_inner,
    normal_direction::AbstractVector,
    x,
    t,
    operator_type::Trixi.Divergence,
    equations_parabolic::ElectrostaticGradientEquation2D,
)
    return flux_inner
end

@inline function (bc::OhmicContactBC)(
    flux_inner,
    u_inner,
    normal_direction::AbstractVector,
    x,
    t,
    operator_type::Trixi.Gradient,
    equations_parabolic::MeanModeGradientEquation2D,
)
    return u_inner
end

@inline function (bc::OhmicContactBC)(
    flux_inner,
    u_inner,
    normal_direction::AbstractVector,
    x,
    t,
    operator_type::Trixi.Divergence,
    equations_parabolic::MeanModeGradientEquation2D,
)
    return flux_inner
end


# ======================================================================================================================
# Diffuse and Specular Target Computation
# ======================================================================================================================

"""
    diffuse_target!(target, state, unit_normal) -> target

Compute diffuse scattering target state in-place.
"""
@inline function diffuse_target!(
    target::AbstractVector{Float64},
    state::AbstractVector{Float64},
    unit_normal::SVector{2, Float64},
)::AbstractVector{Float64}
    c1, s1 = normal_cos_sin(unit_normal)
    N = length(state)
    Mloc = (N - 1) ÷ 2
    @inbounds C = state[1]
    if Mloc >= 1
        cm = c1
        sm = s1
        @inbounds C += 0.5 * pi * (cm * state[cosine_index(1)] + sm * state[sine_index(1)])
        @inbounds for m in 2:Mloc
            cm_next = c1 * cm - s1 * sm
            sm_next = s1 * cm + c1 * sm
            cm = cm_next
            sm = sm_next
            if iseven(m)
                k = m ÷ 2
                sign = isodd(k) ? -1.0 : 1.0
                coeff = -2.0 * sign / (4.0 * k * k - 1.0)
                C += coeff * (cm * state[cosine_index(m)] + sm * state[sine_index(m)])
            end
        end
    end
    @inbounds for i in 1:N
        target[i] = 0.0
    end
    @inbounds target[1] = C
    return target
end

"""
    specular_target!(target, state, unit_normal) -> target

Compute specular reflection target state in-place.
"""
@inline function specular_target!(
    target::AbstractVector{Float64},
    state::AbstractVector{Float64},
    unit_normal::SVector{2, Float64},
)::AbstractVector{Float64}
    c1, s1 = normal_cos_sin(unit_normal)
    N = length(state)
    Mloc = (N - 1) ÷ 2
    @inbounds target[1] = state[1]
    cm = c1
    sm = s1
    @inbounds for m in 1:Mloc
        a = state[cosine_index(m)]
        b = state[sine_index(m)]
        a_rot = cm * a + sm * b
        b_rot = -sm * a + cm * b
        sign = isodd(m) ? -1.0 : 1.0
        a_rot_ref = sign * a_rot
        b_rot_ref = -sign * b_rot
        target[cosine_index(m)] = cm * a_rot_ref - sm * b_rot_ref
        target[sine_index(m)] = sm * a_rot_ref + cm * b_rot_ref
        if m < Mloc
            cm_next = c1 * cm - s1 * sm
            sm_next = s1 * cm + c1 * sm
            cm = cm_next
            sm = sm_next
        end
    end
    return target
end


# ======================================================================================================================
# Kinetic Boundary Conditions via Characteristic Projection
# ======================================================================================================================

# Here we implement the projector computation and application functions needed for the boundary conditions.  Notice that
# a discrete angle method may be more natural for these kinetic BCs.  However, in that case, we would need a large
# number of degrees of freedom even in the hydro, diffusive regimes.  So we work instead in this harmonic basis.

"""
    incoming_projector(Ax, Ay, unit_normal; tol=0.0) -> SparseMatrixCSC{Float64}

Compute incoming projector for flux Jacobian A = Ax*nx + Ay*ny as sparse matrix.
Incoming modes are eigenvectors with eigenvalues <= -tol.
"""
function incoming_projector(
    Ax::AbstractMatrix,
    Ay::AbstractMatrix,
    unit_normal::SVector{2, Float64},
    ;tol::Float64=0.0,
    monopole_indices::AbstractVector{<:Integer}=Int[1],
)::SparseMatrixCSC{Float64, Int}
    A = unit_normal[1] .* Ax .+ unit_normal[2] .* Ay
    n = size(A, 1)
    D = ones(Float64, n)
    Dinv = ones(Float64, n)
    @inbounds for index in monopole_indices
        D[Int(index)] = sqrt(2.0)
        Dinv[Int(index)] = 1 / sqrt(2.0)
    end

    Dx = Diagonal(D)
    Dxinv = Diagonal(Dinv)
    A_s = Dxinv * A * Dx
    eig = eigen(Symmetric(A_s))
    mask = eig.values .<= -tol
    
    if !any(mask)
        return spzeros(Float64, n, n)
    end

    Q_in = eig.vectors[:, mask]
    P_s = Q_in * transpose(Q_in)
    P_dense = Dx * P_s * Dxinv
    P_sparse = sparse(P_dense)
    droptol!(P_sparse, 1e-12)
    return P_sparse
end

function incoming_projector(
    equations::FermiHarmonics2D,
    unit_normal::SVector{2, Float64};
    tol::Float64 = 0.0,
)::SparseMatrixCSC{Float64, Int}
    return incoming_projector(equations.Ax, equations.Ay, unit_normal; tol=tol, monopole_indices=[1])
end

function incoming_projector(
    equations::MultiBandFermiHarmonics2D,
    unit_normal::SVector{2, Float64};
    tol::Float64 = 0.0,
)::SparseMatrixCSC{Float64, Int}
    monopole_indices = [band_offset(equations, band_index) + 1 for band_index in 1:band_count(equations)]
    return incoming_projector(equations.Ax, equations.Ay, unit_normal; tol=tol, monopole_indices=monopole_indices)
end

"""
    apply_projector!(out, state, target, P_in) -> out

Apply incoming projector: out = state + P_in * (target - state). Sparse matrix implementation.
"""
@inline function apply_projector!(
    out::AbstractVector{Float64},
    state::AbstractVector{Float64},
    target::AbstractVector{Float64},
    P_in::AbstractSparseMatrix,
)::AbstractVector{Float64}
    N = length(state)
    @inbounds for i in 1:N
        out[i] = state[i]
    end
    colptr = P_in.colptr
    rowval = P_in.rowval
    nzval = P_in.nzval
    @inbounds for j in 1:N
        dj = target[j] - state[j]
        for idx in colptr[j]:(colptr[j + 1] - 1)
            out[rowval[idx]] += nzval[idx] * dj
        end
    end
    return out
end


# ======================================================================================================================
# Projector Cache Initialization
# ======================================================================================================================

# In this section we implement the caching of projector matrices per face for P4est meshes.
# The main function is `build_projectors` which constructs the sparse projector matrices
# for all faces assigned to a given boundary condition type. These are stored in the
# `BCProjectorCache` struct defined above.

"""
    build_projectors(equations, tol, mesh, dg, cache, boundary_indexing)

Build sparse projector matrix for each boundary face in a P4est mesh.
Returns Dict{Int, SparseMatrixCSC} mapping global boundary_index to projector.
"""
function build_projectors(
    equations::FermiHarmonics2D,
    tol::Float64,
    mesh::Trixi.P4estMesh{2},
    solver,
    cache,
    boundary_indexing::Vector{Int},
)::Dict{Int, SparseMatrixCSC{Float64, Int}}
    n_nodes = Trixi.nnodes(solver)
    contravariant_vectors = cache.elements.contravariant_vectors
    boundaries = cache.boundaries
    projectors = Dict{Int, SparseMatrixCSC{Float64, Int}}()
    # Loop over all global boundary indices assigned to this BC type
    for global_idx in boundary_indexing
        element = boundaries.neighbor_ids[global_idx]
        node_indices = boundaries.node_indices[global_idx]
        direction = Trixi.indices2direction(node_indices)
        # Use first node to compute normal
        if direction == 1 || direction == 2
            i_index = (direction == 1) ? 1 : n_nodes
            j_index = 1
        else  # direction == 3 || direction == 4
            i_index = 1
            j_index = (direction == 3) ? 1 : n_nodes
        end
        normal_direction = Trixi.get_normal_direction(
            direction, contravariant_vectors, i_index, j_index, element
        )
        unit_n = unit_normal(
            SVector(Float64(normal_direction[1]), Float64(normal_direction[2]))
        )
        projectors[global_idx] = incoming_projector(equations, unit_n; tol = tol)
    end
    return projectors
end

function build_projectors(
    equations::MultiBandFermiHarmonics2D,
    tol::Float64,
    mesh::Trixi.P4estMesh{2},
    solver,
    cache,
    boundary_indexing::Vector{Int},
)::Dict{Int, SparseMatrixCSC{Float64, Int}}
    n_nodes = Trixi.nnodes(solver)
    contravariant_vectors = cache.elements.contravariant_vectors
    boundaries = cache.boundaries
    projectors = Dict{Int, SparseMatrixCSC{Float64, Int}}()
    for global_idx in boundary_indexing
        element = boundaries.neighbor_ids[global_idx]
        node_indices = boundaries.node_indices[global_idx]
        direction = Trixi.indices2direction(node_indices)
        if direction == 1 || direction == 2
            i_index = (direction == 1) ? 1 : n_nodes
            j_index = 1
        else
            i_index = 1
            j_index = (direction == 3) ? 1 : n_nodes
        end
        normal_direction = Trixi.get_normal_direction(
            direction, contravariant_vectors, i_index, j_index, element
        )
        unit_n = unit_normal(
            SVector(Float64(normal_direction[1]), Float64(normal_direction[2]))
        )
        projectors[global_idx] = incoming_projector(equations, unit_n; tol = tol)
    end
    return projectors
end

function build_nonlinear_faces(
    equations::FermiHarmonics2D,
    tol::Float64,
    mesh::Trixi.P4estMesh{2},
    solver,
    cache,
    boundary_indexing::Vector{Int},
)::Dict{Int, Any}
    n_nodes = Trixi.nnodes(solver)
    contravariant_vectors = cache.elements.contravariant_vectors
    boundaries = cache.boundaries
    nonlinear_faces = Dict{Int, Any}()
    for global_idx in boundary_indexing
        element = boundaries.neighbor_ids[global_idx]
        node_indices = boundaries.node_indices[global_idx]
        direction = Trixi.indices2direction(node_indices)
        if direction == 1 || direction == 2
            i_index = direction == 1 ? 1 : n_nodes
            j_index = 1
        else
            i_index = 1
            j_index = direction == 3 ? 1 : n_nodes
        end
        normal_direction = Trixi.get_normal_direction(
            direction, contravariant_vectors, i_index, j_index, element
        )
        unit_n = unit_normal(
            SVector(Float64(normal_direction[1]), Float64(normal_direction[2]))
        )
        nonlinear_faces[global_idx] = build_nonlinear_face_data(equations, unit_n, tol)
    end
    return nonlinear_faces
end

function build_nonlinear_faces(
    equations::FermiAngles2D,
    tol::Float64,
    mesh::Trixi.P4estMesh{2},
    solver,
    cache,
    boundary_indexing::Vector{Int},
)::Dict{Int, Any}
    n_nodes = Trixi.nnodes(solver)
    contravariant_vectors = cache.elements.contravariant_vectors
    boundaries = cache.boundaries
    nonlinear_faces = Dict{Int, Any}()
    for global_idx in boundary_indexing
        element = boundaries.neighbor_ids[global_idx]
        node_indices = boundaries.node_indices[global_idx]
        direction = Trixi.indices2direction(node_indices)
        if direction == 1 || direction == 2
            i_index = direction == 1 ? 1 : n_nodes
            j_index = 1
        else
            i_index = 1
            j_index = direction == 3 ? 1 : n_nodes
        end
        normal_direction = Trixi.get_normal_direction(
            direction, contravariant_vectors, i_index, j_index, element
        )
        unit_n = unit_normal(
            SVector(Float64(normal_direction[1]), Float64(normal_direction[2]))
        )
        nonlinear_faces[global_idx] = build_nonlinear_face_data(equations, unit_n, tol)
    end
    return nonlinear_faces
end

"""
    init_projector_cache!(semi)

Initialize the BC projector cache for all boundary conditions in the semidiscretization.

The projector cache depends only on mesh geometry, solver degree (polydeg), and boundary
locations - NOT on physics parameters. This function checks if the cache is already 
initialized and reuses it across multiple solves, which is critical for parameter sweeps.

# Performance note
Building projectors is expensive (eigendecomposition of flux Jacobians). Reusing the cache
can save significant time when running multiple cases with the same mesh/solver but 
different physics parameters (gamma_mr, gamma_mc, etc.).
"""
function init_projector_cache!(
    semi::Trixi.SemidiscretizationHyperbolic{<:Any, <:AbstractFermiTransportEquations2D},
)::Trixi.SemidiscretizationHyperbolic{<:Any, <:AbstractFermiTransportEquations2D}
    boundary_conditions = semi.boundary_conditions
    nvars = Trixi.nvariables(semi.equations)
    desired_signature = transport_is_nonlinear(semi.equations) ?
        (semi.equations.transport, nvars, nonlinear_data(semi.equations).theta_count) :
        (semi.equations.transport, nvars, 0)
    
    # Reset caches if number of variables changed (e.g., adaptive harmonics in sweeps)
    for bc in boundary_conditions.boundary_condition_types
        if bc.cache.signature != desired_signature
            empty!(bc.cache.projectors)
            empty!(bc.cache.nonlinear_faces)
            bc.cache.initialized = false
        end
        bc.cache.nvars = nvars
        bc.cache.signature = desired_signature
    end

    # Check if cache already initialized - reuse if possible
    all_initialized = all(bc.cache.initialized for bc in boundary_conditions.boundary_condition_types)
    if all_initialized
        @debug "BC projector cache already initialized, reusing" mesh = typeof(semi.mesh)
        return semi
    end
    
    # Build projectors for any uninitialized boundary conditions
    for (bc, boundary_indexing) in zip(
        boundary_conditions.boundary_condition_types,
        boundary_conditions.boundary_indices,
    )
        if !bc.cache.initialized
            if transport_is_nonlinear(semi.equations) || semi.equations isa FermiAngles2D
                new_faces = build_nonlinear_faces(
                    semi.equations, bc.tol, semi.mesh, semi.solver, semi.cache, boundary_indexing
                )
                merge!(bc.cache.nonlinear_faces, new_faces)
            else
                new_projectors = build_projectors(
                    semi.equations, bc.tol, semi.mesh, semi.solver, semi.cache, boundary_indexing
                )
                merge!(bc.cache.projectors, new_projectors)
            end
            bc.cache.initialized = true
        end
    end
    @debug "Initialized BC projector cache" mesh = typeof(semi.mesh)
    return semi
end

function init_projector_cache!(
    semi::Trixi.SemidiscretizationHyperbolicParabolic{<:Any, <:AbstractFermiTransportEquations2D},
)::Trixi.SemidiscretizationHyperbolicParabolic{<:Any, <:AbstractFermiTransportEquations2D}
    boundary_conditions = semi.boundary_conditions
    nvars = Trixi.nvariables(semi.equations)
    desired_signature = transport_is_nonlinear(semi.equations) ?
        (semi.equations.transport, nvars, nonlinear_data(semi.equations).theta_count) :
        (semi.equations.transport, nvars, 0)

    for bc in boundary_conditions.boundary_condition_types
        if bc.cache.signature != desired_signature
            empty!(bc.cache.projectors)
            empty!(bc.cache.nonlinear_faces)
            bc.cache.initialized = false
        end
        bc.cache.nvars = nvars
        bc.cache.signature = desired_signature
    end

    all_initialized = all(bc.cache.initialized for bc in boundary_conditions.boundary_condition_types)
    if all_initialized
        @debug "BC projector cache already initialized, reusing" mesh = typeof(semi.mesh)
        return semi
    end

    for (bc, boundary_indexing) in zip(
        boundary_conditions.boundary_condition_types,
        boundary_conditions.boundary_indices,
    )
        if !bc.cache.initialized
            if transport_is_nonlinear(semi.equations) || semi.equations isa FermiAngles2D
                new_faces = build_nonlinear_faces(
                    semi.equations, bc.tol, semi.mesh, semi.solver, semi.cache, boundary_indexing
                )
                merge!(bc.cache.nonlinear_faces, new_faces)
            else
                new_projectors = build_projectors(
                    semi.equations, bc.tol, semi.mesh, semi.solver, semi.cache, boundary_indexing
                )
                merge!(bc.cache.projectors, new_projectors)
            end
            bc.cache.initialized = true
        end
    end
    @debug "Initialized BC projector cache" mesh = typeof(semi.mesh)
    return semi
end

# ======================================================================================================================
# Helper Functions for Normal Vector Computation
# ======================================================================================================================

"""
    unit_normal(normal::SVector{2, Float64}) -> SVector{2, Float64}
    
Normalize a normal vector to unit length.
"""
@inline function unit_normal(normal::SVector{2, Float64})::SVector{2, Float64}
    nx, ny = normal
    nrm = hypot(nx, ny)
    if nrm < 1e-14
        return SVector(1.0, 0.0)
    else
        return SVector(nx / nrm, ny / nrm)
    end
end

"""
    _normal_cos_sin(normal) -> (cosine, sine)

Compute cosine and sine of the angle defined by the normal vector.
"""
@inline function normal_cos_sin(normal::SVector{2, Float64})::Tuple{Float64, Float64}
    nx, ny = normal
    nrm = hypot(nx, ny)
    if nrm < 1e-14
        return 1.0, 0.0
    else
        return nx / nrm, ny / nrm
    end
end
