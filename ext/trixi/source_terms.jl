# source terms for 2D FermiHarmonics equations



# ======================================================================================================================
# Source Terms
# ======================================================================================================================

# The primary interface for source terms in the equations.
# This is currently identical to the physical scattering model.

"""
    source_terms(u, x, t, equations) -> SVector

Compute physical scattering source terms.
"""
@inline function source_terms(u, x, t, equations::AbstractFermiTransportEquations2D)::SVector
    return physical_sources(u, x, t, equations)
end

@inline function source_terms(u, gradients, x, t,
                              equations_parabolic::ElectrostaticGradientEquation2D)::SVector
    equations = equations_parabolic.equations_hyperbolic
    n = length(u)
    out = MVector{n, Float64}(undef)
    electrostatic_force_sources!(out, u, gradients, equations)
    return SVector(out)
end

# ======================================================================================================================
# Physical Scattering Sources
# ======================================================================================================================

# The primary source term in the FermiHarmonics equation comes from physical scattering
# processes. Scattering is differentiated by harmonic mode: monopole (charge conservation)
# is undamped, dipole experiences momentum-relaxing scattering, and higher harmonics
# experience full scattering (both momentum-relaxing and momentum-conserving).

"""
    physical_sources(u, x, t, equations) -> SVector

Source terms from a BGK-type approximation to the collision integral.  We do not damp 
``a_0`` (density). The momentum modes ``a_1, b_1`` are damped at the momentum-relaxing rate 
``\\gamma_mr``, while higher harmonics are damped at the full scattering rate 
``\\gamma_mr + \\gamma_mc``.
"""
@inline function physical_sources(u, x, t, equations::AbstractFermiTransportEquations2D)::SVector
    n = length(u)
    out = MVector{n, Float64}(undef)
    collision_sources!(out, u, equations)
    return SVector(out)
end

@inline function physical_sources(u, x, t, equations::MultiBandFermiHarmonics2D)::SVector
    n = length(u)
    out = MVector{n, Float64}(undef)
    fill!(out, 0.0)

    local_nvars = band_nvars(equations)
    M = equations.max_harmonic
    @inbounds for (band_index, band) in enumerate(equations.bands)
        offset = band_offset(equations, band_index)
        out[offset + 1] = 0.0

        if local_nvars >= 3
            out[offset + cosine_index(1)] = -band.gamma_mr * Float64(u[offset + cosine_index(1)])
            out[offset + sine_index(1)] = -band.gamma_mr * Float64(u[offset + sine_index(1)])
        end

        if M >= 2
            gamma_hi = band.gamma_mr + band.gamma_mc
            for m in 2:M
                out[offset + cosine_index(m)] = -gamma_hi * Float64(u[offset + cosine_index(m)])
                out[offset + sine_index(m)] = -gamma_hi * Float64(u[offset + sine_index(m)])
            end
        end
    end

    if equations.gamma_drag > 0.0 && band_count(equations) >= 2 && local_nvars >= 3
        nb = band_count(equations)
        if nb == 2
            # Original 2-band formula preserved exactly (normalization: w1²+w2²)
            band1 = equations.bands[1]
            band2 = equations.bands[2]
            w1 = band_momentum_weight(band1)
            w2 = band_momentum_weight(band2)
            norm_sq = w1^2 + w2^2
            if norm_sq > 0.0
                drag_scale = equations.gamma_drag / norm_sq

                a1_1 = Float64(u[band_global_cosine_index(equations, 1, 1)])
                a1_2 = Float64(u[band_global_cosine_index(equations, 2, 1)])
                relative_a1 = w2 * a1_1 - w1 * a1_2
                drag_a1_1 = -drag_scale * w2 * relative_a1
                drag_a1_2 = drag_scale * w1 * relative_a1
                out[band_global_cosine_index(equations, 1, 1)] += drag_a1_1
                out[band_global_cosine_index(equations, 2, 1)] += drag_a1_2

                b1_1 = Float64(u[band_global_sine_index(equations, 1, 1)])
                b1_2 = Float64(u[band_global_sine_index(equations, 2, 1)])
                relative_b1 = w2 * b1_1 - w1 * b1_2
                drag_b1_1 = -drag_scale * w2 * relative_b1
                drag_b1_2 = drag_scale * w1 * relative_b1
                out[band_global_sine_index(equations, 1, 1)] += drag_b1_1
                out[band_global_sine_index(equations, 2, 1)] += drag_b1_2
            end
        else
            # N≥3 bands: mean-field drag — each band's dipole relaxes toward the
            # momentum-weighted mean velocity. Conserves total weighted momentum.
            # Note: normalization uses w_total (differs from 2-band formula above).
            ws = ntuple(i -> band_momentum_weight(equations.bands[i]), nb)
            w_total = sum(ws)
            if w_total > 0.0
                p_mean_x = sum(ws[i] * Float64(u[band_global_cosine_index(equations, i, 1)])
                               for i in 1:nb) / w_total
                p_mean_y = sum(ws[i] * Float64(u[band_global_sine_index(equations, i, 1)])
                               for i in 1:nb) / w_total
                for i in 1:nb
                    out[band_global_cosine_index(equations, i, 1)] +=
                        -equations.gamma_drag * (Float64(u[band_global_cosine_index(equations, i, 1)]) - p_mean_x)
                    out[band_global_sine_index(equations, i, 1)] +=
                        -equations.gamma_drag * (Float64(u[band_global_sine_index(equations, i, 1)]) - p_mean_y)
                end
            end
        end
    end

    return SVector(out)
end

@inline function nonlinear_bgk_sources(u, equations::AbstractFermiTransportEquations2D)::SVector
    return nonlinear_bgk_sources(u, Val(:dispatch), equations)
end

@inline function nonlinear_bgk_sources(u, ::Val{:dispatch}, equations::FermiAngles2D)::SVector
    n = length(u)
    out = MVector{n, Float64}(undef)
    drift_equilibrium = MVector{n, Float64}(undef)
    isotropic_equilibrium = MVector{n, Float64}(undef)

    mu, velocity = nonlinear_collision_is_two_rate_bgk(equations) ?
        recover_mu_u_two_rate(u, equations) :
        recover_mu_u(u, equations)
    local_equilibrium_state!(drift_equilibrium, mu, velocity, equations)
    isotropic_equilibrium_state!(isotropic_equilibrium, mu, equations)

    gamma_mr = equations.gamma_mr
    gamma_mc = equations.gamma_mc
    @inbounds for i in 1:n
        out[i] = -gamma_mr * (u[i] - isotropic_equilibrium[i]) -
                 gamma_mc * (u[i] - drift_equilibrium[i])
    end
    return SVector(out)
end

@inline nonlinear_mode_rate(m::Int, equations::FermiHarmonics2D) =
    mode_rate(mode_profile(equations.model.collision), m)

@inline function nonlinear_bgk_sources(u, ::Val{:dispatch}, equations::FermiHarmonics2D)::SVector
    t0 = nonlinear_timing_enabled() ? time_ns() : UInt64(0)
    n = length(u)
    out = MVector{n, Float64}(undef)
    mu, velocity = recover_mu_u(u, equations)

    gamma_mr = equations.gamma_mr
    gamma_mode_2 = n > 3 ? nonlinear_mode_rate(2, equations) : 0.0
    ux, uy = velocity
    eq2c = 0.5 * equations.mass * (ux * ux - uy * uy)
    eq2s = equations.mass * ux * uy
    max_harmonic = (n - 1) ÷ 2
    @inbounds begin
        out[1] = 0.0
        if n >= 3
            out[2] = -gamma_mr * Float64(u[2])
            out[3] = -gamma_mr * Float64(u[3])
        end
        if max_harmonic >= 2
            ci = cosine_index(2)
            si = sine_index(2)
            out[ci] = -(gamma_mr + gamma_mode_2) * Float64(u[ci]) + gamma_mode_2 * eq2c
            out[si] = -(gamma_mr + gamma_mode_2) * Float64(u[si]) + gamma_mode_2 * eq2s
        end
        for m in 3:max_harmonic
            gamma_mode = nonlinear_mode_rate(m, equations)
            ci = cosine_index(m)
            si = sine_index(m)
            out[ci] = -(gamma_mr + gamma_mode) * Float64(u[ci])
            out[si] = -(gamma_mr + gamma_mode) * Float64(u[si])
        end
    end
    if nonlinear_timing_enabled()
        record_nonlinear_timing!(:bgk, time_ns() - t0)
    end
    return SVector(out)
end

function collision_sources!(out::AbstractVector{Float64}, u, equations::NonlinearFermiHarmonics2D)
    copyto!(out, nonlinear_bgk_sources(u, equations))
    return out
end

function collision_sources!(out::AbstractVector{Float64}, u, equations::LinearFermiHarmonics2D)
    n = length(u)
    profile = mode_profile(equations.model.collision)
    @inbounds begin
        out[1] = 0.0
        if n >= 3
            gamma_mr = equations.gamma_mr
            out[2] = -gamma_mr * Float64(u[2])
            out[3] = -gamma_mr * Float64(u[3])
        end
        if n > 3
            max_harmonic = (n - 1) ÷ 2
            for m in 2:max_harmonic
                gamma_mode = equations.gamma_mr + mode_rate(profile, m)
                ci = cosine_index(m)
                si = sine_index(m)
                out[ci] = -gamma_mode * Float64(u[ci])
                out[si] = -gamma_mode * Float64(u[si])
            end
        end
    end
    return out
end

function collision_sources!(out::AbstractVector{Float64}, u, equations::MultiBandFermiHarmonics2D)
    copyto!(out, physical_sources(u, nothing, 0.0, equations))
    return out
end

function collision_sources!(out::AbstractVector{Float64}, u, equations::FermiAngles2D)
    copyto!(out, nonlinear_bgk_sources(u, equations))
    return out
end
