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
@inline function source_terms(u, x, t, equations::FermiHarmonics2D)::SVector
    return physical_sources(u, x, t, equations)
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

Source terms from a tomographic approximation to the collision integral. We do not damp
``a_0`` (density). The momentum modes ``a_1, b_1`` are damped at the momentum-relaxing rate
``\\gamma_{mr}``. Higher even harmonics are damped at ``\\gamma_{mr} + \\gamma_{ee}``, while
higher odd harmonics receive the extra capped enhancement
``\\min(\\gamma_3 m^4, \\gamma_{ee})``.
"""
@inline function physical_sources(u, x, t, equations::FermiHarmonics2D)::SVector
    n = length(u)
    out = MVector{n, Float64}(undef)
    omega_c = equations.omega_c
    @inbounds begin
        # Monopole: no damping (charge conservation)
        out[1] = 0.0
        
        # Dipole: momentum-relaxing scattering only
        if n >= 3
            gamma_mr = equations.gamma_mr
            out[2] = -gamma_mr * u[2] - omega_c * u[3]
            out[3] = -gamma_mr * u[3] + omega_c * u[2]
        end
        
        # Higher harmonics: tomographic odd/even scattering rates
        if n > 3
            gamma_even = equations.gamma_mr + equations.gamma_ee
            max_harmonic = (n - 1) ÷ 2
            for m in 2:max_harmonic
                ci = cosine_index(m)
                si = sine_index(m)
                gamma_mode = if iseven(m)
                    gamma_even
                else
                    gamma_even + min(equations.gamma_3 * m^4, equations.gamma_ee)
                end
                out[ci] = -gamma_mode * u[ci] - m * omega_c * u[si]
                out[si] = -gamma_mode * u[si] + m * omega_c * u[ci]
            end
        end
    end
    return SVector(out)
end
