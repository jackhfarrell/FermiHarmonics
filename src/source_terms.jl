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

Source terms from a BGK-type approximation to the collision integral.  We do not damp 
``a_0`` (density). The momentum modes ``a_1, b_1`` are damped at the momentum-relaxing rate 
``\\gamma_mr``, while higher harmonics are damped at the full scattering rate 
``\\gamma_mr + \\gamma_mc``.
"""
@inline function physical_sources(u, x, t, equations::FermiHarmonics2D)::SVector
    n = length(u)
    out = MVector{n, Float64}(undef)
    @inbounds begin
        # Monopole: no damping (charge conservation)
        out[1] = 0.0
        
        # Dipole: momentum-relaxing scattering only
        if n >= 3
            gamma_mr = equations.gamma_mr
            out[2] = -gamma_mr * u[2]
            out[3] = -gamma_mr * u[3]
        end
        
        # Higher harmonics: full scattering (momentum-relaxing + momentum-conserving)
        if n > 3
            gamma_hi = equations.gamma_mr + equations.gamma_mc
            max_harmonic = (n - 1) ÷ 2
            for m in 2:max_harmonic
                ci = cosine_index(m)
                si = sine_index(m)
                out[ci] = -gamma_hi * u[ci]
                out[si] = -gamma_hi * u[si]
            end
        end
    end
    return SVector(out)
end
