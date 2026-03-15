"""
    blg_reference_setup(; half_height=0.25, p_scatter=1.0)

Return the shared BLG-oriented dimensionless reference convention used by the
straight-channel nonlinear demos and analysis. The solver API remains
dimensionless; this helper only centralizes the canonical choices.
"""
function blg_reference_setup(; half_height::Real=0.25, p_scatter::Real=1.0)
    channel_length = 1.0
    mu0 = 1.0
    mass = 2.0
    gamma_mr = 0.0
    gamma_mc = 0.0
    vF = zero_state_speed(mu0, mass)
    probe_x = 0.0
    probe_y = 0.0
    probe_offset_fraction = 0.3
    half_length = channel_length / 2

    return (
        convention_name = "blg_reference_dimensionless",
        convention_version = 1,
        channel_length = channel_length,
        half_length = half_length,
        half_height = Float64(half_height),
        mu0 = mu0,
        mass = mass,
        vF = vF,
        gamma_mr = gamma_mr,
        gamma_mc = gamma_mc,
        p_scatter = Float64(p_scatter),
        transport = :parabolic_nonlinear,
        probe_x = probe_x,
        probe_y = probe_y,
        left_probe_x = -probe_offset_fraction * channel_length,
        right_probe_x = probe_offset_fraction * channel_length,
    )
end

