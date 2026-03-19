# Linearized Boltzmann equations in harmonic basis



# ======================================================================================================================
# Main Equations Type: FermiHarmonics2D
# ======================================================================================================================

# The FermiHarmonics2D type implements a 2D linearized Boltzmann equation where the
# momentum-dependence is expanded in a basis of circular harmonics. The equations are
# purely advective with source terms from physical scattering.

"""
    FermiHarmonics2D{NVARS} <: Trixi.AbstractEquations{2, NVARS}

Linearized 2D Boltzmann system in harmonic form:
```math
\\partial_t u + A_x \\partial_x u + A_y \\partial_y u = S(u;\\gamma_{mr},\\gamma_{mc})```.
```
"""
struct FermiHarmonics2D{NVARS, TNonlinear} <: Trixi.AbstractEquations{2, NVARS}
    gamma_mr::Float64
    gamma_mc::Float64
    max_speed::Float64
    Ax::Matrix{Float64}
    Ay::Matrix{Float64}
    transport::Symbol
    collision_model::Symbol
    mu0::Float64
    mass::Float64
    electrostatic_coupling::Float64
    theta_oversample::Int
    nonlinear_data::TNonlinear
end

"""
    FermiHarmonics2D(nvars; gamma_mr, gamma_mc, max_harmonic=0, transport=:linear,
                     collision_model=nothing, mu0=nothing, mass=nothing,
                     chi=0.0, theta_oversample=2)

Construct `FermiHarmonics2D`.

Parameters:
- `nvars`: number of state variables, must be odd (`1 + 2M`).
- `gamma_mr`: momentum-relaxing scattering rate.
- `gamma_mc`: momentum-conserving scattering rate.
- `max_harmonic`: optional explicit harmonic cutoff; if set, must satisfy `nvars == 1 + 2*max_harmonic`.
- `transport`: `:linear` or `:parabolic_nonlinear`.
- `collision_model`: defaults to `:linear_mrt` for linear transport and `:exact_bgk` for nonlinear transport.
- `mu0`, `mass`, `theta_oversample`: required nonlinear-transport parameters.
- `chi`: electrostatic coupling in the self-consistent relation `phi = chi * a0 / 2`.

Returns:
- `FermiHarmonics2D{nvars}` equations object.
"""
function FermiHarmonics2D(
    nvars::Integer;
    gamma_mr::Real,
    gamma_mc::Real,
    max_harmonic::Integer = 0,
    transport::Symbol = :linear,
    collision_model::Union{Nothing, Symbol} = nothing,
    mu0::Union{Nothing, Real} = nothing,
    mass::Union{Nothing, Real} = nothing,
    chi::Real = 0.0,
    theta_oversample::Integer = 2,
)
    nvars_int = Int(nvars)
    nvars_int >= 1 || throw(ArgumentError("nvars must be >= 1"))
    isodd(nvars_int) || throw(ArgumentError("nvars must be odd (1 + 2M)"))
    M = (nvars_int - 1) ÷ 2
    if max_harmonic > 0 && max_harmonic != M
        throw(ArgumentError("max_harmonic ($max_harmonic) must match (nvars - 1) ÷ 2 = $M"))
    end

    validate_transport_mode(transport, mu0, mass, theta_oversample)
    if transport === :linear && Float64(chi) != 0.0
        throw(ArgumentError("chi must be 0 for :linear transport"))
    end
    collision_model_value = validate_collision_model(transport, collision_model)

    nonlinear_transport_data = nothing
    mu0_value = isnothing(mu0) ? NaN : Float64(mu0)
    mass_value = isnothing(mass) ? NaN : Float64(mass)
    vF = 1.0
    if transport === :parabolic_nonlinear
        vF = zero_state_speed(mu0_value, mass_value)
        nonlinear_transport_data = create_nonlinear_transport_data(M, Int(theta_oversample))
    end

    Ax, Ay = streaming_matrices(M, vF)
    max_speed = vF

    return FermiHarmonics2D{nvars_int, typeof(nonlinear_transport_data)}(
        Float64(gamma_mr),
        Float64(gamma_mc),
        max_speed,
        Ax,
        Ay,
        transport,
        collision_model_value,
        mu0_value,
        mass_value,
        Float64(chi),
        Int(theta_oversample),
        nonlinear_transport_data,
    )
end


# ======================================================================================================================
# Streaming Matrices and Flux Function
# ======================================================================================================================

# The streaming matrix encodes the linear transport of distributions. In the circular
# harmonic basis, streaming couples adjacent harmonic modes. The flux Jacobians Ax and Ay
# are precomputed and stored in the equations object for efficient evaluation.

"""
    streaming_matrices(M, vF=1.0) -> (Ax, Ay)

Construct ``(1+2M) \\times (1+2M)`` transport matrices with state layout
``\\mathbf{u}=(A_0,A_1,B_1,\\dots,A_M,B_M)^T``. Flux: ``\\mathbf{F}=(n_x A_x+n_y A_y)\\mathbf{u}``
with coupling coefficients ``(A_x)_{A_m,A_{m\\pm1}}=\\tfrac{1}{2}v_F``,
``(A_y)_{A_m,B_{m\\pm1}}=\\pm\\tfrac{1}{2}v_F``, and monopole coupling
``(A_x)_{A_0,A_1}=(A_y)_{A_0,B_1}=v_F``.
"""
function streaming_matrices(M::Int, vF::Float64 = 1.0)
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

"""
    harmonics_flux!(out, state, normal) -> out

Compute flux ``\\mathbf{F}=(n_x A_x+n_y A_y)\\mathbf{u}`` in-place for a given normal direction.
"""
@inline function harmonics_flux!(out::AbstractVector{Float64},
                                 state::AbstractVector{Float64},
                                 normal::SVector{2, Float64})
    n_vars = length(state)
    max_harmonic_local = (n_vars - 1) ÷ 2
    normal_x, normal_y = normal
    vF = 1.0
    @inbounds begin
        out[cosine_index(0)] = (max_harmonic_local >= 1) ?
            (normal_x * vF * state[cosine_index(1)] +
             normal_y * vF * state[sine_index(1)]) : 0.0
        for m in 1:max_harmonic_local
            out[cosine_index(m)] =
                normal_x * (0.5 * vF) * state[cosine_index(m - 1)] +
                (m + 1 <= max_harmonic_local ? normal_x * (0.5 * vF) *
                 state[cosine_index(m + 1)] : 0.0) +
                (m - 1 >= 1 ? normal_y * (-0.5 * vF) *
                 state[sine_index(m - 1)] : 0.0) +
                (m + 1 <= max_harmonic_local ? normal_y * (0.5 * vF) *
                 state[sine_index(m + 1)] : 0.0)
            out[sine_index(m)] =
                (m - 1 >= 1 ? normal_x * (0.5 * vF) *
                 state[sine_index(m - 1)] : 0.0) +
                (m + 1 <= max_harmonic_local ? normal_x * (0.5 * vF) *
                 state[sine_index(m + 1)] : 0.0) +
                normal_y * (0.5 * vF) * state[cosine_index(m - 1)] +
                (m + 1 <= max_harmonic_local ? normal_y * (-0.5 * vF) *
                 state[cosine_index(m + 1)] : 0.0)
        end
    end
    return out
end

"""
    cosine_index(m::Int) -> Int

Index of the cosine harmonic mode ``A_m`` in the state vector.
"""
@inline cosine_index(m::Int) = (m == 0) ? 1 : 2m

"""
    sine_index(m::Int) -> Int

Index of the sine harmonic mode ``B_m`` in the state vector.
"""
@inline sine_index(m::Int) = 2m + 1


# ======================================================================================================================
# Show Methods for Pretty Printing
# ======================================================================================================================

function Base.show(io::IO, equations::FermiHarmonics2D{NVARS}) where {NVARS}
    max_harmonic = (NVARS - 1) ÷ 2
    print(io, "FermiHarmonics2D{$NVARS}(")
    print(io, "max_harmonic=$max_harmonic, ")
    print(io, "γ_mr=$(equations.gamma_mr), ")
    print(io, "γ_mc=$(equations.gamma_mc), ")
    print(io, "transport=$(equations.transport), ")
    print(io, "collision_model=$(equations.collision_model)")
    if transport_is_nonlinear(equations)
        print(io, ", mu0=$(equations.mu0), mass=$(equations.mass), chi=$(equations.electrostatic_coupling), theta_oversample=$(equations.theta_oversample)")
    end
    print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", equations::FermiHarmonics2D{NVARS}) where {NVARS}
    if get(io, :compact, false)
        show(io, equations)
    else
        max_harmonic = (NVARS - 1) ÷ 2
        Trixi.summary_header(io, "FermiHarmonics2D{$NVARS}")
        Trixi.summary_line(io, "max harmonic", max_harmonic)
        Trixi.summary_line(io, "γ_mr (momentum-relaxing)", equations.gamma_mr)
        Trixi.summary_line(io, "γ_mc (momentum-conserving)", equations.gamma_mc)
        Trixi.summary_line(io, "transport", equations.transport)
        Trixi.summary_line(io, "collision model", equations.collision_model)
        if transport_is_nonlinear(equations)
            Trixi.summary_line(io, "mu0", equations.mu0)
            Trixi.summary_line(io, "mass", equations.mass)
            Trixi.summary_line(io, "chi", equations.electrostatic_coupling)
            Trixi.summary_line(io, "theta oversample", equations.theta_oversample)
            Trixi.summary_line(io, "linearized vF", equations.max_speed)
        end
        Trixi.summary_footer(io)
    end
end
