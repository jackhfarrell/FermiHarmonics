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
struct FermiHarmonics2D{NVARS} <: Trixi.AbstractEquations{2, NVARS}
    gamma_mr::Float64
    gamma_mc::Float64
    max_speed::Float64
    Ax::Matrix{Float64}
    Ay::Matrix{Float64}
end

"""
    FermiHarmonics2D(nvars; gamma_mr, gamma_mc, max_harmonic=0)

Construct `FermiHarmonics2D`.

Parameters:
- `nvars`: number of state variables, must be odd (`1 + 2M`).
- `gamma_mr`: momentum-relaxing scattering rate.
- `gamma_mc`: momentum-conserving scattering rate.
- `max_harmonic`: optional explicit harmonic cutoff; if set, must satisfy `nvars == 1 + 2*max_harmonic`.

Returns:
- `FermiHarmonics2D{nvars}` equations object.
"""
function FermiHarmonics2D(
    nvars::Integer;
    gamma_mr::Real,
    gamma_mc::Real,
    max_harmonic::Integer = 0,
)
    nvars_int = Int(nvars)
    nvars_int >= 1 || throw(ArgumentError("nvars must be >= 1"))
    isodd(nvars_int) || throw(ArgumentError("nvars must be odd (1 + 2M)"))
    M = (nvars_int - 1) ÷ 2
    if max_harmonic > 0 && max_harmonic != M
        throw(ArgumentError("max_harmonic ($max_harmonic) must match (nvars - 1) ÷ 2 = $M"))
    end

    vF = 1.0
    Ax, Ay = streaming_matrices(M, vF)
    # Canonical LLF speed for this kinetic model: max |v · n| = vF (for unit normals).
    max_speed = vF

    return FermiHarmonics2D{nvars_int}(
        Float64(gamma_mr),
        Float64(gamma_mc),
        max_speed,
        Ax,
        Ay,
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
    print(io, "γ_mc=$(equations.gamma_mc)")
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
        Trixi.summary_footer(io)
    end
end
