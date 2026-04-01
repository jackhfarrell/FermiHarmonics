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

mutable struct NonlinearTimingThreadStats
    flux_ns::UInt64
    flux_calls::UInt64
    bgk_ns::UInt64
    bgk_calls::UInt64
    boundary_ns::UInt64
    boundary_calls::UInt64
    diffuse_ns::UInt64
    diffuse_calls::UInt64
    speed_ns::UInt64
    speed_calls::UInt64
    gradient_ns::UInt64
    gradient_calls::UInt64
end

NonlinearTimingThreadStats() = NonlinearTimingThreadStats(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

const NONLINEAR_TIMING_ENABLED = Ref(false)
const NONLINEAR_TIMING_STATS = Ref(Vector{NonlinearTimingThreadStats}())

@inline function nonlinear_timing_stats()
    stats = NONLINEAR_TIMING_STATS[]
    nthreads = Threads.nthreads()
    if length(stats) != nthreads
        stats = [NonlinearTimingThreadStats() for _ in 1:nthreads]
        NONLINEAR_TIMING_STATS[] = stats
    end
    return stats
end

@inline nonlinear_timing_enabled() = NONLINEAR_TIMING_ENABLED[]

function enable_nonlinear_timing!()
    nonlinear_timing_stats()
    NONLINEAR_TIMING_ENABLED[] = true
    return nothing
end

function disable_nonlinear_timing!()
    NONLINEAR_TIMING_ENABLED[] = false
    return nothing
end

function reset_nonlinear_timing!()
    stats = nonlinear_timing_stats()
    for stat in stats
        stat.flux_ns = 0
        stat.flux_calls = 0
        stat.bgk_ns = 0
        stat.bgk_calls = 0
        stat.boundary_ns = 0
        stat.boundary_calls = 0
        stat.diffuse_ns = 0
        stat.diffuse_calls = 0
        stat.speed_ns = 0
        stat.speed_calls = 0
        stat.gradient_ns = 0
        stat.gradient_calls = 0
    end
    return nothing
end

@inline function record_nonlinear_timing!(category::Symbol, elapsed_ns::UInt64)
    stat = nonlinear_timing_stats()[Threads.threadid()]
    if category === :flux
        stat.flux_ns += elapsed_ns
        stat.flux_calls += 1
    elseif category === :bgk
        stat.bgk_ns += elapsed_ns
        stat.bgk_calls += 1
    elseif category === :boundary
        stat.boundary_ns += elapsed_ns
        stat.boundary_calls += 1
    elseif category === :diffuse
        stat.diffuse_ns += elapsed_ns
        stat.diffuse_calls += 1
    elseif category === :speed
        stat.speed_ns += elapsed_ns
        stat.speed_calls += 1
    elseif category === :gradient
        stat.gradient_ns += elapsed_ns
        stat.gradient_calls += 1
    else
        error("unknown nonlinear timing category: $category")
    end
    return nothing
end

function nonlinear_timing_snapshot()
    totals = Dict(
        :flux => (ns=UInt64(0), calls=UInt64(0)),
        :bgk => (ns=UInt64(0), calls=UInt64(0)),
        :boundary => (ns=UInt64(0), calls=UInt64(0)),
        :diffuse => (ns=UInt64(0), calls=UInt64(0)),
        :speed => (ns=UInt64(0), calls=UInt64(0)),
        :gradient => (ns=UInt64(0), calls=UInt64(0)),
    )
    for stat in nonlinear_timing_stats()
        totals[:flux] = (ns=totals[:flux].ns + stat.flux_ns, calls=totals[:flux].calls + stat.flux_calls)
        totals[:bgk] = (ns=totals[:bgk].ns + stat.bgk_ns, calls=totals[:bgk].calls + stat.bgk_calls)
        totals[:boundary] = (ns=totals[:boundary].ns + stat.boundary_ns, calls=totals[:boundary].calls + stat.boundary_calls)
        totals[:diffuse] = (ns=totals[:diffuse].ns + stat.diffuse_ns, calls=totals[:diffuse].calls + stat.diffuse_calls)
        totals[:speed] = (ns=totals[:speed].ns + stat.speed_ns, calls=totals[:speed].calls + stat.speed_calls)
        totals[:gradient] = (ns=totals[:gradient].ns + stat.gradient_ns, calls=totals[:gradient].calls + stat.gradient_calls)
    end
    return totals
end

function print_nonlinear_timing_summary(io::IO=stdout)
    snapshot = nonlinear_timing_snapshot()
    total_ns = zero(UInt64)
    for key in (:flux, :bgk, :boundary, :diffuse, :speed, :gradient)
        entry = snapshot[key]
        total_ns += entry.ns
        mean_ns = entry.calls == 0 ? 0.0 : Float64(entry.ns) / Float64(entry.calls)
        println(io, rpad(string(key), 10), ": calls=", entry.calls,
                ", total=", round(Float64(entry.ns) * 1.0e-9; digits=4), " s",
                ", mean=", round(mean_ns * 1.0e-6; digits=4), " ms")
    end
    println(io, "total instrumented time: ", round(Float64(total_ns) * 1.0e-9; digits=4), " s")
    return nothing
end

mutable struct BCProjectorCache
    state_buffers::Vector{Vector{Float64}}
    target_buffers::Vector{Vector{Float64}}
    out_buffers::Vector{Vector{Float64}}
    projectors::Dict{Int, SparseMatrixCSC{Float64, Int}}
    nonlinear_faces::Dict{Int, Any}
    initialized::Bool
    nvars::Int
    signature::Tuple{Symbol, Int, Int}
end

BCProjectorCache() = BCProjectorCache(
    [Float64[] for _ in 1:Threads.nthreads()],
    [Float64[] for _ in 1:Threads.nthreads()],
    [Float64[] for _ in 1:Threads.nthreads()],
    Dict{Int, SparseMatrixCSC{Float64, Int}}(),
    Dict{Int, Any}(),
    false,
    0,
    (:unset, 0, 0),
)

struct NonlinearBoundaryFaceData
    unit_normal::SVector{2, Float64}
    incoming_mask::BitVector
    stencil_indices::Matrix{Int}
    stencil_weights::Matrix{Float64}
    projections::Vector{Float64}
    incoming_weight::Float64
    sample_to_harmonics::Union{Nothing, Matrix{Float64}}
end

"""
    AbstractBoundaryCondition

Abstract type for boundary conditions in kinetic transport.

All boundary condition subtypes must implement:
- A `cache::BCProjectorCache` field for thread-safe workspace
- A `tol::Float64` field for numerical tolerance

Subtypes:
- `AbstractWallBC` — Wall-type boundaries (specular/diffuse reflection)
- `AbstractContactBC` — Contact-type boundaries (carrier injection/extraction)
"""
abstract type AbstractBoundaryCondition end

"""
    AbstractWallBC <: AbstractBoundaryCondition

Wall-type boundary conditions. All subtypes have:
- `p_scatter::Float64` — scattering probability

Subtypes:
- `MaxwellWallBC` — specular + diffuse reflection
"""
abstract type AbstractWallBC <: AbstractBoundaryCondition end

"""
    AbstractContactBC <: AbstractBoundaryCondition

Contact-type boundary conditions. All subtypes have:
- `p_ohmic_absorb::Float64` — absorption probability at contact
- `bias::Float64` (OhmicContactBC) or `target_outward_flux::Float64` (CurrentContactBC)

Subtypes:
- `OhmicContactBC` — voltage-controlled contact
- `CurrentContactBC` — current-controlled contact
"""
abstract type AbstractContactBC <: AbstractBoundaryCondition end

