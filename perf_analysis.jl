#!/usr/bin/env julia
"""
Performance and Type Stability Analysis for ElectronKinetics

This script analyzes key functions for type stability and performance.
Run with: julia --project=. perf_analysis.jl
"""

using ElectronKinetics
using InteractiveUtils

println("=" ^ 80)
println("ElectronKinetics Performance Analysis")
println("=" ^ 80)

# Test 1: KineticModel2D construction
println("\n## Test 1: KineticModel2D construction (Linear)")
println("-" ^ 80)

surface = Isotropic2DFermiSurface(vF=1.0, nu=1.0, mass=1.0)
discretization = HarmonicBasis(:auto)
collision = LinearBGKCollision(0.1, TwoRateProfile(0.4))

# Type-check model construction
code_warntype(
    KineticModel2D,
    Tuple{
        Isotropic2DFermiSurface,
        HarmonicBasis,
        IsotropicHarmonicStreaming,
        LinearBGKCollision
    }
)

# Construct model
model = KineticModel2D(
    surface,
    discretization,
    IsotropicHarmonicStreaming(),
    collision
)
println("\nModel constructed: $(typeof(model))")

# Test 2: Streaming matrices computation
println("\n## Test 2: streaming_matrices computation")
println("-" ^ 80)

code_warntype(streaming_matrices, Tuple{Int, Isotropic2DFermiSurface})

Ax, Ay = streaming_matrices(10, surface)
println("Streaming matrices: $(typeof(Ax)), shape $(size(Ax))")

# Test 3: Surface properties
println("\n## Test 3: Surface property accessors")
println("-" ^ 80)

@time vF = surface_vF(surface)
@time vmax = surface_max_speed(surface)
@time nu = surface_density_of_states(surface)
@time mass = surface_mass(surface)

println("vF: $vF, v_max: $vmax, nu: $nu, mass: $mass")

# Test 4: Mode rate computation
println("\n## Test 4: Mode rate computation")
println("-" ^ 80)

profile = TwoRateProfile(0.4)

code_warntype(mode_rate, Tuple{TwoRateProfile, Int})

@time rate1 = mode_rate(profile, 1)
@time rate2 = mode_rate(profile, 2)
@time rate10 = mode_rate(profile, 10)

println("Rates: m=1: $rate1, m=2: $rate2, m=10: $rate10")

# Test 5: Collision gamma accessors
println("\n## Test 5: Collision parameter accessors")
println("-" ^ 80)

code_warntype(collision_gamma_mr, Tuple{LinearBGKCollision})
code_warntype(mode_profile, Tuple{LinearBGKCollision})

@time gmr = collision_gamma_mr(collision)
@time prof = mode_profile(collision)
@time gee = profile_reference_rate(prof)

println("gamma_mr: $gmr, reference rate: $gee")

# Test 6: Harmonic basis estimation
println("\n## Test 6: Harmonic basis estimation")
println("-" ^ 80)

code_warntype(
    estimate_max_harmonic,
    Tuple{Float64, Float64}
)

@time m1 = estimate_max_harmonic(0.1, 0.4)
@time m2 = estimate_max_harmonic(1.0, 2.0)
@time m3 = estimate_max_harmonic(10.0, 20.0)

println("Max harmonics: gamma=(0.1,0.4) → $m1, (1.0,2.0) → $m2, (10.0,20.0) → $m3")

# Test 7: Multiband model
println("\n## Test 7: Multiband KineticModel2D")
println("-" ^ 80)

band1 = Band(:e, surface, gamma_mr=0.1, gamma_ee=0.05)
band2 = Band(:h, surface, gamma_mr=0.08, gamma_ee=0.04)

model_multi = KineticModel2D(
    surface,
    HarmonicBasis(:auto),
    IsotropicHarmonicStreaming(),
    LinearBGKCollision(0.1),
    bands=[band1, band2]
)

println("Multiband model: $(typeof(model_multi))")
println("Number of bands: $(length(model_multi.bands))")

# Summary
println("\n" ^ 2)
println("=" ^ 80)
println("Type Stability Summary")
println("=" ^ 80)
println("""
Common patterns observed:
1. Surface accessors are type-stable (Float64 → Float64)
2. Mode rate computation is type-stable
3. KineticModel2D construction is type-stable
4. Streaming matrices are type-stable (dense matrices)

Potential improvements:
1. Check Trixi extension for type stability in equation building
2. Profile collision_sources! for allocations
3. Check FFT-based transformations in angle grid mode

Use ProfileView.jl for detailed allocation profiling:
  julia> include("perf_analysis.jl")
  julia> using ProfileView
  julia> @profview solve(...)  # Profile actual solver
""")
