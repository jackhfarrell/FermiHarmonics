using Test
using FermiHarmonics
using Trixi

@testset "FermiHarmonics smoke tests" begin
    params = SolveParams()
    @test params.max_harmonic >= params.min_harmonic
    @test params.max_harmonic_auto >= params.min_harmonic

    @test estimate_max_harmonic(0.0, 0.0; min_harmonic=4, max_harmonic=100) == 100
    @test estimate_max_harmonic(0.0, 500.0; min_harmonic=4, max_harmonic=100) == 4

    eq = FermiHarmonics2D(9; gamma_mr=0.1, gamma_mc=1.0, max_harmonic=4)
    @test typeof(eq) <: Trixi.AbstractEquations{2, 9}
end

include("symmetry_consistency.jl")
