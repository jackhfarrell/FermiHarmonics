using Test
using LinearAlgebra
using StaticArrays
using Random
using ElectronKinetics, Trixi

const TrixiExt = Base.get_extension(ElectronKinetics, :ElectronKineticsTrixiExt)

function reflection_x_matrix(M::Int)
    nvars = 1 + 2 * M
    R = Matrix{Float64}(I, nvars, nvars)
    R[ElectronKinetics.cosine_index(0), ElectronKinetics.cosine_index(0)] = 1.0
    for m in 1:M
        sign_a = isodd(m) ? -1.0 : 1.0
        sign_b = -sign_a
        R[ElectronKinetics.cosine_index(m), ElectronKinetics.cosine_index(m)] = sign_a
        R[ElectronKinetics.sine_index(m), ElectronKinetics.sine_index(m)] = sign_b
    end
    return R
end

function reflection_y_matrix(M::Int)
    nvars = 1 + 2 * M
    R = Matrix{Float64}(I, nvars, nvars)
    for m in 1:M
        R[ElectronKinetics.cosine_index(m), ElectronKinetics.cosine_index(m)] = 1.0
        R[ElectronKinetics.sine_index(m), ElectronKinetics.sine_index(m)] = -1.0
    end
    return R
end

@testset "Symmetry and boundary-operator consistency" begin
    M = 8
    nvars = 1 + 2 * M
    config = SolverConfig(; min_harmonic=M, max_harmonic_auto=M)
    model = KineticModel2D(
        Isotropic2DFermiSurface(),
        HarmonicBasis(M),
        IsotropicHarmonicStreaming(),
        LinearBGKCollision(0.0, ConstantModeRateProfile(0.0)),
    )
    eq, _ = TrixiExt.build_equations(model, config)
    Rx = reflection_x_matrix(M)
    Ry = reflection_y_matrix(M)

    @testset "Streaming reflection identities" begin
        # x-mirror: x -> -x and theta -> pi - theta
        @test opnorm(Rx * eq.Ax * Rx + eq.Ax, Inf) <= 1e-12
        @test opnorm(Rx * eq.Ay * Rx - eq.Ay, Inf) <= 1e-12

        # y-mirror: y -> -y and theta -> -theta
        @test opnorm(Ry * eq.Ax * Ry - eq.Ax, Inf) <= 1e-12
        @test opnorm(Ry * eq.Ay * Ry + eq.Ay, Inf) <= 1e-12
    end

    @testset "Incoming projector idempotence" begin
        unit_n = TrixiExt.unit_normal(SVector(0.6, -0.8))
        P_in = TrixiExt.incoming_projector(eq, unit_n; tol=0.0)
        P_dense = Matrix(P_in)

        @test rank(P_dense) >= 1
        @test norm(P_dense * P_dense - P_dense, Inf) <= 1e-10
    end

    @testset "Specular target invariances" begin
        rng = MersenneTwister(4)
        state = randn(rng, nvars)
        unit_n = TrixiExt.unit_normal(SVector(0.3, 0.7))

        target_plus = zeros(Float64, nvars)
        target_minus = zeros(Float64, nvars)
        target_twice = zeros(Float64, nvars)
        TrixiExt.specular_target!(target_plus, state, unit_n)
        TrixiExt.specular_target!(target_minus, state, -unit_n)
        TrixiExt.specular_target!(target_twice, target_plus, unit_n)

        # Changing n -> -n should leave the geometric mirror operator unchanged.
        @test isapprox(target_plus, target_minus; atol=1e-12, rtol=1e-12)
        # Specular reflection is an involution.
        @test isapprox(target_twice, state; atol=1e-12, rtol=1e-12)
    end

    @testset "Incoming-only BC update" begin
        rng = MersenneTwister(7)
        unit_n = TrixiExt.unit_normal(SVector(-0.4, 0.9))
        P_in = TrixiExt.incoming_projector(eq.Ax, eq.Ay, unit_n; tol=0.0)

        state = randn(rng, nvars)
        target = randn(rng, nvars)
        out = similar(state)
        TrixiExt.apply_projector!(out, state, target, P_in)
        delta = out - state

        expected = Matrix(P_in) * (target - state)
        @test norm(delta - expected, Inf) <= 1e-11
        @test norm((Matrix(P_in) * delta) - delta, Inf) <= 1e-11
    end
end
