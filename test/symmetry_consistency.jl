using Test
using LinearAlgebra
using StaticArrays
using Random

function reflection_x_matrix(M::Int)
    nvars = 1 + 2 * M
    R = Matrix{Float64}(I, nvars, nvars)
    R[FermiHarmonics.cosine_index(0), FermiHarmonics.cosine_index(0)] = 1.0
    for m in 1:M
        sign_a = isodd(m) ? -1.0 : 1.0
        sign_b = -sign_a
        R[FermiHarmonics.cosine_index(m), FermiHarmonics.cosine_index(m)] = sign_a
        R[FermiHarmonics.sine_index(m), FermiHarmonics.sine_index(m)] = sign_b
    end
    return R
end

function reflection_y_matrix(M::Int)
    nvars = 1 + 2 * M
    R = Matrix{Float64}(I, nvars, nvars)
    for m in 1:M
        R[FermiHarmonics.cosine_index(m), FermiHarmonics.cosine_index(m)] = 1.0
        R[FermiHarmonics.sine_index(m), FermiHarmonics.sine_index(m)] = -1.0
    end
    return R
end

@testset "Symmetry and boundary-operator consistency" begin
    M = 8
    nvars = 1 + 2 * M
    eq = FermiHarmonics2D(nvars; gamma_mr=0.0, gamma_mc=0.0, max_harmonic=M)
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

    @testset "Projector audit and reciprocity" begin
        unit_n = FermiHarmonics.unit_normal(SVector(0.6, -0.8))
        audit = FermiHarmonics.audit_characteristic_projectors(eq.Ax, eq.Ay, unit_n; tol=0.0)

        @test audit.nincoming + audit.noutgoing + audit.ngrazing == nvars
        @test audit.ngrazing >= 1
        @test audit.idempotence_in <= 1e-10
        @test audit.idempotence_out <= 1e-10
        @test audit.idempotence_grazing <= 1e-10
        @test audit.decomposition_error <= 1e-10
        @test audit.orthogonality_error <= 1e-10
        @test audit.opposite_incoming_outgoing_error <= 1e-10
        @test audit.opposite_outgoing_incoming_error <= 1e-10
        @test audit.opposite_grazing_error <= 1e-10
    end

    @testset "Specular target invariances" begin
        rng = MersenneTwister(4)
        state = randn(rng, nvars)
        unit_n = FermiHarmonics.unit_normal(SVector(0.3, 0.7))

        target_plus = zeros(Float64, nvars)
        target_minus = zeros(Float64, nvars)
        target_twice = zeros(Float64, nvars)
        FermiHarmonics.specular_target!(target_plus, state, unit_n)
        FermiHarmonics.specular_target!(target_minus, state, -unit_n)
        FermiHarmonics.specular_target!(target_twice, target_plus, unit_n)

        # Changing n -> -n should leave the geometric mirror operator unchanged.
        @test isapprox(target_plus, target_minus; atol=1e-12, rtol=1e-12)
        # Specular reflection is an involution.
        @test isapprox(target_twice, state; atol=1e-12, rtol=1e-12)
    end

    @testset "Incoming-only BC update" begin
        rng = MersenneTwister(7)
        unit_n = FermiHarmonics.unit_normal(SVector(-0.4, 0.9))
        P_in = FermiHarmonics.incoming_projector(eq.Ax, eq.Ay, unit_n; tol=0.0)
        proj = FermiHarmonics.characteristic_projectors(eq.Ax, eq.Ay, unit_n; tol=0.0)

        state = randn(rng, nvars)
        target = randn(rng, nvars)
        out = similar(state)
        FermiHarmonics.apply_projector!(out, state, target, P_in)
        delta = out - state

        expected = Matrix(P_in) * (target - state)
        @test norm(delta - expected, Inf) <= 1e-11
        @test norm(Matrix(proj.outgoing) * delta, Inf) <= 1e-11
        @test norm(Matrix(proj.grazing) * delta, Inf) <= 1e-11
    end

    @testset "Boundary-node indexing helper" begin
        @test FermiHarmonics.boundary_node_ij(1, 2, 5) == (1, 2)
        @test FermiHarmonics.boundary_node_ij(2, 2, 5) == (5, 2)
        @test FermiHarmonics.boundary_node_ij(3, 4, 5) == (4, 1)
        @test FermiHarmonics.boundary_node_ij(4, 4, 5) == (4, 5)
    end
end
