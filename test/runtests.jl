using Test
using FermiHarmonics
using Trixi
using StaticArrays
using LinearAlgebra
using HDF5
using FFTW

@testset "BLG reference convention" begin
    reference = blg_reference_setup()
    @test reference.convention_name == "blg_reference_dimensionless"
    @test reference.channel_length == 1.0
    @test reference.mu0 == 1.0
    @test reference.mass == 2.0
    @test reference.gamma_mr == 0.0
    @test reference.gamma_mc == 0.0
    @test reference.vF ≈ 1.0 atol=1e-12 rtol=1e-12
    @test reference.left_probe_x == -0.3
    @test reference.right_probe_x == 0.3

    geo_path = normpath(joinpath(@__DIR__, "..", "demo", "mesh", "straight_channel.geo"))
    geo_contents = read(geo_path, String)
    @test occursin("length_x = 1.0;", geo_contents)
end

@testset "Mesh-native analysis export" begin
    mesh_path = normpath(joinpath(@__DIR__, "..", "projects", "nonlinearities", "mesh", "tesla_valve.inp"))
    boundary_conditions = Dict(
        :walls => MaxwellWallBC(1.0),
        :inlet => OhmicContactBC(0.05),
        :outlet => OhmicContactBC(-0.05),
    )
    params = SolveParams(;
        polydeg=1,
        tspan_end=0.01,
        residual_tol=1e-3,
        cfl=0.2,
        log_every=10_000,
        min_harmonic=2,
        max_harmonic_auto=4,
    )
    sol, semi = solve(
        mesh_path,
        boundary_conditions,
        params,
        0.0,
        0.5;
        transport=:parabolic_nonlinear,
        max_harmonic=2,
        mu0=1.0,
        mass=2.0,
        chi=0.2,
        name="test_mesh_native_export",
    )

    mesh_data = FermiHarmonics.compute_mesh_native_analysis(sol.u[end], semi; refine=2)
    @test all(isfinite, mesh_data.x)
    @test all(isfinite, mesh_data.y)
    @test all(isfinite, mesh_data.n)
    @test all(isfinite, mesh_data.jx)
    @test all(isfinite, mesh_data.jy)
    @test size(mesh_data.triangles, 2) == 3
    @test minimum(mesh_data.triangles) >= 0
    @test maximum(mesh_data.triangles) < length(mesh_data.x)

    sample_index = findfirst(i -> isfinite(mesh_data.n[i]) && isfinite(mesh_data.jx[i]) && isfinite(mesh_data.jy[i]), eachindex(mesh_data.n))
    @test sample_index !== nothing
    obs = evaluate_observables(sol, semi, mesh_data.x[sample_index], mesh_data.y[sample_index])
    @test obs.in_domain
    @test obs.n ≈ mesh_data.n[sample_index] atol=1e-9 rtol=1e-9
    @test obs.jx ≈ mesh_data.jx[sample_index] atol=1e-9 rtol=1e-9
    @test obs.jy ≈ mesh_data.jy[sample_index] atol=1e-9 rtol=1e-9

    mktempdir() do dir
        mesh_native_path = joinpath(dir, "tesla_mesh_native.h5")
        FermiHarmonics.save_mesh_native_analysis(sol, semi, mesh_native_path; refine=2)
        @test isfile(mesh_native_path)
        h5open(mesh_native_path, "r") do f
            @test haskey(f, "x")
            @test haskey(f, "y")
            @test haskey(f, "triangles")
            @test haskey(f, "n")
            @test haskey(f, "jx")
            @test haskey(f, "jy")
            @test read(attributes(f)["grid_type"]) == "mesh_native_triangles"
            @test read(attributes(f)["time"]) ≈ sol.t[end] atol=1e-12 rtol=1e-12
        end

        reduced_cartesian_path = joinpath(dir, "tesla_cartesian_reduced.h5")
        FermiHarmonics.save_for_analysis(sol, semi, reduced_cartesian_path; nvisnodes=24, observables=[:n, :jx, :jy])
        h5open(reduced_cartesian_path, "r") do f
            @test haskey(f, "n")
            @test haskey(f, "jx")
            @test haskey(f, "jy")
            @test haskey(f, "x")
            @test haskey(f, "y")
            @test haskey(f, "mask")
            @test !haskey(f, "a1")
            @test !haskey(f, "b1")
            @test read(attributes(f)["saved_observables"]) == "n,jx,jy"
        end

        reduced_mesh_native_path = joinpath(dir, "tesla_mesh_native_reduced.h5")
        FermiHarmonics.save_mesh_native_analysis(sol, semi, reduced_mesh_native_path; refine=2, observables=[:n, :jx, :jy])
        h5open(reduced_mesh_native_path, "r") do f
            @test haskey(f, "x")
            @test haskey(f, "y")
            @test haskey(f, "triangles")
            @test haskey(f, "n")
            @test haskey(f, "jx")
            @test haskey(f, "jy")
            @test !haskey(f, "a1")
            @test !haskey(f, "b1")
            @test read(attributes(f)["saved_observables"]) == "n,jx,jy"
        end

        output_png = joinpath(dir, "tesla_mesh_native.png")
        run(`python3 demo/plot_mesh_native_streamlines.py $mesh_native_path --output $output_png --stream-grid 80`)
        @test isfile(output_png)
        @test filesize(output_png) > 0
    end
end

@testset "FermiHarmonics smoke tests" begin
    params = SolveParams()
    @test params.max_harmonic >= params.min_harmonic
    @test params.max_harmonic_auto >= params.min_harmonic

    @test estimate_max_harmonic(0.0, 0.0; min_harmonic=4, max_harmonic=100) == 100
    @test estimate_max_harmonic(0.0, 500.0; min_harmonic=4, max_harmonic=100) == 4

    eq = FermiHarmonics2D(9; gamma_mr=0.1, gamma_mc=1.0, max_harmonic=4)
    @test typeof(eq) <: Trixi.AbstractEquations{2, 9}
end

@testset "Linear multiband harmonics" begin
    bands = [
        BandSpec(
            name=:light,
            vF=1.0,
            nu=1.5,
            mass=1.0,
            charge=-1.0,
            gamma_mr=0.0,
            gamma_mc=0.2,
        ),
        BandSpec(
            name=:heavy,
            vF=1.0,
            nu=2.0,
            mass=3.0,
            charge=1.0,
            gamma_mr=0.0,
            gamma_mc=0.4,
        ),
    ]

    eq = MultiBandFermiHarmonics2D(2; bands=bands, gamma_drag=0.7)
    @test typeof(eq) <: Trixi.AbstractEquations{2, 10}
    @test Trixi.nvariables(eq) == 10
    @test FermiHarmonics.band_count(eq) == 2
    @test FermiHarmonics.band_nvars(eq) == 5
    @test Trixi.varnames(Trixi.cons2cons, eq) ==
          ("light_a0", "light_a1", "light_b1", "light_a2", "light_b2",
           "heavy_a0", "heavy_a1", "heavy_b1", "heavy_a2", "heavy_b2")

    w1 = FermiHarmonics.band_momentum_weight(bands[1])
    w2 = FermiHarmonics.band_momentum_weight(bands[2])

    monopole_state = zeros(10)
    monopole_state[1] = 2.0
    monopole_state[6] = -3.0
    @test FermiHarmonics.physical_sources(monopole_state, nothing, 0.0, eq) ≈ zeros(10) atol=1e-12 rtol=1e-12

    total_momentum_state = zeros(10)
    total_momentum_state[2] = w1
    total_momentum_state[7] = w2
    @test FermiHarmonics.physical_sources(total_momentum_state, nothing, 0.0, eq) ≈ zeros(10) atol=1e-12 rtol=1e-12

    relative_state = zeros(10)
    relative_state[2] = w2
    relative_state[7] = -w1
    relative_source = FermiHarmonics.physical_sources(relative_state, nothing, 0.0, eq)
    @test relative_source[2] ≈ -eq.gamma_drag * relative_state[2] atol=1e-12 rtol=1e-12
    @test relative_source[7] ≈ -eq.gamma_drag * relative_state[7] atol=1e-12 rtol=1e-12

    obs_state = zeros(10)
    obs_state[1] = 1.2
    obs_state[2] = 0.3
    obs_state[3] = -0.1
    obs_state[6] = -0.5
    obs_state[7] = 0.25
    obs_state[8] = 0.4
    obs = FermiHarmonics.multiband_observables(obs_state, eq)
    expected_n = bands[1].nu * obs_state[1] + bands[2].nu * obs_state[6]
    expected_jx = bands[1].charge * bands[1].nu * bands[1].vF * obs_state[2] +
                  bands[2].charge * bands[2].nu * bands[2].vF * obs_state[7]
    expected_jy = bands[1].charge * bands[1].nu * bands[1].vF * obs_state[3] +
                  bands[2].charge * bands[2].nu * bands[2].vF * obs_state[8]
    @test obs.n ≈ expected_n atol=1e-12 rtol=1e-12
    @test obs.jx ≈ expected_jx atol=1e-12 rtol=1e-12
    @test obs.jy ≈ expected_jy atol=1e-12 rtol=1e-12
    @test obs.bands.light.n ≈ bands[1].nu * obs_state[1] atol=1e-12 rtol=1e-12
    @test obs.bands.heavy.jx ≈ bands[2].charge * bands[2].nu * bands[2].vF * obs_state[7] atol=1e-12 rtol=1e-12

    unit_n = SVector(1.0, 0.0)
    P_in_multi = FermiHarmonics.incoming_projector(eq, unit_n)
    wall_out = zeros(10)
    wall_target = zeros(10)
    wall_state = [0.7, 0.2, -0.1, 0.05, 0.03, -0.4, -0.3, 0.25, -0.07, 0.02]
    FermiHarmonics.maxwell_wall!(wall_out, wall_state, unit_n, P_in_multi, 0.35, wall_target, eq)

    expected_wall = zeros(10)
    expected_target = zeros(10)
    for band_index in 1:2
        band_eq = FermiHarmonics2D(5;
            gamma_mr=bands[band_index].gamma_mr,
            gamma_mc=bands[band_index].gamma_mc,
            max_harmonic=2,
        )
        offset = (band_index - 1) * 5
        FermiHarmonics.maxwell_wall!(
            @view(expected_wall[(offset + 1):(offset + 5)]),
            @view(wall_state[(offset + 1):(offset + 5)]),
            unit_n,
            FermiHarmonics.incoming_projector(band_eq, unit_n),
            0.35,
            @view(expected_target[(offset + 1):(offset + 5)]),
            band_eq,
        )
    end
    @test wall_out ≈ expected_wall atol=1e-12 rtol=1e-12

    contact_out = zeros(10)
    contact_target = zeros(10)
    FermiHarmonics.ohmic_contact!(contact_out, wall_state, unit_n, P_in_multi, 1.0, 0.2, contact_target, eq)
    @test contact_target[1] ≈ contact_target[6] atol=1e-12 rtol=1e-12
    @test contact_target[1] ≈ 0.2 atol=1e-12 rtol=1e-12

    mesh_path = normpath(joinpath(@__DIR__, "..", "demo", "mesh", "straight_channel.inp"))
    boundary_conditions = Dict(
        :walls => MaxwellWallBC(1.0),
        :inlet => OhmicContactBC(0.02),
        :outlet => OhmicContactBC(-0.02),
    )
    params = SolveParams(;
        polydeg=1,
        tspan_end=0.05,
        residual_tol=1e-3,
        cfl=0.2,
        log_every=10_000,
        min_harmonic=1,
        max_harmonic_auto=2,
    )

    band1_sol, band1_semi = solve(
        mesh_path,
        boundary_conditions,
        params,
        bands[1].gamma_mr,
        bands[1].gamma_mc;
        max_harmonic=1,
        name="test_multiband_single_light",
    )
    band2_sol, band2_semi = solve(
        mesh_path,
        boundary_conditions,
        params,
        bands[2].gamma_mr,
        bands[2].gamma_mc;
        max_harmonic=1,
        name="test_multiband_single_heavy",
    )
    multi_sol, multi_semi = solve(
        mesh_path,
        boundary_conditions,
        params,
        bands;
        max_harmonic=1,
        gamma_drag=0.0,
        name="test_multiband_uncoupled",
    )

    obs1 = evaluate_observables(band1_sol, band1_semi, 0.0, 0.0)
    obs2 = evaluate_observables(band2_sol, band2_semi, 0.0, 0.0)
    obs_multi = evaluate_observables(multi_sol, multi_semi, 0.0, 0.0)
    @test obs1.in_domain && obs2.in_domain && obs_multi.in_domain
    @test obs_multi.n ≈ bands[1].nu * obs1.a0 + bands[2].nu * obs2.a0 atol=1e-6 rtol=1e-6
    @test obs_multi.jx ≈ bands[1].charge * bands[1].nu * bands[1].vF * obs1.a1 +
                         bands[2].charge * bands[2].nu * bands[2].vF * obs2.a1 atol=1e-6 rtol=1e-6
    @test obs_multi.jy ≈ bands[1].charge * bands[1].nu * bands[1].vF * obs1.b1 +
                         bands[2].charge * bands[2].nu * bands[2].vF * obs2.b1 atol=1e-6 rtol=1e-6

    drag_sol, drag_semi = solve(
        mesh_path,
        boundary_conditions,
        params,
        bands;
        max_harmonic=1,
        gamma_drag=0.3,
        name="test_multiband_drag",
    )
    drag_obs = evaluate_observables(drag_sol, drag_semi, 0.0, 0.0)
    @test drag_obs.in_domain
    @test isfinite(drag_obs.n)
    @test isfinite(drag_obs.jx)
    @test isfinite(drag_obs.jy)
    @test hasproperty(drag_obs.bands, :light)
    @test hasproperty(drag_obs.bands, :heavy)

    grids = FermiHarmonics.compute_analysis_grids(drag_sol.u[end], drag_semi; nvisnodes=12)
    @test haskey(grids.bands, :light)
    @test haskey(grids.bands, :heavy)
    ix = 6
    iy = 6
    @test isfinite(grids.bands[:light].n[ix, iy])
    @test isfinite(grids.bands[:heavy].jx[ix, iy])

    mktempdir() do dir
        cartesian_path = joinpath(dir, "multiband_cartesian.h5")
        FermiHarmonics.save_for_analysis(
            drag_sol,
            drag_semi,
            cartesian_path;
            nvisnodes=12,
            observables=[:n, :jx, :jy, :light_n, :light_jx, :heavy_n, :heavy_jy],
        )
        h5open(cartesian_path, "r") do f
            @test haskey(f, "n")
            @test haskey(f, "jx")
            @test haskey(f, "jy")
            @test haskey(f, "light_n")
            @test haskey(f, "light_jx")
            @test haskey(f, "heavy_n")
            @test haskey(f, "heavy_jy")
            @test read(attributes(f)["saved_observables"]) == "n,jx,jy,light_n,light_jx,heavy_n,heavy_jy"
        end
    end

    @test_throws ArgumentError solve(
        mesh_path,
        boundary_conditions,
        params,
        bands;
        transport=:parabolic_nonlinear,
    )
end

@testset "Quadratic nonlinear transport utilities" begin
    function reference_harmonic_state_to_samples(state, eq)
        ntheta = FermiHarmonics.nonlinear_data(eq).theta_count
        max_harmonic = (length(state) - 1) ÷ 2
        spectrum = zeros(ComplexF64, ntheta)
        spectrum[1] = ComplexF64(0.5 * ntheta * Float64(state[1]), 0.0)
        for m in 1:max_harmonic
            coeff = 0.5 * ComplexF64(Float64(state[FermiHarmonics.cosine_index(m)]), -Float64(state[FermiHarmonics.sine_index(m)]))
            scaled = ntheta * coeff
            spectrum[m + 1] = scaled
            spectrum[ntheta - m + 1] = conj(scaled)
        end
        return ifft(spectrum)
    end

    function reference_harmonic_theta_derivative_to_samples(state, eq)
        ntheta = FermiHarmonics.nonlinear_data(eq).theta_count
        max_harmonic = (length(state) - 1) ÷ 2
        spectrum = zeros(ComplexF64, ntheta)
        for m in 1:max_harmonic
            coeff = 0.5 * ComplexF64(Float64(state[FermiHarmonics.cosine_index(m)]), -Float64(state[FermiHarmonics.sine_index(m)]))
            derivative_coeff = ComplexF64(-imag(coeff) * m, real(coeff) * m)
            scaled = ntheta * derivative_coeff
            spectrum[m + 1] = scaled
            spectrum[ntheta - m + 1] = conj(scaled)
        end
        return ifft(spectrum)
    end

    function reference_electrostatic_force_sources(state, gradients, eq)
        phi = reference_harmonic_state_to_samples(state, eq)
        dtheta = reference_harmonic_theta_derivative_to_samples(state, eq)
        data = FermiHarmonics.nonlinear_data(eq)
        v0 = eq.max_speed
        chi = eq.electrostatic_coupling
        p0 = eq.mass * v0
        grad_phi0_x = 0.5 * Float64(gradients[1][1])
        grad_phi0_y = 0.5 * Float64(gradients[2][1])
        work = Vector{ComplexF64}(undef, length(phi))
        for j in eachindex(work)
            phi_j = real(phi[j])
            dphi_dtheta = real(dtheta[j])
            cos_theta = data.cos_theta[j]
            sin_theta = data.sin_theta[j]
            p_hat_grad_phi0 = cos_theta * grad_phi0_x + sin_theta * grad_phi0_y
            theta_hat_grad_phi0 = -sin_theta * grad_phi0_x + cos_theta * grad_phi0_y
            source = 0.0
            if chi != 0.0
                source -= chi * v0 * p_hat_grad_phi0
                source -= chi * v0 * (0.5 / eq.mu0) * phi_j * p_hat_grad_phi0
                source += (chi / p0) * theta_hat_grad_phi0 * dphi_dtheta
            end
            work[j] = ComplexF64(source, 0.0)
        end
        out = zeros(Float64, length(state))
        FermiHarmonics.samples_to_harmonics!(out, copy(work), eq)
        return out
    end

    function reference_full_nonlinear_force_sources(state, gradients, eq)
        phi = reference_harmonic_state_to_samples(state, eq)
        gradx = reference_harmonic_state_to_samples(gradients[1], eq)
        grady = reference_harmonic_state_to_samples(gradients[2], eq)
        dtheta = reference_harmonic_theta_derivative_to_samples(state, eq)
        data = FermiHarmonics.nonlinear_data(eq)
        v0 = eq.max_speed
        chi = eq.electrostatic_coupling
        p0 = eq.mass * v0
        grad_phi0_x = 0.5 * Float64(gradients[1][1])
        grad_phi0_y = 0.5 * Float64(gradients[2][1])
        work = Vector{ComplexF64}(undef, length(phi))
        @inbounds for j in eachindex(work)
            phi_j = real(phi[j])
            dphi_dtheta = real(dtheta[j])
            dphi_dx = real(gradx[j])
            dphi_dy = real(grady[j])
            cos_theta = data.cos_theta[j]
            sin_theta = data.sin_theta[j]
            p_hat_grad_phi = cos_theta * dphi_dx + sin_theta * dphi_dy
            p_hat_grad_phi0 = cos_theta * grad_phi0_x + sin_theta * grad_phi0_y
            theta_hat_grad_phi0 = -sin_theta * grad_phi0_x + cos_theta * grad_phi0_y
            source = -(v0 * (0.5 / eq.mu0)) * phi_j * p_hat_grad_phi
            if chi != 0.0
                source -= chi * v0 * p_hat_grad_phi0
                source -= chi * v0 * (0.5 / eq.mu0) * phi_j * p_hat_grad_phi0
                source += (chi / p0) * theta_hat_grad_phi0 * dphi_dtheta
            end
            work[j] = ComplexF64(source, 0.0)
        end
        out = zeros(Float64, length(state))
        FermiHarmonics.samples_to_harmonics!(out, copy(work), eq)
        return out
    end

    eq = FermiHarmonics2D(
        9;
        gamma_mr=0.1,
        gamma_mc=1.0,
        max_harmonic=4,
        transport=:parabolic_nonlinear,
        mu0=1.0,
        mass=2.0,
    )
    eq_linear = FermiHarmonics2D(9; gamma_mr=0.1, gamma_mc=1.0, max_harmonic=4)
    @test eq.collision_model === :quadratic_bgk
    @test eq.gamma3 ≈ eq.gamma_mc atol=1e-12 rtol=1e-12
    @test FermiHarmonics.nonlinear_data(eq).theta_count == 16
    @test eq.max_speed ≈ 1.0 atol=1e-12 rtol=1e-12
    @test eq.timestep_speed ≈ 1.0 atol=1e-12 rtol=1e-12
    @test FermiHarmonics.nonlinear_bias_scale(eq) ≈ 1.0 atol=1e-12 rtol=1e-12

    eq_chi = FermiHarmonics2D(
        9;
        gamma_mr=0.1,
        gamma_mc=1.0,
        max_harmonic=4,
        transport=:parabolic_nonlinear,
        mu0=1.0,
        mass=2.0,
        chi=10.0,
    )
    @test eq_chi.max_speed ≈ eq.max_speed atol=1e-12 rtol=1e-12
    @test eq_chi.timestep_speed ≈ 11.0 atol=1e-12 rtol=1e-12
    @test FermiHarmonics.nonlinear_bias_scale(eq_chi) ≈ 11.0 atol=1e-12 rtol=1e-12
    @test FermiHarmonics.nonlinear_electrochemical_bias(0.1, eq_chi) ≈ 1.1 atol=1e-12 rtol=1e-12

    cache = FermiHarmonics.get_nonlinear_cache(eq)
    sample_state = zeros(Float64, 9)
    sample_state[1] = 0.08
    sample_state[2] = 0.03
    sample_state[3] = -0.02
    sample_state[4] = 0.01
    sample_state[5] = 0.006
    reference_samples = reference_harmonic_state_to_samples(sample_state, eq)
    FermiHarmonics.harmonic_state_to_spectrum!(cache.spectrum, sample_state, eq)
    FermiHarmonics.harmonic_spectrum_to_samples!(cache.samples, cache.spectrum, eq)
    @test cache.samples ≈ reference_samples atol=1e-12 rtol=1e-12
    reference_derivative = reference_harmonic_theta_derivative_to_samples(sample_state, eq)
    FermiHarmonics.harmonic_spectrum_to_theta_derivative_samples!(cache.theta_derivative_samples, cache.spectrum, eq)
    @test cache.theta_derivative_samples ≈ reference_derivative atol=1e-12 rtol=1e-12

    gradients = (copy(sample_state), -0.5 .* sample_state)
    gradients[1][1] = 0.04
    gradients[2][1] = -0.03
    source_ref = reference_electrostatic_force_sources(sample_state, gradients, eq_chi)
    source_new = similar(source_ref)
    FermiHarmonics.electrostatic_force_sources!(source_new, sample_state, gradients, eq_chi)
    @test source_new ≈ source_ref atol=1e-12 rtol=1e-12
    source_reference_helper = similar(source_ref)
    FermiHarmonics.electrostatic_force_sources_reference!(source_reference_helper, sample_state, gradients, eq_chi)
    @test source_reference_helper ≈ reference_full_nonlinear_force_sources(sample_state, gradients, eq_chi) atol=1e-12 rtol=1e-12
    alloc_source_sparse = @allocated FermiHarmonics.electrostatic_force_sources!(source_new, sample_state, gradients, eq_chi)
    alloc_source_reference = @allocated FermiHarmonics.electrostatic_force_sources_reference!(source_reference_helper, sample_state, gradients, eq_chi)
    @test alloc_source_sparse <= alloc_source_reference
    FermiHarmonics.prepare_harmonic_gradient_theta_work!(
        cache.samples,
        cache.theta_derivative_samples,
        cache.gradx_samples,
        cache.grady_samples,
        cache.spectrum,
        sample_state,
        gradients,
        eq,
    )
    alloc_theta = @allocated FermiHarmonics.prepare_harmonic_gradient_theta_work!(
        cache.samples,
        cache.theta_derivative_samples,
        cache.gradx_samples,
        cache.grady_samples,
        cache.spectrum,
        sample_state,
        gradients,
        eq,
    )
    @test alloc_theta <= 64

    state = zeros(Float64, 9)
    state[1] = 0.04
    state[2] = 0.02
    state[3] = -0.01
    state[4] = 0.008
    state[5] = -0.006

    flux_zero_x = Trixi.flux(zeros(9), 1, eq)
    flux_zero_y = Trixi.flux(zeros(9), 2, eq)
    @test flux_zero_x ≈ zeros(9)
    @test flux_zero_y ≈ zeros(9)

    flux_x = Trixi.flux(state, 1, eq)
    flux_linear_x = Trixi.flux(state, 1, eq_linear)
    @test norm(flux_x - flux_linear_x) < 2.0e-3

    expected_density = eq.mass * (eq.mu0 + 0.5 * state[1]) / (2.0 * pi)
    @test FermiHarmonics.nonlinear_density(state, eq) ≈ expected_density atol=1e-12 rtol=1e-12
    @test FermiHarmonics.derived_harmonics(state, eq) == (state[1], state[2], state[3])

    mu_target = 1.08
    velocity_target = SVector(0.05, -0.02)
    equilibrium_state = zeros(Float64, 9)
    FermiHarmonics.local_equilibrium_state!(equilibrium_state, mu_target, velocity_target, eq)
    recovered_mu, recovered_velocity = FermiHarmonics.recover_mu_u(equilibrium_state, eq)
    @test recovered_mu ≈ mu_target atol=1e-12 rtol=1e-12
    @test recovered_velocity ≈ velocity_target atol=1e-12 rtol=1e-12
    @test equilibrium_state[FermiHarmonics.cosine_index(2)] ≈
          0.5 * eq.mass * (velocity_target[1]^2 - velocity_target[2]^2) atol=1e-12 rtol=1e-12
    @test equilibrium_state[FermiHarmonics.sine_index(2)] ≈
          eq.mass * velocity_target[1] * velocity_target[2] atol=1e-12 rtol=1e-12

    reconstructed_state = zeros(Float64, 9)
    FermiHarmonics.local_equilibrium_state!(reconstructed_state, equilibrium_state, eq)
    @test reconstructed_state ≈ equilibrium_state atol=1e-12 rtol=1e-12

    perturbed_state = copy(equilibrium_state)
    perturbed_state[4] += 0.01
    perturbed_state[5] -= 0.008
    target_density = FermiHarmonics.nonlinear_density(perturbed_state, eq)
    target_current = FermiHarmonics.nonlinear_current(perturbed_state, eq)
    recovered_mu_perturbed, recovered_velocity_perturbed = FermiHarmonics.recover_mu_u(perturbed_state, eq)
    matched_state = zeros(Float64, 9)
    FermiHarmonics.local_equilibrium_state!(matched_state, recovered_mu_perturbed, recovered_velocity_perturbed, eq)
    @test FermiHarmonics.nonlinear_density(matched_state, eq) ≈ target_density atol=1e-12 rtol=1e-12
    @test collect(FermiHarmonics.nonlinear_current(matched_state, eq)) ≈ collect(target_current) atol=5e-4 rtol=5e-3

    high_eq = FermiHarmonics2D(
        9;
        gamma_mr=0.0,
        gamma_mc=20.0,
        max_harmonic=4,
        transport=:parabolic_nonlinear,
        mu0=1.0,
        mass=2.0,
    )
    high_state = copy(equilibrium_state)
    high_state[4] += 1.0e-3
    high_state[5] -= 7.0e-4
    high_mu, high_velocity = FermiHarmonics.recover_mu_u(high_state, high_eq)
    @test isfinite(high_mu)
    @test all(isfinite, high_velocity)
    high_source = FermiHarmonics.physical_sources(high_state, nothing, 0.0, high_eq)
    @test all(isfinite, high_source)

    @test_throws ArgumentError FermiHarmonics2D(
        9;
        gamma_mr=0.1,
        gamma_mc=0.2,
        max_harmonic=4,
        transport=:parabolic_nonlinear,
        collision_model=:exact_bgk,
        mu0=1.0,
        mass=2.0,
    )
    @test_throws ArgumentError FermiHarmonics2D(
        3;
        gamma_mr=0.1,
        gamma_mc=0.2,
        transport=:linear,
        chi=0.1,
    )
    @test_throws ArgumentError FermiHarmonics2D(
        9;
        gamma_mr=0.1,
        gamma_mc=0.2,
        gamma3=-0.1,
        max_harmonic=4,
        transport=:parabolic_nonlinear,
        mu0=1.0,
        mass=2.0,
    )
end

@testset "Exact nonlinear angle reference" begin
    eq = FermiAngles2D(
        32;
        gamma_mr=0.1,
        gamma_mc=1.0,
        mu0=1.0,
        mass=2.0,
    )
    eq_chi = FermiAngles2D(
        32;
        gamma_mr=0.1,
        gamma_mc=1.0,
        mu0=1.0,
        mass=2.0,
        chi=10.0,
    )
    eq_fast = FermiAngles2D(
        32;
        gamma_mr=0.1,
        gamma_mc=1.0,
        collision_model=:two_rate_bgk,
        mu0=1.0,
        mass=2.0,
    )
    @test eq.collision_model === :exact_bgk
    @test eq_fast.collision_model === :two_rate_bgk
    @test eq.max_speed ≈ 1.0 atol=1e-12 rtol=1e-12
    @test eq.timestep_speed ≈ 1.0 atol=1e-12 rtol=1e-12
    @test eq_chi.max_speed ≈ eq.max_speed atol=1e-12 rtol=1e-12
    @test eq_chi.timestep_speed ≈ 11.0 atol=1e-12 rtol=1e-12

    data = FermiHarmonics.nonlinear_data(eq)
    state = @. 0.08 * cos(data.theta) - 0.05 * sin(data.theta) + 0.03 * cos(2.0 * data.theta)

    flux_x = collect(Trixi.flux(state, 1, eq))
    @test flux_x ≈ [data.cos_theta[j] * FermiHarmonics.parabolic_shifted_flux(state[j], eq) for j in eachindex(state)] atol=1e-12 rtol=1e-12

    a0, a1, b1 = FermiHarmonics.derived_harmonics(state, eq)
    @test a0 ≈ 0.0 atol=1e-10 rtol=1e-10
    @test a1 ≈ 0.08 atol=1e-10 rtol=1e-10
    @test b1 ≈ -0.05 atol=1e-10 rtol=1e-10

    weak_state = @. 0.02 * cos(data.theta) - 0.01 * sin(data.theta) + 0.005 * cos(2.0 * data.theta)
    source_exact = FermiHarmonics.physical_sources(weak_state, nothing, 0.0, eq)
    source_fast = FermiHarmonics.physical_sources(weak_state, nothing, 0.0, eq_fast)
    @test all(isfinite, source_fast)
    @test norm(source_exact - source_fast) < 1.0e-3
    alloc_exact = @allocated FermiHarmonics.physical_sources(weak_state, nothing, 0.0, eq)
    alloc_fast = @allocated FermiHarmonics.physical_sources(weak_state, nothing, 0.0, eq_fast)
    @test alloc_fast <= alloc_exact

    @test_throws ArgumentError FermiAngles2D(
        7;
        gamma_mr=0.1,
        gamma_mc=0.2,
        mu0=1.0,
        mass=1.0,
    )
    @test_throws ArgumentError FermiAngles2D(
        10;
        gamma_mr=0.1,
        gamma_mc=0.2,
        collision_model=:quadratic_bgk,
        mu0=1.0,
        mass=1.0,
    )
    @test_throws DomainError Trixi.flux(fill(-2.5, 32), 1, FermiAngles2D(
        32;
        gamma_mr=0.1,
        gamma_mc=0.2,
        mu0=1.0,
        mass=1.0,
    ))
end

@testset "Electrostatic self-consistent force" begin
    function expected_electrostatic_force_source(state, gradients, eq)
        ntheta = FermiHarmonics.nonlinear_data(eq).theta_count
        phi = Vector{ComplexF64}(undef, ntheta)
        dtheta = Vector{ComplexF64}(undef, ntheta)
        work = Vector{ComplexF64}(undef, ntheta)
        FermiHarmonics.harmonic_state_to_samples!(phi, state, eq)
        FermiHarmonics.harmonic_theta_derivative_to_samples!(dtheta, state, eq)
        data = FermiHarmonics.nonlinear_data(eq)
        v0 = eq.max_speed
        p0 = eq.mass * v0
        grad_phi0_x = 0.5 * gradients[1][1]
        grad_phi0_y = 0.5 * gradients[2][1]
        @inbounds for j in eachindex(work)
            phi_j = real(phi[j])
            p_hat_grad_phi0 = data.cos_theta[j] * grad_phi0_x + data.sin_theta[j] * grad_phi0_y
            theta_hat_grad_phi0 = -data.sin_theta[j] * grad_phi0_x + data.cos_theta[j] * grad_phi0_y
            source = 0.0
            if eq.electrostatic_coupling != 0.0
                source -= eq.electrostatic_coupling * v0 * p_hat_grad_phi0
                source -= eq.electrostatic_coupling * v0 * (0.5 / eq.mu0) * phi_j * p_hat_grad_phi0
                source += (eq.electrostatic_coupling / p0) * theta_hat_grad_phi0 * real(dtheta[j])
            end
            work[j] = ComplexF64(source, 0.0)
        end
        out = zeros(Float64, length(state))
        FermiHarmonics.samples_to_harmonics!(out, work, eq)
        return out
    end

    eq = FermiHarmonics2D(
        9;
        gamma_mr=0.0,
        gamma_mc=0.3,
        max_harmonic=4,
        transport=:parabolic_nonlinear,
        mu0=1.0,
        mass=2.0,
        chi=0.6,
    )
    eq_parabolic = FermiHarmonics.ElectrostaticGradientEquation2D(eq)

    uniform_state = zeros(Float64, 9)
    zero_gradients = (zeros(9), zeros(9))
    uniform_force = FermiHarmonics.source_terms(uniform_state, zero_gradients, nothing, 0.0, eq_parabolic)
    @test uniform_force ≈ zeros(9) atol=1.0e-12 rtol=1.0e-12

    driven_state = zeros(Float64, 9)
    driven_state[1] = 0.04
    driven_state[2] = 0.03
    driven_state[3] = -0.02
    driven_state[4] = 0.01
    gradients = (vcat(0.8, zeros(8)), vcat(-0.5, zeros(8)))
    force_source = FermiHarmonics.source_terms(driven_state, gradients, nothing, 0.0, eq_parabolic)
    @test all(isfinite, force_source)
    @test norm(force_source) > 0.0
    @test collect(force_source) ≈ expected_electrostatic_force_source(driven_state, gradients, eq) atol=1.0e-12 rtol=1.0e-12

    eq_zero = FermiHarmonics2D(
        9;
        gamma_mr=0.0,
        gamma_mc=0.3,
        max_harmonic=4,
        transport=:parabolic_nonlinear,
        mu0=1.0,
        mass=2.0,
        chi=0.0,
    )
    eq_zero_parabolic = FermiHarmonics.ElectrostaticGradientEquation2D(eq_zero)
    force_zero = FermiHarmonics.source_terms(driven_state, gradients, nothing, 0.0, eq_zero_parabolic)
    @test all(isfinite, force_zero)
    @test force_zero ≈ zeros(9) atol=1.0e-12 rtol=1.0e-12
end

@testset "Nonlinear BGK source terms" begin
    eq_bgk = FermiHarmonics2D(
        9;
        gamma_mr=0.0,
        gamma_mc=0.8,
        gamma3=0.05,
        max_harmonic=4,
        transport=:parabolic_nonlinear,
        mu0=1.0,
        mass=2.0,
    )
    equilibrium_state = zeros(Float64, 9)
    FermiHarmonics.local_equilibrium_state!(equilibrium_state, 1.06, SVector(0.05, -0.02), eq_bgk)
    source_eq = FermiHarmonics.physical_sources(equilibrium_state, nothing, 0.0, eq_bgk)
    @test source_eq ≈ zeros(9) atol=5e-12 rtol=5e-12
    recovered_mu, recovered_velocity = FermiHarmonics.recover_mu_u(equilibrium_state, eq_bgk)
    drift_equilibrium = zeros(Float64, 9)
    FermiHarmonics.local_equilibrium_state!(drift_equilibrium, recovered_mu, recovered_velocity, eq_bgk)
    @test FermiHarmonics.nonlinear_density(drift_equilibrium, eq_bgk) ≈
          FermiHarmonics.nonlinear_density(equilibrium_state, eq_bgk) atol=1e-12 rtol=1e-12
    @test collect(FermiHarmonics.nonlinear_current(drift_equilibrium, eq_bgk)) ≈
          collect(FermiHarmonics.nonlinear_current(equilibrium_state, eq_bgk)) atol=1e-11 rtol=1e-11

    perturbed = copy(equilibrium_state)
    perturbed[5] += 0.02
    source_perturbed = FermiHarmonics.physical_sources(perturbed, nothing, 0.0, eq_bgk)
    @test norm(source_perturbed) > 0.0
    @test all(isfinite, source_perturbed)
    @test source_perturbed[1] ≈ 0.0 atol=1e-12 rtol=1e-12
    @test source_perturbed[2] ≈ 0.0 atol=1e-12 rtol=1e-12
    @test source_perturbed[3] ≈ 0.0 atol=1e-12 rtol=1e-12

    even_mode_state = zeros(Float64, 9)
    even_mode_state[FermiHarmonics.cosine_index(2)] = 0.03
    even_mode_source = FermiHarmonics.physical_sources(even_mode_state, nothing, 0.0, eq_bgk)
    @test even_mode_source[FermiHarmonics.cosine_index(2)] ≈ -eq_bgk.gamma_mc * even_mode_state[FermiHarmonics.cosine_index(2)] atol=1e-12 rtol=1e-12

    odd_mode_state = zeros(Float64, 9)
    odd_mode_state[FermiHarmonics.cosine_index(3)] = 0.03
    odd_mode_source = FermiHarmonics.physical_sources(odd_mode_state, nothing, 0.0, eq_bgk)
    @test odd_mode_source[FermiHarmonics.cosine_index(3)] ≈
          -min(eq_bgk.gamma_mc, eq_bgk.gamma3 * 3^4) * odd_mode_state[FermiHarmonics.cosine_index(3)] atol=1e-12 rtol=1e-12

    eq_two_rate = FermiHarmonics2D(
        9;
        gamma_mr=0.2,
        gamma_mc=0.5,
        max_harmonic=4,
        transport=:parabolic_nonlinear,
        mu0=1.0,
        mass=2.0,
    )
    two_rate_state = zeros(Float64, 9)
    FermiHarmonics.local_equilibrium_state!(two_rate_state, 1.04, SVector(0.04, 0.01), eq_two_rate)
    two_rate_source = FermiHarmonics.physical_sources(two_rate_state, nothing, 0.0, eq_two_rate)
    @test norm(two_rate_source) > 0.0
    @test two_rate_source[2] ≈ -eq_two_rate.gamma_mr * two_rate_state[2] atol=1e-12 rtol=1e-12
    @test two_rate_source[3] ≈ -eq_two_rate.gamma_mr * two_rate_state[3] atol=1e-12 rtol=1e-12
end

@testset "Quadratic nonlinear regression vs exact reference" begin
    eq_quad = FermiHarmonics2D(
        9;
        gamma_mr=0.1,
        gamma_mc=0.4,
        max_harmonic=4,
        transport=:parabolic_nonlinear,
        mu0=1.0,
        mass=2.0,
    )
    ntheta = FermiHarmonics.nonlinear_data(eq_quad).theta_count
    eq_exact = FermiAngles2D(
        ntheta;
        gamma_mr=0.1,
        gamma_mc=0.4,
        mu0=1.0,
        mass=2.0,
    )

    state = zeros(Float64, 9)
    state[1] = 0.004
    state[2] = 0.0025
    state[3] = -0.0015
    state[4] = 0.001
    state[5] = 0.0008

    samples = Vector{ComplexF64}(undef, ntheta)
    FermiHarmonics.harmonic_state_to_samples!(samples, state, eq_quad)
    sample_state = real.(samples)

    flux_quad = zeros(Float64, 9)
    FermiHarmonics.nonlinear_flux!(flux_quad, state, SVector(1.0, 0.0), eq_quad)
    flux_exact_samples = zeros(Float64, ntheta)
    FermiHarmonics.nonlinear_flux!(flux_exact_samples, sample_state, SVector(1.0, 0.0), eq_exact)
    flux_exact_harmonics = zeros(Float64, 9)
    FermiHarmonics.samples_to_harmonics!(flux_exact_harmonics, ComplexF64.(collect(flux_exact_samples)), eq_quad)
    @test norm(flux_quad - flux_exact_harmonics) < 5.0e-5

    equilibrium_quad = zeros(Float64, 9)
    FermiHarmonics.local_equilibrium_state!(equilibrium_quad, 1.03, SVector(0.03, -0.02), eq_quad)
    equilibrium_exact_samples = zeros(Float64, ntheta)
    FermiHarmonics.local_equilibrium_state!(equilibrium_exact_samples, 1.03, SVector(0.03, -0.02), eq_exact)
    equilibrium_exact_harmonics = zeros(Float64, 9)
    FermiHarmonics.samples_to_harmonics!(equilibrium_exact_harmonics, ComplexF64.(collect(equilibrium_exact_samples)), eq_quad)
    @test equilibrium_quad ≈ equilibrium_exact_harmonics atol=2.0e-4 rtol=2.0e-3

    reference_state = zeros(Float64, 9)
    FermiHarmonics.local_equilibrium_state!(reference_state, 1.03, SVector(0.03, -0.02), eq_quad)
    reference_samples = Vector{ComplexF64}(undef, ntheta)
    FermiHarmonics.harmonic_state_to_samples!(reference_samples, reference_state, eq_quad)
    reference_sample_state = real.(reference_samples)

    recovered_mu_quad, recovered_velocity_quad = FermiHarmonics.recover_mu_u(reference_state, eq_quad)
    recovered_mu_exact, recovered_velocity_exact = FermiHarmonics.recover_mu_u(reference_sample_state, eq_exact)
    @test recovered_mu_quad ≈ recovered_mu_exact atol=2.0e-4 rtol=2.0e-3

    recovered_quad_state = zeros(Float64, 9)
    FermiHarmonics.local_equilibrium_state!(recovered_quad_state, recovered_mu_quad, recovered_velocity_quad, eq_quad)
    recovered_exact_samples = zeros(Float64, ntheta)
    FermiHarmonics.local_equilibrium_state!(recovered_exact_samples, recovered_mu_exact, recovered_velocity_exact, eq_exact)
    recovered_exact_harmonics = zeros(Float64, 9)
    FermiHarmonics.samples_to_harmonics!(recovered_exact_harmonics, ComplexF64.(collect(recovered_exact_samples)), eq_quad)
    @test FermiHarmonics.nonlinear_density(recovered_quad_state, eq_quad) ≈
          FermiHarmonics.nonlinear_density(reference_state, eq_quad) atol=1.0e-12 rtol=1.0e-12
    @test collect(FermiHarmonics.nonlinear_current(recovered_quad_state, eq_quad)) ≈
          collect(FermiHarmonics.nonlinear_current(reference_state, eq_quad)) atol=1.0e-11 rtol=1.0e-11
    @test FermiHarmonics.nonlinear_density(recovered_exact_samples, eq_exact) ≈
          FermiHarmonics.nonlinear_density(reference_sample_state, eq_exact) atol=1.0e-12 rtol=1.0e-12
    @test all(isfinite, recovered_exact_samples)
    @test all(isfinite, recovered_velocity_exact)
end

@testset "Nonlinear boundary conditions" begin
    function reference_harmonic_state_to_samples(state, eq)
        ntheta = FermiHarmonics.nonlinear_data(eq).theta_count
        max_harmonic = (length(state) - 1) ÷ 2
        spectrum = zeros(ComplexF64, ntheta)
        spectrum[1] = ComplexF64(0.5 * ntheta * Float64(state[1]), 0.0)
        for m in 1:max_harmonic
            coeff = 0.5 * ComplexF64(Float64(state[FermiHarmonics.cosine_index(m)]), -Float64(state[FermiHarmonics.sine_index(m)]))
            scaled = ntheta * coeff
            spectrum[m + 1] = scaled
            spectrum[ntheta - m + 1] = conj(scaled)
        end
        return ifft(spectrum)
    end

    function reference_nonlinear_boundary_samples(state, unit_normal, incoming_value, specular_weight, eq, tol)
        state_samples = reference_harmonic_state_to_samples(state, eq)
        result_samples = copy(state_samples)
        face_data = FermiHarmonics.build_nonlinear_face_data(eq, unit_normal, tol)
        diffuse_weight = 1.0 - specular_weight
        for j in eachindex(result_samples)
            if face_data.incoming_mask[j]
                specular_value = specular_weight > 0.0 ? real(FermiHarmonics.apply_specular_stencil(state_samples, face_data, j)) : 0.0
                incoming_sample = diffuse_weight * incoming_value + specular_weight * specular_value
                result_samples[j] = ComplexF64(incoming_sample, 0.0)
            end
        end
        out = zeros(Float64, length(state))
        FermiHarmonics.samples_to_harmonics!(out, result_samples, eq)
        return out
    end

    function reference_nonlinear_boundary_flux(state, normal, unit_normal, incoming_value, specular_weight, eq, tol)
        state_samples = reference_harmonic_state_to_samples(state, eq)
        face_data = FermiHarmonics.build_nonlinear_face_data(eq, unit_normal, tol)
        diffuse_weight = 1.0 - specular_weight
        flux_samples = Vector{ComplexF64}(undef, length(state_samples))
        scale = hypot(normal[1], normal[2])
        for j in eachindex(state_samples)
            phi_trace = real(state_samples[j])
            if face_data.incoming_mask[j]
                specular_value = specular_weight > 0.0 ? real(FermiHarmonics.apply_specular_stencil(state_samples, face_data, j)) : 0.0
                phi_trace = diffuse_weight * incoming_value + specular_weight * specular_value
            end
            directional = scale * face_data.projections[j]
            flux_samples[j] = ComplexF64(
                directional * FermiHarmonics.quadratic_shifted_flux(phi_trace, eq),
                0.0,
            )
        end
        out = zeros(Float64, length(state))
        FermiHarmonics.samples_to_harmonics!(out, flux_samples, eq)
        return out
    end

    eq = FermiHarmonics2D(
        9;
        gamma_mr=0.0,
        gamma_mc=0.0,
        max_harmonic=4,
        transport=:parabolic_nonlinear,
        mu0=1.0,
        mass=2.0,
    )
    eq_linear = FermiHarmonics2D(9; gamma_mr=0.0, gamma_mc=0.0, max_harmonic=4)
    unit_normal = SVector(1.0, 0.0)
    P_in = FermiHarmonics.incoming_projector(eq, unit_normal; tol=1.0e-12)
    P_in_linear = FermiHarmonics.incoming_projector(eq_linear, unit_normal; tol=1.0e-12)
    @test Matrix(P_in) ≈ Matrix(P_in_linear) atol=1.0e-12 rtol=1.0e-12

    state = zeros(Float64, 9)
    out = similar(state)
    scratch = similar(state)

    FermiHarmonics.maxwell_wall!(out, state, unit_normal, P_in, 1.0, scratch)
    @test out ≈ zeros(9)

    FermiHarmonics.ohmic_contact!(out, state, unit_normal, P_in, 1.0, 1.0, scratch)
    out_a0, out_a1, out_b1 = FermiHarmonics.derived_harmonics(out, eq)
    @test out_a0 > 0.0
    @test out_a1 < 0.0
    @test isfinite(out_b1)

    out_top = similar(state)
    out_bottom = similar(state)
    P_top = FermiHarmonics.incoming_projector(eq, SVector(0.0, -1.0); tol=1.0e-12)
    P_bottom = FermiHarmonics.incoming_projector(eq, SVector(0.0, 1.0); tol=1.0e-12)
    FermiHarmonics.ohmic_contact!(out_top, state, SVector(0.0, -1.0), P_top, 1.0, 1.0, scratch)
    FermiHarmonics.ohmic_contact!(out_bottom, state, SVector(0.0, 1.0), P_bottom, 1.0, 1.0, scratch)
    top_a0, top_a1, top_b1 = FermiHarmonics.derived_harmonics(out_top, eq)
    bottom_a0, bottom_a1, bottom_b1 = FermiHarmonics.derived_harmonics(out_bottom, eq)
    @test top_a0 ≈ bottom_a0 atol=1e-12 rtol=1e-12
    @test top_a1 ≈ bottom_a1 atol=1e-12 rtol=1e-12
    @test top_b1 ≈ -bottom_b1 atol=1e-12 rtol=1e-12

    nonlinear_state = zeros(Float64, 9)
    nonlinear_state[1] = 0.08
    nonlinear_state[2] = 0.04
    nonlinear_state[3] = -0.015
    nonlinear_state[4] = 0.01
    nonlinear_bc = OhmicContactBC(0.12)
    linear_bc = OhmicContactBC(0.12)
    trace_flux = (u_inner, u_outer, normal_direction, equations) -> SVector{length(u_outer), Float64}(u_outer)
    nonlinear_trace = nonlinear_bc(nonlinear_state, [1.0, 0.0], nothing, 0.0, trace_flux, eq)
    linear_trace = linear_bc(nonlinear_state, [1.0, 0.0], nothing, 0.0, trace_flux, eq_linear)
    expected_nonlinear_trace = similar(nonlinear_state)
    FermiHarmonics.nonlinear_ohmic_contact!(
        expected_nonlinear_trace,
        nonlinear_state,
        SVector(1.0, 0.0),
        1.0,
        0.12,
        similar(nonlinear_state),
        eq,
        1.0e-12,
    )
    @test collect(nonlinear_trace) ≈ expected_nonlinear_trace atol=1.0e-12 rtol=1.0e-12
    @test norm(collect(nonlinear_trace) - collect(linear_trace)) > 1.0e-6

    wall_trace = zeros(Float64, 9)
    FermiHarmonics.nonlinear_maxwell_wall!(
        wall_trace,
        nonlinear_state,
        unit_normal,
        0.7,
        similar(nonlinear_state),
        eq,
        1.0e-12,
    )
    state_samples = reference_harmonic_state_to_samples(nonlinear_state, eq)
    wall_face_data = FermiHarmonics.build_nonlinear_face_data(eq, unit_normal, 1.0e-12)
    wall_incoming = FermiHarmonics.nonlinear_diffuse_incoming_value(state_samples, wall_face_data, eq, 1.0e-12)
    wall_trace_ref = reference_nonlinear_boundary_samples(
        nonlinear_state,
        unit_normal,
        wall_incoming,
        0.3,
        eq,
        1.0e-12,
    )
    @test wall_trace ≈ wall_trace_ref atol=1.0e-12 rtol=1.0e-12

    wall_flux = zeros(Float64, 9)
    FermiHarmonics.nonlinear_maxwell_wall_flux!(
        wall_flux,
        nonlinear_state,
        unit_normal,
        unit_normal,
        0.7,
        similar(nonlinear_state),
        eq,
        1.0e-12,
    )
    wall_flux_ref = reference_nonlinear_boundary_flux(
        nonlinear_state,
        unit_normal,
        unit_normal,
        wall_incoming,
        0.3,
        eq,
        1.0e-12,
    )
    @test wall_flux ≈ wall_flux_ref atol=1.0e-12 rtol=1.0e-12

    chi_eq = FermiHarmonics2D(
        9;
        gamma_mr=0.0,
        gamma_mc=0.5,
        max_harmonic=4,
        transport=:parabolic_nonlinear,
        mu0=1.0,
        mass=2.0,
        chi=10.0,
    )
    zero_state = zeros(Float64, 9)
    incoming_value = FermiHarmonics.nonlinear_ohmic_incoming_value(
        zero_state,
        SVector(1.0, 0.0),
        1.0,
        0.1,
        similar(zero_state),
        chi_eq,
        1.0e-12,
    )
    expected_incoming = FermiHarmonics.nonlinear_electrochemical_bias(0.1, chi_eq) /
                        (1.0 + chi_eq.electrostatic_coupling * (7.0 / 16.0))
    @test incoming_value ≈ expected_incoming atol=1e-12 rtol=1e-12
end

@testset "Solve smoke tests" begin
    mesh_path = normpath(joinpath(@__DIR__, "..", "projects", "square_bells_ucsb", "mesh", "square_bells.inp"))
    boundary_conditions = Dict(
        :walls => MaxwellWallBC(1.0),
        :contact_top => OhmicContactBC(-0.1),
        :contact_bottom => OhmicContactBC(0.1),
    )
    params = SolveParams(;
        polydeg=1,
        tspan_end=0.02,
        residual_tol=1e-3,
        cfl=0.4,
        log_every=10_000,
        min_harmonic=2,
        max_harmonic_auto=4,
    )

    sol_linear, semi_linear = solve(
        mesh_path,
        boundary_conditions,
        params,
        0.0,
        0.5;
        max_harmonic=2,
        name="test_linear",
    )
    @test length(sol_linear.u[end]) == length(Trixi.wrap_array(sol_linear.u[end], semi_linear))

    sol_nonlinear, semi_nonlinear = solve(
        mesh_path,
        boundary_conditions,
        params,
        0.0,
        0.5;
        transport=:parabolic_nonlinear,
        max_harmonic=2,
        mu0=1.0,
        mass=2.0,
        name="test_nonlinear",
    )
    @test semi_nonlinear.equations.transport === :parabolic_nonlinear
    @test semi_nonlinear.equations.collision_model === :quadratic_bgk
    @test semi_nonlinear.equations isa FermiHarmonics2D
    @test length(sol_nonlinear.u[end]) == length(Trixi.wrap_array(sol_nonlinear.u[end], semi_nonlinear))

    sol_exact, semi_exact = solve(
        mesh_path,
        boundary_conditions,
        params,
        0.0,
        0.5;
        transport=:parabolic_nonlinear,
        collision_model=:exact_bgk,
        n_angles=16,
        mu0=1.0,
        mass=2.0,
        name="test_exact_nonlinear",
    )
    @test semi_exact.equations.collision_model === :exact_bgk
    @test semi_exact.equations isa FermiAngles2D

    sol_two_rate, semi_two_rate = solve(
        mesh_path,
        boundary_conditions,
        params,
        0.0,
        0.5;
        transport=:parabolic_nonlinear,
        collision_model=:two_rate_bgk,
        n_angles=16,
        mu0=1.0,
        mass=2.0,
        name="test_two_rate_nonlinear",
    )
    @test semi_two_rate.equations.collision_model === :two_rate_bgk
    @test semi_two_rate.equations isa FermiAngles2D

    linear_probe = evaluate_observables(sol_linear, semi_linear, 0.0, 0.0)
    @test linear_probe.in_domain
    @test linear_probe.jx ≈ linear_probe.a1 atol=1e-10 rtol=1e-10
    @test linear_probe.jy ≈ linear_probe.b1 atol=1e-10 rtol=1e-10

    nonlinear_probe = evaluate_observables(sol_nonlinear, semi_nonlinear, 0.0, 0.0)
    @test nonlinear_probe.in_domain
    @test isfinite(nonlinear_probe.n)
    @test isfinite(nonlinear_probe.jx)
    @test isfinite(nonlinear_probe.jy)
    @test Trixi.varnames(FermiHarmonics.analysis_variables, semi_nonlinear.equations) == ("n", "jx", "jy")
    @test semi_nonlinear.equations.collision_model === :quadratic_bgk
    exact_probe = evaluate_observables(sol_exact, semi_exact, 0.0, 0.0)
    @test exact_probe.in_domain
    @test isfinite(exact_probe.n)
    two_rate_probe = evaluate_observables(sol_two_rate, semi_two_rate, 0.0, 0.0)
    @test two_rate_probe.in_domain
    @test isfinite(two_rate_probe.n)

    @test_throws ArgumentError solve(
        mesh_path,
        boundary_conditions,
        params,
        0.0,
        0.5;
        transport=:parabolic_nonlinear,
        collision_model=:exact_bgk,
        mu0=1.0,
        mass=2.0,
        name="missing_exact_angles",
    )
    @test_throws ArgumentError solve(
        mesh_path,
        boundary_conditions,
        params,
        0.0,
        0.5;
        transport=:parabolic_nonlinear,
        collision_model=:two_rate_bgk,
        mu0=1.0,
        mass=2.0,
        name="missing_two_rate_angles",
    )
    @test_throws ArgumentError solve(
        mesh_path,
        boundary_conditions,
        params,
        0.0,
        0.5;
        transport=:parabolic_nonlinear,
        n_angles=16,
        mu0=1.0,
        mass=2.0,
        name="invalid_quadratic_angles",
    )
    @test_throws ArgumentError solve(
        mesh_path,
        boundary_conditions,
        params,
        0.0,
        0.5;
        transport=:parabolic_nonlinear,
        collision_model=:exact_bgk,
        max_harmonic=2,
        n_angles=16,
        mu0=1.0,
        mass=2.0,
        name="invalid_exact_max_harmonic",
    )
    @test_throws ArgumentError solve(
        mesh_path,
        boundary_conditions,
        params,
        0.0,
        0.5;
        transport=:parabolic_nonlinear,
        collision_model=:two_rate_bgk,
        max_harmonic=2,
        n_angles=16,
        mu0=1.0,
        mass=2.0,
        name="invalid_two_rate_max_harmonic",
    )
    @test_throws ArgumentError solve(
        mesh_path,
        boundary_conditions,
        params,
        0.0,
        0.5;
        transport=:parabolic_nonlinear,
        collision_model=:linear_mrt,
        max_harmonic=2,
        mu0=1.0,
        mass=2.0,
        name="invalid_collision_model",
    )
    @test_throws ArgumentError solve(
        mesh_path,
        boundary_conditions,
        params,
        0.0,
        0.5;
        transport=:parabolic_nonlinear,
        collision_model=:exact_bgk,
        n_angles=16,
        mu0=1.0,
        mass=2.0,
        u0_override=zeros(length(sol_exact.u[end]) + 1),
        name="invalid_warm_start",
    )

    stability_params = SolveParams(;
        polydeg=1,
        tspan_end=0.08,
        residual_tol=1e-3,
        cfl=0.35,
        log_every=10_000,
        min_harmonic=2,
        max_harmonic_auto=4,
    )
    sol_stable, semi_stable = solve(
        mesh_path,
        boundary_conditions,
        stability_params,
        0.0,
        20.0;
        transport=:parabolic_nonlinear,
        max_harmonic=2,
        mu0=1.0,
        mass=2.0,
        name="test_nonlinear_stability",
    )
    @test semi_stable.equations.collision_model === :quadratic_bgk
    @test sol_stable.t[end] ≈ stability_params.tspan_end atol=1.0e-12 rtol=1.0e-12
    status_stable = FermiHarmonics.solve_status(sol_stable, semi_stable, stability_params)
    @test status_stable.hit_final_time
    @test status_stable.stop_reason === :final_time
    @test !status_stable.converged

    sol_gamma3, semi_gamma3 = solve(
        mesh_path,
        boundary_conditions,
        params,
        0.0,
        0.5;
        transport=:parabolic_nonlinear,
        max_harmonic=2,
        gamma3=0.2,
        mu0=1.0,
        mass=2.0,
        name="test_nonlinear_gamma3",
    )
    @test semi_gamma3.equations.gamma3 ≈ 0.2 atol=1e-12 rtol=1e-12

    sol_chi, semi_chi = solve(
        mesh_path,
        boundary_conditions,
        params,
        0.0,
        0.5;
        transport=:parabolic_nonlinear,
        max_harmonic=2,
        mu0=1.0,
        mass=2.0,
        chi=0.3,
        name="test_nonlinear_chi_hyperbolic",
    )
    @test semi_chi isa Trixi.SemidiscretizationHyperbolic
    @test semi_chi.equations.electrostatic_coupling ≈ 0.3 atol=1e-12 rtol=1e-12
    @test sol_chi.t[end] > 0.0

    @test_throws ArgumentError solve(
        mesh_path,
        boundary_conditions,
        params,
        0.0,
        0.5;
        transport=:parabolic_nonlinear,
        max_harmonic=2,
        gamma3=-0.1,
        mu0=1.0,
        mass=2.0,
        name="invalid_gamma3",
    )
end

@testset "Tesla cluster sweep helpers" begin
    ENV["TESLA_SWEEP_PARTITION"] = "test-partition"
    ENV["TESLA_SWEEP_QOS"] = "test-qos"
    ENV["TESLA_SWEEP_TIME"] = "02:00:00"
    include(normpath(joinpath(@__DIR__, "..", "projects", "nonlinear", "scripts", "run_tesla_valve_cluster_bias_sweep.jl")))
    include(normpath(joinpath(@__DIR__, "..", "projects", "nonlinear", "scripts", "submit_tesla_valve_cluster_bias_sweep.jl")))

    biases = tesla_cluster_bias_values()
    @test length(biases) == 30
    @test first(biases) ≈ 0.01 atol=1e-12 rtol=1e-12
    @test last(biases) ≈ 0.5 atol=1e-12 rtol=1e-12
    @test all(diff(biases) .> 0.0)

    forward_bcs = tesla_cluster_boundary_conditions("forward", 0.2, 1.0)
    reverse_bcs = tesla_cluster_boundary_conditions("reverse", 0.2, 1.0)
    @test forward_bcs[:inlet].bias ≈ 0.1 atol=1e-12 rtol=1e-12
    @test forward_bcs[:outlet].bias ≈ -0.1 atol=1e-12 rtol=1e-12
    @test reverse_bcs[:inlet].bias ≈ -0.1 atol=1e-12 rtol=1e-12
    @test reverse_bcs[:outlet].bias ≈ 0.1 atol=1e-12 rtol=1e-12

    commands = tesla_cluster_submission_commands(; dry_run=true)
    @test Set(keys(commands)) == Set(("forward", "reverse"))
    forward_cmd = sprint(show, commands["forward"])
    reverse_cmd = sprint(show, commands["reverse"])
    @test occursin("JULIA_NUM_THREADS=16", forward_cmd)
    @test occursin("TESLA_DIRECTION=forward", forward_cmd)
    @test occursin("--cpus-per-task=16", forward_cmd)
    @test occursin("--mem=8G", forward_cmd)
    @test occursin("--partition=test-partition", forward_cmd)
    @test occursin("--qos=test-qos", forward_cmd)
    @test occursin("--time=02:00:00", forward_cmd)
    @test occursin("TESLA_DIRECTION=reverse", reverse_cmd)
    @test occursin("--array=1-1", reverse_cmd)
end
