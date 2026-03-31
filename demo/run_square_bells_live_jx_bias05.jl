using GLMakie
using ElectronKinetics, Trixi

function trapz(x_values, y_values)
    length(x_values) == length(y_values) || throw(ArgumentError("trapz inputs must have the same length"))
    length(x_values) >= 2 || return 0.0
    total = 0.0
    @inbounds for i in 1:(length(x_values) - 1)
        total += 0.5 * (x_values[i + 1] - x_values[i]) * (y_values[i + 1] + y_values[i])
    end
    return total
end

function integrated_horizontal_current(solution_vector, semi, target_y; nvisnodes::Int=301)
    grids = ElectronKinetics.compute_analysis_grids(solution_vector, semi; nvisnodes=nvisnodes)
    y_index = argmin(abs.(grids.y .- target_y))
    line_mask = vec(grids.mask[:, y_index])
    x_line = collect(grids.x[line_mask])
    jy_line = collect(grids.jy[line_mask, y_index])
    length(x_line) >= 2 || return (current=0.0, average=0.0, y=grids.y[y_index], num_points=length(x_line))
    order = sortperm(x_line)
    x_sorted = x_line[order]
    jy_sorted = jy_line[order]
    integrated = trapz(x_sorted, jy_sorted)
    span = maximum(x_sorted) - minimum(x_sorted)
    average = iszero(span) ? jy_sorted[1] : integrated / span
    return (current=integrated, average=average, y=grids.y[y_index], num_points=length(x_sorted))
end

function main()
    project_root = normpath(joinpath(@__DIR__, ".."))
    mesh_path = joinpath(project_root, "projects", "square_bells_ucsb", "mesh", "square_bells.inp")

    reference = ElectronKinetics.blg_reference_setup()
    mu0       = reference.mu0
    mass      = reference.mass
    gamma_mr  = reference.gamma_mr
    gamma_ee  = 100.0
    target_current = 1.0
    chi       = 10.0
    p_scatter = reference.p_scatter

    boundary_conditions = Dict(
        :walls          => MaxwellWallBC(p_scatter),
        :contact_top    => CurrentContactBC( target_current / 2),
        :contact_bottom => CurrentContactBC(-target_current / 2),
    )

    params = SolveParams(;
        polydeg           = 3,
        tspan_end         = 1e8,
        residual_tol      = 1e-3,
        cfl               = 0.8,
        log_every         = 500,
        min_harmonic      = 4,
        max_harmonic_auto = 20,
    )

    sol, semi = ElectronKinetics.solve(
        mesh_path,
        boundary_conditions,
        params,
        gamma_mr,
        gamma_ee;
        transport          = :parabolic_nonlinear,
        max_harmonic       = :auto,
        mu0                = mu0,
        mass               = mass,
        chi                = chi,
        visualize     = true,
        viz_field     = :jx,
        viz_colormap  = :RdBu,
        name               = "square_bells_jx_current_drive10",
    )

    measured = integrated_horizontal_current(sol.u[end], semi, 0.0)
    @info "Current-drive summary" target_current measured_midline_avg_jy=measured.average measured_midline_integrated_jy=measured.current measured_line_y=measured.y measured_points=measured.num_points

    return nothing
end

main()
