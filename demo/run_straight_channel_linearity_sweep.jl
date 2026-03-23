using FermiHarmonics
using Plots
using TOML

function get_env_float(name::AbstractString, default::Float64)
    value = get(ENV, name, nothing)
    return isnothing(value) ? default : parse(Float64, value)
end

function gamma_label(gamma_value::Real)
    return replace(string(round(Float64(gamma_value); sigdigits=6)), "." => "p")
end

function mesh_has_nodesets(mesh_path::AbstractString, names)
    isfile(mesh_path) || return false
    contents = read(mesh_path, String)
    return all(occursin("*NSET,NSET=$(name)", contents) for name in names)
end

function ensure_straight_channel_mesh(geo_path::AbstractString, mesh_path::AbstractString)
    required_nodesets = ("inlet", "outlet", "walls")
    needs_rebuild = !isfile(mesh_path) ||
                    mtime(mesh_path) < mtime(geo_path) ||
                    !mesh_has_nodesets(mesh_path, required_nodesets)
    needs_rebuild || return mesh_path

    mkpath(dirname(mesh_path))
    cmd = Cmd([
        "gmsh",
        "-2",
        geo_path,
        "-string",
        "Mesh.SaveGroupsOfNodes=1;",
        "-format",
        "inp",
        "-o",
        mesh_path,
    ])
    @info "Building straight-channel mesh" geo_path mesh_path
    run(cmd)
    return mesh_path
end

function write_probe_table(path::AbstractString, rows)
    open(path, "w") do io
        println(
            io,
            "bias,jx,jx_over_bias,a1,a1_over_bias,a0,integrated_jx,integrated_jx_over_bias," *
            "left_avg_a0,right_avg_a0,measured_avg_drop,integrated_jx_over_measured_avg_drop",
        )
        for row in rows
            println(io, join((
                row.bias,
                row.jx,
                row.jx_over_bias,
                row.a1,
                row.a1_over_bias,
                row.a0,
                row.integrated_jx,
                row.integrated_jx_over_bias,
                row.left_avg_a0,
                row.right_avg_a0,
                row.measured_avg_drop,
                row.integrated_jx_over_measured_avg_drop,
            ), ","))
        end
    end
    return path
end

function write_reference_metadata(path::AbstractString, reference, geo_path::AbstractString, mesh_path::AbstractString;
                                  analysis_case::AbstractString, gamma_mc::Float64)
    metadata = Dict(
        "convention_name" => reference.convention_name,
        "convention_version" => reference.convention_version,
        "channel_length" => reference.channel_length,
        "half_height" => reference.half_height,
        "mu0" => reference.mu0,
        "mass" => reference.mass,
        "vF" => reference.vF,
        "gamma_mr" => reference.gamma_mr,
        "gamma_mc" => reference.gamma_mc,
        "p_scatter" => reference.p_scatter,
        "transport" => String(reference.transport),
        "probe_x" => reference.probe_x,
        "probe_y" => reference.probe_y,
        "left_probe_x" => reference.left_probe_x,
        "right_probe_x" => reference.right_probe_x,
        "mesh_geo" => basename(geo_path),
        "mesh_file" => basename(mesh_path),
        "analysis_case" => analysis_case,
        "sweep_gamma_mc" => gamma_mc,
    )
    open(path, "w") do io
        TOML.print(io, metadata)
    end
    return path
end

function trapz(x_values, y_values)
    length(x_values) == length(y_values) || throw(ArgumentError("trapz inputs must have the same length"))
    length(x_values) >= 2 || return 0.0

    total = 0.0
    @inbounds for i in 1:(length(x_values) - 1)
        total += 0.5 * (x_values[i + 1] - x_values[i]) * (y_values[i + 1] + y_values[i])
    end
    return total
end

function integrated_cross_section_current(solution_vector, semi, target_x; nvisnodes::Int=301)
    grids = FermiHarmonics.compute_analysis_grids(solution_vector, semi; nvisnodes=nvisnodes)
    current_grid = something(grids.jx, grids.a1)
    x_index = argmin(abs.(grids.x .- target_x))
    line_mask = vec(grids.mask[x_index, :])
    y_line = collect(grids.y[line_mask])
    jx_line = collect(current_grid[x_index, line_mask])
    integrated = trapz(y_line, jx_line)
    return (
        current = integrated,
        x = grids.x[x_index],
        num_points = length(y_line),
    )
end

function cross_section_average_a0(solution_vector, semi, target_x; nvisnodes::Int=301)
    grids = FermiHarmonics.compute_analysis_grids(solution_vector, semi; nvisnodes=nvisnodes)
    x_index = argmin(abs.(grids.x .- target_x))
    line_mask = vec(grids.mask[x_index, :])
    y_line = collect(grids.y[line_mask])
    a0_line = collect(grids.a0[x_index, line_mask])
    span = maximum(y_line) - minimum(y_line)
    average = iszero(span) ? a0_line[1] : trapz(y_line, a0_line) / span
    return (
        average = average,
        x = grids.x[x_index],
        num_points = length(y_line),
    )
end

function fit_small_drive_breakdown(rows; max_fit_bias::Float64=0.2154434690031884)
    fit_rows = [row for row in rows if row.bias <= max_fit_bias]
    length(fit_rows) >= 2 || throw(ArgumentError("Need at least two low-drive points for the fit"))

    x = [row.bias^2 for row in fit_rows]
    y = [row.integrated_jx_over_bias for row in fit_rows]
    design = hcat(ones(length(x)), x)
    coeffs = design \ y
    c1 = coeffs[1]
    c3 = coeffs[2]

    return (
        c1 = c1,
        c3 = c3,
        max_fit_bias = max_fit_bias,
        fit_rows = fit_rows,
    )
end

function write_fit_summary(path::AbstractString, fit)
    open(path, "w") do io
        println(io, "Model: integrated_jx_over_bias = c1 + c3 * bias^2")
        println(io, "fit_bias_max=$(fit.max_fit_bias)")
        println(io, "c1=$(fit.c1)")
        println(io, "c3=$(fit.c3)")
    end
    return path
end

function main()
    project_root = normpath(joinpath(@__DIR__, ".."))
    geo_path = joinpath(project_root, "demo", "mesh", "straight_channel.geo")
    mesh_path = joinpath(project_root, "demo", "mesh", "straight_channel.inp")
    reference = FermiHarmonics.blg_reference_setup()
    mu0 = reference.mu0
    mass = reference.mass
    p_scatter = reference.p_scatter
    gamma_mr = reference.gamma_mr
    gamma_mc = get_env_float("STRAIGHT_CHANNEL_GAMMA_MC", reference.gamma_mc)
    probe_x = reference.probe_x
    probe_y = reference.probe_y
    left_probe_x = reference.left_probe_x
    right_probe_x = reference.right_probe_x
    analysis_case = iszero(gamma_mc - reference.gamma_mc) ? "reference" : "comparison"
    output_dir_default = analysis_case == "reference" ?
        joinpath(project_root, "demo", "data_straight_channel_linearity") :
        joinpath(project_root, "demo", "data_straight_channel_linearity_gamma_mc_" * gamma_label(gamma_mc))
    output_dir = get(ENV, "STRAIGHT_CHANNEL_OUTPUT_DIR", output_dir_default)
    mkpath(output_dir)
    ensure_straight_channel_mesh(geo_path, mesh_path)

    bias_values = collect(10.0 .^ range(log10(0.01), log10(1.0), length=7))

    params = SolveParams(;
        polydeg = 2,
        tspan_end = 40.0,
        residual_tol = 1e-4,
        cfl = 0.5,
        log_every = 200,
        min_harmonic = 4,
        max_harmonic_auto = 20,
    )

    rows = NamedTuple[]
    u0 = nothing

    for bias in bias_values
        boundary_conditions = Dict(
            :walls => MaxwellWallBC(p_scatter),
            :inlet => OhmicContactBC(bias / 2),
            :outlet => OhmicContactBC(-bias / 2),
        )

        run_prefix = analysis_case == "reference" ? "reference_straight_channel_bias_" :
            "comparison_gamma_mc_" * gamma_label(gamma_mc) * "_bias_"
        run_name = run_prefix * replace(string(round(bias; sigdigits=4)), "." => "p")
        @info "Running straight-channel sweep case" bias gamma_mr gamma_mc mu0 mass

        sol, semi = FermiHarmonics.solve(
            mesh_path,
            boundary_conditions,
            params,
            gamma_mr,
            gamma_mc;
            transport = :parabolic_nonlinear,
            max_harmonic = :auto,
            mu0 = mu0,
            mass = mass,
            u0_override = u0,
            visualize = false,
            name = run_name,
        )

        probe = FermiHarmonics.evaluate_observables(sol, semi, probe_x, probe_y)
        probe.in_domain || error("Probe point ($(probe_x), $(probe_y)) is outside the straight-channel mesh")
        section = integrated_cross_section_current(sol.u[end], semi, probe_x)
        left_avg = cross_section_average_a0(sol.u[end], semi, left_probe_x)
        right_avg = cross_section_average_a0(sol.u[end], semi, right_probe_x)
        measured_avg_drop = left_avg.average - right_avg.average

        push!(rows, (
            bias = bias,
            jx = probe.jx,
            jx_over_bias = probe.jx / bias,
            a1 = probe.a1,
            a1_over_bias = probe.a1 / bias,
            a0 = probe.a0,
            integrated_jx = section.current,
            integrated_jx_over_bias = section.current / bias,
            left_avg_a0 = left_avg.average,
            right_avg_a0 = right_avg.average,
            measured_avg_drop = measured_avg_drop,
            integrated_jx_over_measured_avg_drop = section.current / measured_avg_drop,
        ))

        u0 = copy(sol.u[end])
    end

    csv_path = joinpath(output_dir, "probe_sweep.csv")
    metadata_path = joinpath(output_dir, "probe_sweep_metadata.toml")
    write_probe_table(csv_path, rows)
    write_reference_metadata(metadata_path, reference, geo_path, mesh_path;
                             analysis_case=analysis_case, gamma_mc=gamma_mc)

    baseline = rows[1].jx_over_bias
    fig = plot(
        [row.bias for row in rows],
        [row.jx_over_bias for row in rows];
        xscale = :log10,
        marker = :circle,
        linewidth = 2,
        xlabel = "drive amplitude (bias)",
        ylabel = "jx(x*, y*) / bias",
        label = "probe at (0, 0)",
        title = analysis_case == "reference" ?
            "Reference straight-channel linearity breakdown" :
            "Comparison straight-channel linearity breakdown (gamma_mc = $(round(gamma_mc; sigdigits=4)))",
        legend = :topright,
    )
    hline!(
        fig,
        [baseline];
        linestyle = :dash,
        linewidth = 2,
        color = :black,
        label = "linear-response baseline",
    )

    png_path = joinpath(output_dir, "straight_channel_linearity_breakdown.png")
    savefig(fig, png_path)

    integrated_baseline = rows[1].integrated_jx_over_bias
    fit = fit_small_drive_breakdown(rows)
    integrated_fig = plot(
        [row.bias for row in rows],
        [row.integrated_jx_over_bias for row in rows];
        xscale = :log10,
        marker = :circle,
        linewidth = 2,
        xlabel = "drive amplitude (bias)",
        ylabel = "I(x*) / bias",
        label = "cross-section current at x = 0",
        title = analysis_case == "reference" ?
            "Reference straight-channel integrated current breakdown" :
            "Comparison straight-channel integrated current breakdown (gamma_mc = $(round(gamma_mc; sigdigits=4)))",
        legend = :topright,
    )
    hline!(
        integrated_fig,
        [integrated_baseline];
        linestyle = :dash,
        linewidth = 2,
        color = :black,
        label = "linear-response baseline",
    )
    annotate!(
        integrated_fig,
        (0.45, minimum([row.integrated_jx_over_bias for row in rows]) + 0.0006,
         text("low-drive fit: c3 = $(round(fit.c3; sigdigits=4))", 10, :red)),
    )

    integrated_png_path = joinpath(output_dir, "straight_channel_integrated_current_breakdown.png")
    savefig(integrated_fig, integrated_png_path)

    measured_baseline = rows[1].integrated_jx_over_measured_avg_drop
    measured_fig = plot(
        [row.measured_avg_drop for row in rows],
        [row.integrated_jx_over_measured_avg_drop for row in rows];
        xscale = :log10,
        marker = :circle,
        linewidth = 2,
        xlabel = "measured average drop (\$\\Delta a_0\$)",
        ylabel = "I(x*) / measured average drop",
        label = "cross-section current at x = 0",
        title = analysis_case == "reference" ?
            "Reference straight-channel current vs measured average drop" :
            "Comparison straight-channel current vs measured average drop (gamma_mc = $(round(gamma_mc; sigdigits=4)))",
        legend = :topright,
    )
    hline!(
        measured_fig,
        [measured_baseline];
        linestyle = :dash,
        linewidth = 2,
        color = :black,
        label = "low-drive baseline",
    )

    measured_png_path = joinpath(output_dir, "straight_channel_measured_drop_breakdown.png")
    savefig(measured_fig, measured_png_path)

    fit_x = [row.bias^2 for row in fit.fit_rows]
    fit_y = [row.integrated_jx_over_bias for row in fit.fit_rows]
    fit_line_x = range(0.0, maximum(fit_x) * 1.05, length=120)
    fit_line_y = fit.c1 .+ fit.c3 .* fit_line_x
    cubic_fit_fig = plot(
        fit_x,
        fit_y;
        marker = :circle,
        linewidth = 0,
        markersize = 5,
        xlabel = "bias^2",
        ylabel = "I(x*) / bias",
        label = "low-drive data",
        title = analysis_case == "reference" ?
            "Reference straight-channel low-drive cubic-fit check" :
            "Comparison straight-channel low-drive cubic-fit check (gamma_mc = $(round(gamma_mc; sigdigits=4)))",
        legend = :topright,
    )
    plot!(
        cubic_fit_fig,
        fit_line_x,
        fit_line_y;
        linewidth = 2,
        color = :red,
        label = "c1 + c3 bias^2",
    )
    cubic_fit_png_path = joinpath(output_dir, "straight_channel_integrated_current_cubic_fit.png")
    savefig(cubic_fit_fig, cubic_fit_png_path)

    fit_summary_path = joinpath(output_dir, "straight_channel_integrated_current_fit.txt")
    write_fit_summary(fit_summary_path, fit)

    @info "Straight-channel linearity sweep complete" csv_path metadata_path png_path integrated_png_path cubic_fit_png_path fit_summary_path c1=fit.c1 c3=fit.c3
    return nothing
end

main()
