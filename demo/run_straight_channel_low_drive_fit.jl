using ElectronKinetics, Trixi
using Plots
using TOML

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
    grids = ElectronKinetics.compute_analysis_grids(solution_vector, semi; nvisnodes=nvisnodes)
    current_grid = something(grids.jx, grids.a1)
    x_index = argmin(abs.(grids.x .- target_x))
    line_mask = vec(grids.mask[x_index, :])
    y_line = collect(grids.y[line_mask])
    jx_line = collect(current_grid[x_index, line_mask])
    return trapz(y_line, jx_line)
end

function fit_small_drive(rows)
    x = [row.bias^2 for row in rows]
    y = [row.integrated_jx_over_bias for row in rows]
    design = hcat(ones(length(x)), x)
    coeffs = design \ y
    y_fit = design * coeffs
    residuals = y .- y_fit
    rms = sqrt(sum(residuals .^ 2) / length(residuals))
    return (c1 = coeffs[1], c3 = coeffs[2], rms = rms)
end

function write_rows(path::AbstractString, rows)
    open(path, "w") do io
        println(io, "bias,bias_squared,integrated_jx,integrated_jx_over_bias")
        for row in rows
            println(io, join((row.bias, row.bias^2, row.integrated_jx, row.integrated_jx_over_bias), ","))
        end
    end
    return path
end

function write_reference_metadata(path::AbstractString, reference, geo_path::AbstractString, mesh_path::AbstractString)
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
        "mesh_geo" => basename(geo_path),
        "mesh_file" => basename(mesh_path),
        "analysis_case" => "reference_low_drive",
    )
    open(path, "w") do io
        TOML.print(io, metadata)
    end
    return path
end

function write_fit_summary(path::AbstractString, fit, fit_max_bias::Float64)
    open(path, "w") do io
        println(io, "Model: integrated_jx_over_bias = c1 + c3 * bias^2")
        println(io, "fit_bias_max=$(fit_max_bias)")
        println(io, "c1=$(fit.c1)")
        println(io, "c3=$(fit.c3)")
        println(io, "rms=$(fit.rms)")
    end
    return path
end

function main()
    project_root = normpath(joinpath(@__DIR__, ".."))
    geo_path = joinpath(project_root, "demo", "mesh", "straight_channel.geo")
    mesh_path = joinpath(project_root, "demo", "mesh", "straight_channel.inp")
    output_dir = joinpath(project_root, "demo", "data_straight_channel_low_drive_fit")
    mkpath(output_dir)
    ensure_straight_channel_mesh(geo_path, mesh_path)

    reference = ElectronKinetics.blg_reference_setup()
    mu0 = reference.mu0
    mass = reference.mass
    p_scatter = reference.p_scatter
    gamma_mr = reference.gamma_mr
    gamma_mc = reference.gamma_mc
    probe_x = reference.probe_x

    bias_values = [
        1.0e-3,
        2.0e-3,
        5.0e-3,
        1.0e-2,
        2.0e-2,
        5.0e-2,
        1.0e-1,
    ]
    fit_max_bias = 5.0e-2

    params = SolveParams(;
        polydeg = 2,
        tspan_end = 80.0,
        residual_tol = 1e-7,
        cfl = 0.5,
        log_every = 500,
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

        run_name = "reference_straight_channel_low_drive_" * replace(string(bias), "." => "p")
        @info "Running low-drive fit case" bias gamma_mr gamma_mc residual_tol=params.residual_tol

        sol, semi = ElectronKinetics.solve(
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

        integrated_jx = integrated_cross_section_current(sol.u[end], semi, probe_x)
        push!(rows, (
            bias = bias,
            integrated_jx = integrated_jx,
            integrated_jx_over_bias = integrated_jx / bias,
        ))

        u0 = copy(sol.u[end])
    end

    fit_rows = [row for row in rows if row.bias <= fit_max_bias]
    fit = fit_small_drive(fit_rows)

    csv_path = joinpath(output_dir, "low_drive_fit_data.csv")
    metadata_path = joinpath(output_dir, "low_drive_fit_metadata.toml")
    write_rows(csv_path, rows)
    write_reference_metadata(metadata_path, reference, geo_path, mesh_path)

    fig = plot(
        [row.bias^2 for row in fit_rows],
        [row.integrated_jx_over_bias for row in fit_rows];
        marker = :circle,
        linewidth = 0,
        markersize = 6,
        xlabel = "bias^2",
        ylabel = "I(x*) / bias",
        label = "fit-window data",
        title = "Reference straight-channel low-drive cubic-fit check",
        legend = :topright,
    )
    fit_x = range(0.0, maximum([row.bias^2 for row in fit_rows]) * 1.05, length=200)
    plot!(
        fig,
        fit_x,
        fit.c1 .+ fit.c3 .* fit_x;
        linewidth = 2,
        color = :red,
        label = "c1 + c3 bias^2",
    )

    png_path = joinpath(output_dir, "low_drive_cubic_fit.png")
    savefig(fig, png_path)

    fit_summary_path = joinpath(output_dir, "low_drive_fit_summary.txt")
    write_fit_summary(fit_summary_path, fit, fit_max_bias)

    @info "Low-drive fit complete" csv_path metadata_path png_path fit_summary_path c1=fit.c1 c3=fit.c3 rms=fit.rms
    return nothing
end

main()
