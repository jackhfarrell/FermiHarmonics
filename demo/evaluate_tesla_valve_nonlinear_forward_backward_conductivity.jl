using HDF5
using TOML

function trapz(x_values, y_values)
    length(x_values) == length(y_values) || throw(ArgumentError("trapz inputs must have the same length"))
    length(x_values) >= 2 || return 0.0

    total = 0.0
    @inbounds for i in 1:(length(x_values) - 1)
        total += 0.5 * (x_values[i + 1] - x_values[i]) * (y_values[i + 1] + y_values[i])
    end
    return total
end

function load_cartesian_analysis(path::AbstractString)
    isfile(path) || throw(ArgumentError("analysis file not found: $(path)"))

    return h5open(path, "r") do file
        haskey(file, "x") || throw(ArgumentError("missing dataset 'x' in $(path)"))
        haskey(file, "y") || throw(ArgumentError("missing dataset 'y' in $(path)"))
        haskey(file, "jx") || throw(ArgumentError("missing dataset 'jx' in $(path)"))
        haskey(file, "mask") || throw(ArgumentError("missing dataset 'mask' in $(path)"))

        x = read(file["x"])
        y = read(file["y"])
        jx = read(file["jx"])
        mask_raw = read(file["mask"])
        mask = mask_raw .!= 0

        return (
            x = x,
            y = y,
            jx = jx,
            mask = mask,
        )
    end
end

function integrated_cross_section_current(grids, target_x)
    x_index = argmin(abs.(grids.x .- target_x))
    line_mask = vec(grids.mask[x_index, :])
    y_line = collect(grids.y[line_mask])
    jx_line = collect(grids.jx[x_index, line_mask])

    length(y_line) >= 2 || throw(ArgumentError("cross-section has fewer than 2 valid points at x=$(grids.x[x_index])"))
    order = sortperm(y_line)
    integrated = trapz(y_line[order], jx_line[order])

    return (
        current = integrated,
        x = Float64(grids.x[x_index]),
        num_points = length(y_line),
    )
end

function write_summary_csv(path::AbstractString, rows)
    open(path, "w") do io
        println(io, "direction,bias_magnitude,integrated_jx,conductivity")
        for row in rows
            println(io, join((row.direction, row.bias_magnitude, row.integrated_jx, row.conductivity), ","))
        end
    end
    return path
end

function main()
    default_forward = joinpath(
        @__DIR__,
        "data_nonlinear",
        "tesla_valve_live_mesh_native_gamma_mc02_M20_bias02_chi000_poly3_t2.h5",
    )
    default_backward = joinpath(
        @__DIR__,
        "data_nonlinear",
        "tesla_valve_live_mesh_native_gamma_mc02_M20_bias02_reversed_chi000_poly3_t2.h5",
    )

    forward_path = length(ARGS) >= 1 ? ARGS[1] : get(ENV, "TESLA_FORWARD_FILE", default_forward)
    backward_path = length(ARGS) >= 2 ? ARGS[2] : get(ENV, "TESLA_BACKWARD_FILE", default_backward)
    bias_magnitude = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : parse(Float64, get(ENV, "TESLA_BIAS_MAGNITUDE", "0.2"))
    target_x = length(ARGS) >= 4 ? parse(Float64, ARGS[4]) : parse(Float64, get(ENV, "TESLA_TARGET_X", "0.0"))
    output_dir = length(ARGS) >= 5 ? ARGS[5] : get(ENV, "TESLA_COMPARE_OUTPUT_DIR", dirname(forward_path))
    output_prefix = length(ARGS) >= 6 ? ARGS[6] : get(ENV, "TESLA_COMPARE_PREFIX", "tesla_valve_existing_forward_backward")

    mkpath(output_dir)

    forward_grids = load_cartesian_analysis(forward_path)
    backward_grids = load_cartesian_analysis(backward_path)

    forward_section = integrated_cross_section_current(forward_grids, target_x)
    backward_section = integrated_cross_section_current(backward_grids, target_x)

    forward_conductivity = abs(forward_section.current) / bias_magnitude
    backward_conductivity = abs(backward_section.current) / bias_magnitude
    rectification_ratio = forward_conductivity / backward_conductivity
    relative_contrast = (forward_conductivity - backward_conductivity) /
        (0.5 * (forward_conductivity + backward_conductivity))

    rows = [
        (
            direction="forward",
            bias_magnitude=bias_magnitude,
            integrated_jx=forward_section.current,
            conductivity=forward_conductivity,
        ),
        (
            direction="backward",
            bias_magnitude=bias_magnitude,
            integrated_jx=backward_section.current,
            conductivity=backward_conductivity,
        ),
    ]

    csv_path = joinpath(output_dir, output_prefix * "_conductivity_comparison.csv")
    write_summary_csv(csv_path, rows)

    summary_path = joinpath(output_dir, output_prefix * "_summary.txt")
    open(summary_path, "w") do io
        println(io, "Tesla valve nonlinear forward/backward conductivity comparison")
        println(io, "forward_file=$(abspath(forward_path))")
        println(io, "backward_file=$(abspath(backward_path))")
        println(io, "bias_magnitude=$(bias_magnitude)")
        println(io, "target_x_requested=$(target_x)")
        println(io, "target_x_forward_actual=$(forward_section.x)")
        println(io, "target_x_backward_actual=$(backward_section.x)")
        println(io, "forward_integrated_jx=$(forward_section.current)")
        println(io, "backward_integrated_jx=$(backward_section.current)")
        println(io, "forward_conductivity=$(forward_conductivity)")
        println(io, "backward_conductivity=$(backward_conductivity)")
        println(io, "rectification_ratio_forward_over_backward=$(rectification_ratio)")
        println(io, "relative_conductivity_contrast=$(relative_contrast)")
    end

    metadata_path = joinpath(output_dir, output_prefix * "_metadata.toml")
    metadata = Dict(
        "forward_file" => abspath(forward_path),
        "backward_file" => abspath(backward_path),
        "bias_magnitude" => bias_magnitude,
        "target_x_requested" => target_x,
        "target_x_forward_actual" => forward_section.x,
        "target_x_backward_actual" => backward_section.x,
        "forward_integrated_jx" => forward_section.current,
        "backward_integrated_jx" => backward_section.current,
        "forward_conductivity" => forward_conductivity,
        "backward_conductivity" => backward_conductivity,
        "rectification_ratio_forward_over_backward" => rectification_ratio,
        "relative_conductivity_contrast" => relative_contrast,
    )
    open(metadata_path, "w") do io
        TOML.print(io, metadata)
    end

    @info "Forward/backward conductivity comparison complete (using existing files)" csv_path summary_path metadata_path rectification_ratio relative_contrast
    return nothing
end

main()