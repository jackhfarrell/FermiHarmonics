# Run two rectangle devices (full height and half height) at identical parameters,
# save outputs, and estimate contact resistance vs device length.
#
# Usage:
#   JULIA_NUM_THREADS=4 julia --project=. projects/simple_geometries/scripts/run_rectangle_contact_resistance.jl

using Dates
using DrWatson
using FermiHarmonics
using HDF5
using Plots
using StaticArrays
using Statistics
using Trixi


# ======================================================================================================================
# Configuration
# ======================================================================================================================

project_root = normpath(joinpath(@__DIR__, ".."))
results_root = joinpath(project_root, "data", "contact_resistance_rectangles")

visualize = true
nvisnodes = 400

bias = 1.0
p_scatter = 1.0
gamma_mr = 0.0
gamma_mc = 1e3

params = SolveParams(;
    polydeg = 3,
    tspan_end = 100.0,
    residual_tol = 1e-5,
    cfl = 0.5,
    log_every = 1000,
)

max_harmonic = 3

devices = [
    (
        name = "rectangle",
        mesh_path = joinpath(project_root, "meshes", "rectangle", "rectangle.inp"),
    ),
    (
        name = "rectangle_halfheight",
        mesh_path = joinpath(project_root, "meshes", "rectangle_half", "rectangle_halfheight.inp"),
    ),
]

function make_boundary_conditions(bias, p_scatter)
    return Dict(
        :walls => MaxwellWallBC(p_scatter),
        :contact_top => OhmicContactBC(-bias / 2),
        :contact_bottom => OhmicContactBC(bias / 2),
    )
end


# ======================================================================================================================
# Live Visualization Override (safe against NaN/Inf in plotting)
# ======================================================================================================================

@inline function safe_current_norm_variables(u, equations::FermiHarmonics.FermiHarmonics2D{NVARS}) where {NVARS}
    a1 = length(u) >= 2 && isfinite(u[2]) ? u[2] : 0.0
    b1 = length(u) >= 3 && isfinite(u[3]) ? u[3] : 0.0
    j_norm = hypot(a1, b1)
    return SVector{NVARS, Float64}(ntuple(i -> i == 1 ? j_norm : 0.0, NVARS))
end

function Trixi.varnames(::typeof(safe_current_norm_variables), equations::FermiHarmonics.FermiHarmonics2D{NVARS}) where {NVARS}
    return ntuple(i -> i == 1 ? "j_norm" : "_viz_pad_$(i)", NVARS)
end

function FermiHarmonics.visualization_callback(params::FermiHarmonics.SolveParams, semi, name::AbstractString)
    return Trixi.VisualizationCallback(
        semi;
        interval=params.log_every,
        solution_variables=safe_current_norm_variables,
        variable_names=["j_norm"],
        filename="live_viz_$(name)_j_norm",
        overwrite=true,
        seriescolor=:magma,
    )
end


# ======================================================================================================================
# Analysis Helpers (simple contact resistance proxy)
# ======================================================================================================================

function mean_profile_over_x(field, mask)
    _, ny = size(field)
    profile = fill(NaN, ny)
    @inbounds for j in 1:ny
        valid = mask[:, j] .& .!isnan.(field[:, j])
        if any(valid)
            profile[j] = mean(field[findall(valid), j])
        end
    end
    return profile
end

function band_mean(y, values, y_lo, y_hi)
    in_band = (y .>= y_lo) .& (y .<= y_hi) .& .!isnan.(values)
    return any(in_band) ? mean(values[in_band]) : NaN
end

function analyze_contact_resistance(h5_path, bias)
    x = nothing
    y = nothing
    a0 = nothing
    b1 = nothing
    mask = nothing
    h5open(h5_path, "r") do f
        x = read(f["x"])
        y = read(f["y"])
        a0 = read(f["a0"])
        b1 = read(f["b1"])
        mask = read(f["mask"])
    end

    a0_y = mean_profile_over_x(a0, mask)
    j_y = mean_profile_over_x(abs.(b1), mask)

    y_min = minimum(y)
    y_max = maximum(y)
    length = y_max - y_min
    width = maximum(x) - minimum(x)

    # Bands are fractions of device length.
    # Use near-contact bands for contact jumps and interior bands for bulk resistivity.
    bottom_contact_bulk = band_mean(y, a0_y, y_min + 0.10 * length, y_min + 0.20 * length)
    top_contact_bulk = band_mean(y, a0_y, y_max - 0.20 * length, y_max - 0.10 * length)
    bottom_inner_bulk = band_mean(y, a0_y, y_min + 0.30 * length, y_min + 0.40 * length)
    top_inner_bulk = band_mean(y, a0_y, y_max - 0.40 * length, y_max - 0.30 * length)
    current_bulk = band_mean(y, j_y, y_min + 0.45 * length, y_min + 0.55 * length)

    v_bottom = bias / 2
    v_top = -bias / 2
    jump_bottom = abs(v_bottom - bottom_contact_bulk)
    jump_top = abs(v_top - top_contact_bulk)

    if isnan(current_bulk) || current_bulk <= 0 || isnan(jump_bottom) || isnan(jump_top)
        return (
            device_length = length,
            device_width = width,
            current_bulk = NaN,
            jump_bottom = NaN,
            jump_top = NaN,
            a0_bulk_bottom = bottom_inner_bulk,
            a0_bulk_top = top_inner_bulk,
            r_contact_bottom = NaN,
            r_contact_top = NaN,
            r_contact_total = NaN,
            r_bulk = NaN,
            r_total = NaN,
        )
    end

    r_contact_bottom = jump_bottom / current_bulk
    r_contact_top = jump_top / current_bulk
    r_contact_total = (jump_bottom + jump_top) / current_bulk
    r_bulk = abs(bottom_inner_bulk - top_inner_bulk) / current_bulk
    r_total = bias / current_bulk

    return (
        device_length = length,
        device_width = width,
        current_bulk = current_bulk,
        jump_bottom = jump_bottom,
        jump_top = jump_top,
        a0_bulk_bottom = bottom_inner_bulk,
        a0_bulk_top = top_inner_bulk,
        r_contact_bottom = r_contact_bottom,
        r_contact_top = r_contact_top,
        r_contact_total = r_contact_total,
        r_bulk = r_bulk,
        r_total = r_total,
    )
end

function write_summary_csv(path, rows)
    open(path, "w") do io
        write(io, "device,mesh_file,output_h5,length,width,current_bulk,a0_bulk_bottom,a0_bulk_top,jump_bottom,jump_top,Rc_bottom,Rc_top,Rc_total,R_bulk,R_total\n")
        for row in rows
            write(
                io,
                string(
                    row.device, ",",
                    row.mesh_file, ",",
                    row.output_h5, ",",
                    row.length, ",",
                    row.width, ",",
                    row.current_bulk, ",",
                    row.a0_bulk_bottom, ",",
                    row.a0_bulk_top, ",",
                    row.jump_bottom, ",",
                    row.jump_top, ",",
                    row.r_contact_bottom, ",",
                    row.r_contact_top, ",",
                    row.r_contact_total, ",",
                    row.r_bulk, ",",
                    row.r_total, "\n",
                ),
            )
        end
    end
end

function print_contact_resistance_summary(rows)
    println("\n================ Contact Resistance Summary ================")
    for row in rows
        println(
            "device=$(row.device), L=$(round(row.length, digits=4)), ",
            "Rc_total=$(round(row.r_contact_total, digits=6)), ",
            "Rc_bottom=$(round(row.r_contact_bottom, digits=6)), ",
            "Rc_top=$(round(row.r_contact_top, digits=6)), ",
            "R_bulk=$(round(row.r_bulk, digits=6)), ",
            "R_total=$(round(row.r_total, digits=6))",
        )
    end

    finite_rows = filter(row -> isfinite(row.r_contact_total), rows)
    if length(finite_rows) >= 2
        short = first(finite_rows)
        long = last(finite_rows)
        delta = long.r_contact_total - short.r_contact_total
        pct = short.r_contact_total != 0 ? 100 * delta / short.r_contact_total : NaN
        println("\nInterpretation:")
        println(
            "- Rc_total changes from $(round(short.r_contact_total, digits=6)) at L=$(round(short.length, digits=4)) ",
            "to $(round(long.r_contact_total, digits=6)) at L=$(round(long.length, digits=4))."
        )
        println("- Absolute change: $(round(delta, digits=6)); percent change: $(round(pct, digits=3))%.")
        if delta > 0
            println("- Contact resistance increases with length for these two devices.")
        elseif delta < 0
            println("- Contact resistance decreases with length for these two devices.")
        else
            println("- Contact resistance is unchanged across the two lengths.")
        end
        println("- Compare Rc_total vs R_bulk to judge whether contacts or bulk dominate total resistance.")
    else
        println("\nInterpretation: insufficient finite data to compare Rc trends across device lengths.")
    end
    println("===========================================================\n")
end


# ======================================================================================================================
# Run Simulations and Analyze
# ======================================================================================================================

timestamp = Dates.format(now(), "yyyy-mm-dd_HHMMSS")
run_dir = joinpath(results_root, "run_$(timestamp)")
mkpath(run_dir)

summary_rows = NamedTuple[]

for device in devices
    @info "Running device" device=device.name mesh=basename(device.mesh_path) bias p_scatter gamma_mr gamma_mc
    boundary_conditions = make_boundary_conditions(bias, p_scatter)

    sol, semi = FermiHarmonics.solve(
        device.mesh_path,
        boundary_conditions,
        params,
        gamma_mr,
        gamma_mc;
        max_harmonic = max_harmonic,
        visualize = visualize,
        name = "simple_geometries_$(device.name)",
    )

    file_params = (
        device = device.name,
        bias = bias,
        p_scatter = p_scatter,
        gamma_mr = gamma_mr,
        gamma_mc = gamma_mc,
    )
    output_h5 = joinpath(run_dir, "observables_" * DrWatson.savename(file_params, "h5"))
    save_for_analysis(sol, semi, output_h5; nvisnodes = nvisnodes)

    metrics = analyze_contact_resistance(output_h5, bias)
    if isnan(metrics.r_contact_total)
        @warn "Contact-resistance analysis returned NaNs" device=device.name output_h5
    end

    push!(
        summary_rows,
        (
            device = device.name,
            mesh_file = basename(device.mesh_path),
            output_h5 = basename(output_h5),
            length = metrics.device_length,
            width = metrics.device_width,
            current_bulk = metrics.current_bulk,
            a0_bulk_bottom = metrics.a0_bulk_bottom,
            a0_bulk_top = metrics.a0_bulk_top,
            jump_bottom = metrics.jump_bottom,
            jump_top = metrics.jump_top,
            r_contact_bottom = metrics.r_contact_bottom,
            r_contact_top = metrics.r_contact_top,
            r_contact_total = metrics.r_contact_total,
            r_bulk = metrics.r_bulk,
            r_total = metrics.r_total,
        ),
    )

    @info "Device complete" device=device.name length=metrics.device_length Rc_total=metrics.r_contact_total R_total=metrics.r_total
end

sort!(summary_rows, by = row -> row.length)
summary_csv = joinpath(run_dir, "contact_resistance_summary.csv")
write_summary_csv(summary_csv, summary_rows)
print_contact_resistance_summary(summary_rows)

open(joinpath(run_dir, "run_info.txt"), "w") do io
    write(io, "timestamp=$(timestamp)\n")
    write(io, "bias=$(bias)\n")
    write(io, "p_scatter=$(p_scatter)\n")
    write(io, "gamma_mr=$(gamma_mr)\n")
    write(io, "gamma_mc=$(gamma_mc)\n")
    write(io, "nvisnodes=$(nvisnodes)\n")
end

@info "All runs complete" run_dir summary_csv
