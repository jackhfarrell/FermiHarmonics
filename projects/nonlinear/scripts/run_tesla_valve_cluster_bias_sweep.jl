using FermiHarmonics

function tesla_cluster_mesh_has_nodesets(mesh_path::AbstractString, names)
    isfile(mesh_path) || return false
    contents = read(mesh_path, String)
    return all(occursin("*NSET,NSET=$(name)", contents) for name in names)
end

function tesla_cluster_mesh_has_quad_elements(mesh_path::AbstractString)
    isfile(mesh_path) || return false
    contents = read(mesh_path, String)
    return occursin("type=CPS4", contents) || occursin("type=CPE4", contents)
end

function ensure_tesla_cluster_mesh(geo_path::AbstractString, mesh_path::AbstractString)
    required_nodesets = ("inlet", "outlet", "walls")
    needs_rebuild = !isfile(mesh_path) ||
                    mtime(mesh_path) < mtime(geo_path) ||
                    !tesla_cluster_mesh_has_nodesets(mesh_path, required_nodesets) ||
                    !tesla_cluster_mesh_has_quad_elements(mesh_path)
    needs_rebuild || return mesh_path

    mkpath(dirname(mesh_path))
    cmd = Cmd([
        "gmsh",
        "-2",
        geo_path,
        "-string",
        "Mesh.SaveGroupsOfNodes=1; Mesh.RecombineAll=1; Mesh.Algorithm=8;",
        "-format",
        "inp",
        "-o",
        mesh_path,
    ])
    @info "Building Tesla valve mesh" geo_path mesh_path
    run(cmd)
    return mesh_path
end

function tesla_cluster_bias_values()
    return collect(range(0.01, 0.5, length=30))
end

function tesla_cluster_direction_label(direction::AbstractString)
    direction in ("forward", "reverse") ||
        throw(ArgumentError("direction must be \"forward\" or \"reverse\""))
    return direction
end

function tesla_cluster_boundary_conditions(direction::AbstractString, bias::Real, p_scatter::Real)
    direction_label = tesla_cluster_direction_label(direction)
    half_bias = 0.5 * Float64(bias)
    inlet_bias, outlet_bias = direction_label == "forward" ? (half_bias, -half_bias) : (-half_bias, half_bias)
    return Dict(
        :walls => MaxwellWallBC(Float64(p_scatter)),
        :inlet => OhmicContactBC(inlet_bias),
        :outlet => OhmicContactBC(outlet_bias),
    )
end

function tesla_cluster_bias_label(bias::Real)
    return replace(string(round(Float64(bias); digits=5)), "." => "p")
end

function tesla_cluster_run_name(direction::AbstractString, bias::Real)
    return "tesla_valve_cluster_$(tesla_cluster_direction_label(direction))_gamma_mc100_M10_bias$(tesla_cluster_bias_label(bias))_chi000_poly3_t2"
end

function run_tesla_valve_cluster_bias_sweep(; direction::AbstractString=get(ENV, "TESLA_DIRECTION", "forward"),
                                             output_root::AbstractString=get(ENV, "TESLA_OUTPUT_DIR",
                                                                             joinpath(@__DIR__, "..", "data",
                                                                                      "tesla_valve_cluster_bias_sweep")))
    project_root = normpath(joinpath(@__DIR__, "..", "..", ".."))
    geo_path = joinpath(project_root, "projects", "nonlinearities", "mesh", "tesla_valve.geo")
    mesh_path = joinpath(project_root, "projects", "nonlinearities", "mesh", "tesla_valve.inp")
    ensure_tesla_cluster_mesh(geo_path, mesh_path)

    reference = FermiHarmonics.blg_reference_setup()
    mu0 = reference.mu0
    mass = reference.mass
    p_scatter = reference.p_scatter
    gamma_mr = 0.0
    gamma_mc = 100.0
    chi = 0.0
    max_harmonic = 10
    biases = tesla_cluster_bias_values()
    direction_label = tesla_cluster_direction_label(direction)
    output_dir = joinpath(output_root, direction_label)
    mkpath(output_dir)

    params = SolveParams(;
        polydeg=3,
        tspan_end=parse(Float64, get(ENV, "TESLA_TSPAN_END", "2.0")),
        residual_tol=1e-3,
        cfl=0.8,
        log_every=200,
        min_harmonic=4,
        max_harmonic_auto=20,
    )

    @info "Running Tesla valve cluster bias sweep" direction=direction_label gamma_mr gamma_mc chi max_harmonic nbias=length(biases) output_dir

    u0 = nothing
    for (bias_index, bias) in enumerate(biases)
        boundary_conditions = tesla_cluster_boundary_conditions(direction_label, bias, p_scatter)
        run_name = tesla_cluster_run_name(direction_label, bias)

        @info "Running Tesla valve sweep case" direction=direction_label bias bias_index total_biases=length(biases) warm_start=!isnothing(u0)

        sol, semi = FermiHarmonics.solve(
            mesh_path,
            boundary_conditions,
            params,
            gamma_mr,
            gamma_mc;
            transport=:parabolic_nonlinear,
            max_harmonic=max_harmonic,
            mu0=mu0,
            mass=mass,
            chi=chi,
            visualize=false,
            u0_override=u0,
            name=run_name,
        )

        status = FermiHarmonics.solve_status(sol, semi, params)
        cartesian_path = joinpath(output_dir, "$(run_name).h5")
        mesh_native_path = joinpath(output_dir, "$(run_name)_mesh_native.h5")
        FermiHarmonics.save_for_analysis(sol, semi, cartesian_path; nvisnodes=400, observables=[:n, :jx, :jy])
        FermiHarmonics.save_mesh_native_analysis(sol, semi, mesh_native_path; refine=6, observables=[:n, :jx, :jy])
        @info "Saved Tesla valve sweep outputs" direction=direction_label bias cartesian_path mesh_native_path final_time=sol.t[end] stop_reason=status.stop_reason final_residual=status.final_residual converged=status.converged
        if !status.converged
            @warn "Tesla valve sweep case did not reach residual tolerance before stopping" direction=direction_label bias final_time=status.final_time target_final_time=status.target_final_time final_residual=status.final_residual tolerance=params.residual_tol stop_reason=status.stop_reason
        end

        u0 = copy(sol.u[end])
    end

    return output_dir
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_tesla_valve_cluster_bias_sweep()
end
