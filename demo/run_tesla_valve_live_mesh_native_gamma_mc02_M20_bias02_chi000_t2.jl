using FermiHarmonics

function mesh_has_nodesets(mesh_path::AbstractString, names)
    isfile(mesh_path) || return false
    contents = read(mesh_path, String)
    return all(occursin("*NSET,NSET=$(name)", contents) for name in names)
end

function mesh_has_quad_elements(mesh_path::AbstractString)
    isfile(mesh_path) || return false
    contents = read(mesh_path, String)
    return occursin("type=CPS4", contents) || occursin("type=CPE4", contents)
end

function ensure_tesla_valve_mesh(geo_path::AbstractString, mesh_path::AbstractString)
    required_nodesets = ("inlet", "outlet", "walls")
    needs_rebuild = !isfile(mesh_path) ||
                    mtime(mesh_path) < mtime(geo_path) ||
                    !mesh_has_nodesets(mesh_path, required_nodesets) ||
                    !mesh_has_quad_elements(mesh_path)
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

function main()
    project_root = normpath(joinpath(@__DIR__, ".."))
    geo_path = joinpath(project_root, "projects", "nonlinearities", "mesh", "tesla_valve.geo")
    mesh_path = joinpath(project_root, "projects", "nonlinearities", "mesh", "tesla_valve.inp")
    output_dir = joinpath(@__DIR__, "data_nonlinear")
    mkpath(output_dir)
    ensure_tesla_valve_mesh(geo_path, mesh_path)

    reference = FermiHarmonics.blg_reference_setup()
    mu0 = reference.mu0
    mass = reference.mass
    gamma_mr = reference.gamma_mr
    gamma_mc = 0.2
    bias = 0.2
    chi = 0.0
    max_harmonic = 20
    p_scatter = reference.p_scatter

    boundary_conditions = Dict(
        :walls => MaxwellWallBC(p_scatter),
        :inlet => OhmicContactBC(bias / 2),
        :outlet => OhmicContactBC(-bias / 2),
    )

    params = SolveParams(;
        polydeg=3,
        tspan_end=2.0,
        residual_tol=1e-3,
        cfl=0.8,
        log_every=200,
        min_harmonic=4,
        max_harmonic_auto=20,
    )

    run_name = "tesla_valve_live_mesh_native_gamma_mc02_M20_bias02_chi000_poly3_t2"

    @info "Running Tesla valve nonlinear live mesh-native case" mu0 mass gamma_mr gamma_mc bias chi max_harmonic polydeg=params.polydeg tspan_end=params.tspan_end residual_tol=params.residual_tol

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
        visualize=true,
        visualize_every=2000,
        visualization_mode=:mesh_native,
        name=run_name,
    )

    cartesian_path = joinpath(output_dir, "$(run_name).h5")
    mesh_native_path = joinpath(output_dir, "$(run_name)_mesh_native.h5")
    FermiHarmonics.save_for_analysis(sol, semi, cartesian_path)
    FermiHarmonics.save_mesh_native_analysis(sol, semi, mesh_native_path; refine=6)
    @info "Saved Tesla valve nonlinear analysis output" cartesian_path mesh_native_path final_time=sol.t[end]

    return nothing
end

main()
