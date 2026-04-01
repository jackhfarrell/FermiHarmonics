using UUIDs

const MESH_PROVENANCE_LOCK = ReentrantLock()
const MESH_PROVENANCE = Dict{String, NamedTuple}()
const ACTIVE_TEMP_MESH_DIRS = String[]

@inline normalize_mesh_path(path::AbstractString) = abspath(String(path))

@inline function mesh_build_summary(config::MeshBuildConfig)
    return (
        recombine_all = config.recombine_all,
        algorithm = config.algorithm,
        save_groups_of_nodes = config.save_groups_of_nodes,
        output_mode = config.output_mode,
        output_dir = isnothing(config.output_dir) ? nothing : abspath(config.output_dir),
        prefix = config.prefix,
        gmsh_options = Dict{String, Float64}(string(name) => Float64(value) for (name, value) in config.gmsh_options),
    )
end

function register_mesh_provenance!(mesh_path::AbstractString;
                                   geometry_path::Union{Nothing, AbstractString}=nothing,
                                   mesh_build::Union{Nothing, MeshBuildConfig}=nothing)
    normalized_mesh_path = normalize_mesh_path(mesh_path)
    normalized_geometry_path = isnothing(geometry_path) ? nothing : normalize_mesh_path(geometry_path)
    summary = isnothing(mesh_build) ? nothing : mesh_build_summary(mesh_build)

    lock(MESH_PROVENANCE_LOCK) do
        MESH_PROVENANCE[normalized_mesh_path] = (
            mesh_path = normalized_mesh_path,
            geometry_path = normalized_geometry_path,
            mesh_build = summary,
        )
    end

    return normalized_mesh_path
end

function mesh_provenance(mesh_path::AbstractString)
    normalized_mesh_path = normalize_mesh_path(mesh_path)
    lock(MESH_PROVENANCE_LOCK) do
        return get(MESH_PROVENANCE, normalized_mesh_path, nothing)
    end
end

function gmsh_options_string(options::AbstractDict{String, Float64})
    isempty(options) && return ""
    ordered = sort!(collect(pairs(options)); by=first)
    return join(["$(key)=$(value)" for (key, value) in ordered], ",")
end

function mesh_provenance_attributes(mesh_path::AbstractString)
    provenance = mesh_provenance(mesh_path)
    attrs = Dict{String, Any}(
        "mesh_file" => basename(mesh_path),
    )
    isnothing(provenance) && return attrs

    if !isnothing(provenance.geometry_path)
        attrs["mesh_source_geometry"] = basename(provenance.geometry_path)
        attrs["mesh_source_geometry_path"] = provenance.geometry_path
    end

    if !isnothing(provenance.mesh_build)
        mesh_build = provenance.mesh_build
        attrs["mesh_build_recombine_all"] = Int(mesh_build.recombine_all)
        attrs["mesh_build_algorithm"] = mesh_build.algorithm
        attrs["mesh_build_save_groups_of_nodes"] = Int(mesh_build.save_groups_of_nodes)
        attrs["mesh_build_output_mode"] = String(mesh_build.output_mode)
        attrs["mesh_build_prefix"] = mesh_build.prefix
        if !isnothing(mesh_build.output_dir)
            attrs["mesh_build_output_dir"] = mesh_build.output_dir
        end
        attrs["mesh_build_gmsh_options"] = gmsh_options_string(mesh_build.gmsh_options)
    end

    return attrs
end

function mesh_scratch_root()
    return get(ENV, "SLURM_TMPDIR", get(ENV, "TMPDIR", tempdir()))
end

function mesh_output_path(geo_path::AbstractString, config::MeshBuildConfig)
    stem = splitext(basename(geo_path))[1]
    if config.output_mode === :temporary
        temp_dir = mktempdir(mesh_scratch_root(); prefix=config.prefix)
        push!(ACTIVE_TEMP_MESH_DIRS, temp_dir)
        return joinpath(temp_dir, "$(stem)_$(uuid4()).inp")
    end

    target_dir = isnothing(config.output_dir) ? dirname(geo_path) : abspath(config.output_dir)
    mkpath(target_dir)
    return joinpath(target_dir, "$(stem).inp")
end

function mesh_nset_names(mesh_path::AbstractString)
    names = Set{String}()
    for line in eachline(mesh_path)
        startswith(uppercase(line), "*NSET") || continue
        match_result = match(r"NSET\s*=\s*([^,\s]+)"i, line)
        isnothing(match_result) || push!(names, match_result.captures[1])
    end
    return names
end

function mesh_has_quad_elements(mesh_path::AbstractString)
    contents = read(mesh_path, String)
    return occursin("type=CPS4", contents) || occursin("type=CPE4", contents)
end

function mesh_has_domain_or_surface_nodeset(nset_names::Set{String}, required_boundary_names::Set{String})
    for name in nset_names
        lowercase_name = lowercase(name)
        if lowercase_name in ("domain", "fluid") || startswith(lowercase_name, "surface")
            return true
        end
    end
    return any(!(name in required_boundary_names) for name in nset_names)
end

function validate_generated_mesh(mesh_path::AbstractString; required_boundary_names=String[])
    isfile(mesh_path) || throw(ArgumentError("generated mesh file not found: $mesh_path"))
    mesh_has_quad_elements(mesh_path) ||
        throw(ArgumentError("generated mesh is missing quad elements (expected CPS4 or CPE4): $mesh_path"))

    nset_names = mesh_nset_names(mesh_path)
    isempty(nset_names) && throw(ArgumentError("generated mesh is missing *NSET definitions: $mesh_path"))

    required_names_set = Set(String.(required_boundary_names))
    missing = sort!(collect(setdiff(required_names_set, nset_names)))
    isempty(missing) || throw(ArgumentError(
        "generated mesh is missing required boundary node sets $(join(missing, ", ")): $mesh_path",
    ))

    mesh_has_domain_or_surface_nodeset(nset_names, required_names_set) ||
        throw(ArgumentError("generated mesh is missing a domain/surface nodeset: $mesh_path"))

    return mesh_path
end

const GMSH_UUID = UUID("705231aa-382f-11e9-3f0c-b7cb4346fdeb")

function gmsh_module()
    pkg_id = Base.PkgId(GMSH_UUID, "Gmsh")
    Base.require(pkg_id)
    return Base.loaded_modules[pkg_id]
end

function apply_gmsh_options!(config::MeshBuildConfig)
    gmsh = gmsh_module().gmsh
    Base.invokelatest(gmsh.option.setNumber, "Mesh.SaveGroupsOfNodes", config.save_groups_of_nodes ? 1.0 : 0.0)
    Base.invokelatest(gmsh.option.setNumber, "Mesh.RecombineAll", config.recombine_all ? 1.0 : 0.0)
    Base.invokelatest(gmsh.option.setNumber, "Mesh.Algorithm", Float64(config.algorithm))
    for (name, value) in config.gmsh_options
        Base.invokelatest(gmsh.option.setNumber, name, Float64(value))
    end
    return nothing
end

function synchronize_gmsh_model!()
    gmsh = gmsh_module().gmsh
    try
        Base.invokelatest(gmsh.model.geo.synchronize)
    catch
    end
    try
        Base.invokelatest(gmsh.model.occ.synchronize)
    catch
    end
    return nothing
end

function generate_mesh_from_geo(geo_path::AbstractString; config::MeshBuildConfig=MeshBuildConfig())
    validated_config = validate(config)
    normalized_geo_path = normalize_mesh_path(geo_path)
    isfile(normalized_geo_path) || throw(ArgumentError("geometry file not found: $normalized_geo_path"))

    if validated_config.output_mode !== :temporary
        @warn "MeshBuildConfig.output_mode=$(validated_config.output_mode) is ignored; .geo meshing always writes a unique temporary .inp."
        validated_config = MeshBuildConfig(;
            recombine_all = validated_config.recombine_all,
            algorithm = validated_config.algorithm,
            save_groups_of_nodes = validated_config.save_groups_of_nodes,
            output_mode = :temporary,
            output_dir = validated_config.output_dir,
            prefix = validated_config.prefix,
            gmsh_options = validated_config.gmsh_options,
        )
    end

    output_path = mesh_output_path(normalized_geo_path, validated_config)
    mkpath(dirname(output_path))

    gmsh = gmsh_module().gmsh
    Base.invokelatest(gmsh.initialize, String[], false, false)
    try
        apply_gmsh_options!(validated_config)
        Base.invokelatest(gmsh.open, normalized_geo_path)
        synchronize_gmsh_model!()
        Base.invokelatest(gmsh.model.mesh.generate, 2)
        Base.invokelatest(gmsh.write, output_path)
    finally
        Base.invokelatest(gmsh.finalize)
    end

    validate_generated_mesh(output_path)
    register_mesh_provenance!(output_path; geometry_path=normalized_geo_path, mesh_build=validated_config)
    return output_path
end

function resolve_mesh_path(mesh_path::AbstractString,
                           boundary_conditions::AbstractDict{Symbol, <:Any};
                           mesh_build::MeshBuildConfig=MeshBuildConfig())
    normalized_path = normalize_mesh_path(mesh_path)
    lowercase_path = lowercase(normalized_path)

    if endswith(lowercase_path, ".geo")
        resolved_mesh_path = generate_mesh_from_geo(normalized_path; config=mesh_build)
        validate_generated_mesh(resolved_mesh_path; required_boundary_names=string.(keys(boundary_conditions)))
        return resolved_mesh_path
    elseif endswith(lowercase_path, ".inp")
        isfile(normalized_path) || error("Mesh file not found: $normalized_path")
        register_mesh_provenance!(normalized_path)
        return normalized_path
    end

    throw(ArgumentError("mesh path must point to a .inp or .geo file: $normalized_path"))
end

function resolve_mesh_path(problem::TrixiProblem)
    validate(problem)
    source_path = isnothing(problem.geometry_path) ? something(problem.mesh_path) : problem.geometry_path
    return resolve_mesh_path(source_path, problem.boundary_conditions; mesh_build=problem.mesh_build)
end
