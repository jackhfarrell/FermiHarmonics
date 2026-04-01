# Mesh

This project expects **2D quad meshes** in Abaqus `.inp` format at runtime (loaded with `Trixi.P4estMesh{2}`).
The editable geometry source of truth can now be a Gmsh `.geo` file generated into `.inp` directly from Julia.

Recommended workflow:

1. Define geometry + physical boundary names in a `.geo` file.
2. Generate a quad mesh from Julia with [`generate_mesh_from_geo`](@ref).
3. Let `solve` consume either the generated `.inp` or the `.geo` directly.
4. Use the same physical names in your Julia `boundary_conditions` dictionary.

`generate_mesh_from_geo` uses `Gmsh.jl` directly (not the `gmsh` CLI), so the
Gmsh Julia package must be available in your environment.

## Quad Requirement (Important)

Do not generate triangle-only meshes.

In your `.geo`, force recombination to quads, e.g.

```geo
Mesh.RecombineAll = 1;
Recombine Surface {1};
```

Optional (often helpful for quads):

```geo
Mesh.Algorithm = 8; // Frontal-Delaunay for quads
```

## Boundary Naming

There is no single required set of boundary names at the package level.
Use whatever physical names fit your geometry, as long as:

- they are defined as `Physical Curve("...")` in Gmsh, and
- the same names are used as keys in your Julia `boundary_conditions`.

For the square-bells example, the names are `walls`, `contact_top`, and `contact_bottom`:

```geo
Physical Surface("domain") = {1};
Physical Curve("contact_bottom") = {1};
Physical Curve("contact_top") = {11};
Physical Curve("walls") = {2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15, 16, 17, 18, 19, 20};
```

Reference file:

- `assets/square_bells.geo`

## Build `.inp` From `.geo` In Julia

From Julia:

```julia
using ElectronKinetics

mesh_path = generate_mesh_from_geo(
    "assets/square_bells.geo";
    config=MeshBuildConfig(mesh_scale=3.0),
)
```

`MeshBuildConfig(mesh_scale=3.0)` uses temporary unique output paths, which is safe for
cluster sweeps and parallel parameter scans. `.geo` meshing always writes a
unique temporary `.inp` (the `.geo` is the persistent source of truth).

### MeshBuildConfig Defaults

- `algorithm = 8` (quad-friendly Gmsh algorithm for quasi-structured quads)
- `recombine_all = true` (prefer quadrilateral elements)
- `mesh_scale = 3.0` (coarser meshes via `Mesh.CharacteristicLengthFactor`)

## Map Mesh Boundaries To Boundary Conditions

Boundary-condition keys must match physical curve names in the mesh:

```julia
using ElectronKinetics
using Trixi

boundary_conditions = Dict(
    :walls => MaxwellWallBC(1.0),
    :contact_top => OhmicContactBC(-0.5),
    :contact_bottom => OhmicContactBC(0.5),
)

surface = Isotropic2DFermiSurface()
model = KineticModel2D(
    surface,
    HarmonicBasis(:auto),
    IsotropicHarmonicStreaming(),
    LinearBGKCollision(0.0, TwoRateProfile(50.0)),
)

config = SolverConfig()

problem = TrixiProblem(;
    geometry_path = "assets/square_bells.geo",
    boundary_conditions = boundary_conditions,
)

callbacks = default_callbacks_builder()

sol, semi = solve(problem, model, config; callbacks=callbacks)
```

If names do not match, boundary assignment fails.

## Notes

- Keep physical names stable once sweeps start.

You can also keep using `mesh_path = "…/square_bells.inp"` for prebuilt meshes.
