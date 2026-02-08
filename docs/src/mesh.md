# Mesh

This project expects **2D quad meshes** in Abaqus `.inp` format (loaded with `Trixi.P4estMesh{2}`).

Workflow:

1. Define geometry + physical boundary names in a `.geo` file.
2. Force quad meshing in Gmsh.
3. Export to `.inp`.
4. Use those same names in your Julia `boundary_conditions` dictionary.

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

- `projects/square_bells_ucsb/mesh/square_bells.geo`

## Build `.inp` From `.geo` With Gmsh

From repository root (quad-focused command):

```bash
gmsh -2 projects/square_bells_ucsb/mesh/square_bells.geo \
  -string "Mesh.RecombineAll=1;" \
  -format inp \
  -o projects/square_bells_ucsb/mesh/square_bells.inp
```

## Map Mesh Boundaries To Boundary Conditions

Boundary-condition keys must match physical curve names in the mesh:

```julia
using FermiHarmonics
using Trixi

boundary_conditions = Dict(
    :walls => MaxwellWallBC(1.0),
    :contact_top => OhmicContactBC(-0.5),
    :contact_bottom => OhmicContactBC(0.5),
)

params = SolveParams()

sol, semi = solve(
    "projects/square_bells_ucsb/mesh/square_bells.inp",
    boundary_conditions,
    params,
    0.0,   # gamma_mr
    50.0;  # gamma_mc
    max_harmonic = :auto,
)
```

If names do not match, boundary assignment fails.

## Notes

- Keep physical names stable once sweeps start.
- When running on SLURM, this project copies the mesh to node-local scratch to avoid parallel file-system races.

## TODO

- Add a Julia-side mesh generation helper so `.geo -> .inp` can be called directly from Julia scripts (instead of shelling out to `gmsh` manually).
