# ElectronKinetics

Extensible 2D Fermi-liquid kinetic transport on unstructured meshes.

ElectronKinetics is a Julia framework for space-resolved electron flow in 2D kinetic theory. The v1 refactor separates a backend-agnostic physics core from the current `Trixi.jl` backend, so the package can grow toward richer Fermi surfaces, streaming operators, and collision models while keeping the existing isotropic production kernels as explicit fast paths. In `transport=:parabolic_nonlinear` mode, the solver combines nonlinear streaming with nonlinear BGK collision targets built from the exact drifting local-equilibrium manifold of the parabolic band.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18528662.svg)](https://doi.org/10.5281/zenodo.18528662)

Documentation: <https://fermiflows.jackhfarrell.com>

## Quick Start

### 1. Install Julia dependencies

From the repository root:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

### 2. Run an example

```bash
julia --project=. examples/01_linear_harmonic_transport.jl
```

Enable live visualization (requires GLMakie):

```bash
julia --project=. examples/01_linear_harmonic_transport.jl --live
```

Preview the mesh only (no solve):

```bash
julia --project=. examples/01_linear_harmonic_transport.jl --live --mesh-only
```

## Core Solve Shape

The public API is now built around typed model composition in the core package
plus a Trixi extension for execution:

```julia
using ElectronKinetics
using Trixi

surface = Isotropic2DFermiSurface(; fermi_velocity=1.0, nu=1.0, mass=1.0, charge=-1.0)
model = KineticModel2D(
    surface,
    HarmonicBasis(:auto),
    IsotropicHarmonicStreaming(),
    LinearBGKCollision(0.0, TwoRateProfile(50.0)),
)

problem = TrixiProblem(;
    # Convert the .geo to .inp first, or use a .geo directly with mesh_build.
    geometry_path = "assets/square_bells.geo",
    boundary_conditions = Dict(
        :walls => MaxwellWallBC(1.0),
        :contact_top => OhmicContactBC(-0.5),
        :contact_bottom => OhmicContactBC(0.5),
    ),
    mesh_build = MeshBuildConfig(mesh_scale=3.0),
)

config = SolverConfig(;
    polydeg = 3,
    cfl = 0.8,
    tspan_end = 50.0,
    residual_tol = 1e-5,
    log_every = 200,
    min_harmonic = 4,
    max_harmonic_auto = 150,
)

# MeshBuildConfig defaults: algorithm=8 (quasi-structured quads), mesh_scale=3.0.

callbacks = default_callbacks_builder()

sol, semi = solve(
    problem,
    model,
    config;
    callbacks=callbacks,
    name="quick_start",
)
```

Minimal custom callback example:

```julia
callbacks = (semi, ode, config, live_visualization, name) -> begin
    SciMLBase.CallbackSet(
        Trixi.StepsizeCallback(cfl=config.cfl),
    )
end
```

Note: mesh preview uses quad outlines from the `.inp`, while the live solver
dashboard visualizes a triangulated mesh-native grid.

## Practical I/O + Visualization Notes

Here are the most common “how do I actually use this?” workflows.

### Save Mesh-Native vs Cartesian Outputs

Mesh-native outputs follow the unstructured mesh (triangulated for visualization),
while cartesian outputs resample onto a uniform grid. Both are available:

```julia
using ElectronKinetics
using Trixi

# After solve(...)
save_mesh_native_analysis(sol, semi, "run_mesh_native.h5"; refine=4, observables=[:n, :jx, :jy])
save_for_analysis(sol, semi, "run_cartesian.h5"; nvisnodes=256, observables=[:n, :jx, :jy])
```

Tip: `save_mesh_native_analysis` uses a mesh-native triangulation even if the
original mesh is quad-based.

### Live Visualization (Configurable)

You can control the live dashboard layout and geometry using `LiveVisualizationConfig`:

```julia
live_visualization = LiveVisualizationConfig(;
    geometry_mode=:mesh_native, # or :cartesian
    field=:n,                   # :n, :a0, :jx, :jy, :current_magnitude, ...
    refine=3,                   # mesh-native refinement
    nvisnodes=160,              # cartesian grid resolution
    accepted_step_interval=50,
    min_update_seconds=0.2,
)

callbacks = (semi, ode, config, live_config, name) -> begin
    Trixi.CallbackSet(
        Trixi.StepsizeCallback(cfl=config.cfl),
        Trixi.SteadyStateCallback(abstol=config.residual_tol, reltol=config.residual_reltol),
        ElectronKinetics.visualization_callback(config, semi, name; interval=config.log_every, mode=live_visualization.geometry_mode),
    )
end
```

You can also use the built-in `default_callbacks_builder()` if you just want
the default live dashboard wiring.

### Magnetic Fields (Current Scope)

Uniform perpendicular magnetic fields are supported for **single-band linear
harmonic** runs:

```julia
magnetic_field = MagneticField2D(0.2)
model = KineticModel2D(
    surface,
    HarmonicBasis(:auto),
    IsotropicHarmonicStreaming(),
    LinearBGKCollision(0.0, TwoRateProfile(50.0));
    magnetic_field=magnetic_field,
)
```

If you need magnetic fields in multiband or nonlinear angle runs, we can add
that next.

The physics core can be loaded without `Trixi`; the backend extension activates
when `using Trixi` is present in the environment.

## BLG Reference Convention

For nonlinear straight-channel studies, the repository now standardizes a
dimensionless BLG-oriented reference convention while keeping the solver API
dimensionless:

- `mu0 = 1`
- `fermi_velocity = 1`
- `mass = 2`
- `gamma_mr = 0`
- `gamma_ee = 0`
- straight-channel total length `L = 1`

These values are available through `ElectronKinetics.blg_reference_setup()`.
The choice `mu0 = 1` and `fermi_velocity = 1` implies `mass = 2` because the nonlinear
transport model uses `fermi_velocity = sqrt(2 * mu0 / mass)`.

The physical BLG effective mass `m = 0.03 m_e` is treated as part of the
back-mapping to physical units, not as a direct low-level solver input.
Exact neutrality is also outside the current nonlinear parabolic-band solver,
which requires `mu0 > 0`.

## Citation
If you happen to find this code useful, it would be great if you would cite our upcoming theory/numerics paper as well as the codebase itself!

- Code release (Zenodo): [10.5281/zenodo.18528662](https://doi.org/10.5281/zenodo.18528662)
- Theory/application paper: Farrell & Lucas (2026, to appear)

### BibTeX (Code)

```bibtex
@software{fermiflows_zenodo,
  author = {Farrell, Jack H.},
  title = {{ElectronKinetics}: Extensible 2D Fermi-liquid kinetic transport on unstructured meshes},
  year = {2026},
  publisher = {Zenodo},
  doi = {10.5281/zenodo.18528662},
  url = {https://doi.org/10.5281/zenodo.18528662}
}
```

### BibTeX (Theory Paper)

```bibtex
@article{farrellSimpleDevices2026,
  title = {Simple devices that distinguish hydrodynamic, ballistic, and diffusive transport},
  author = {Farrell, Jack H. and Lucas, Andrew},
  year = {2026},
  journal = {to appear}
}
```
