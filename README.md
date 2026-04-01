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

## Core Solve Shape

The public API is now built around typed model composition in the core package
plus a Trixi extension for execution:

```julia
using ElectronKinetics
using Trixi

surface = Isotropic2DFermiSurface(; vF=1.0, nu=1.0, mass=1.0, charge=-1.0)
model = KineticModel2D(
    surface,
    HarmonicBasis(:auto),
    IsotropicHarmonicStreaming(),
    LinearBGKCollision(0.0, TwoRateProfile(50.0)),
)

problem = TrixiProblem(;
    # Convert the .geo to .inp first, e.g.:
    # gmsh -2 assets/square_bells.geo -format inp -o assets/square_bells.inp
    mesh_path = "assets/square_bells.inp",
    boundary_conditions = Dict(
        :walls => MaxwellWallBC(1.0),
        :contact_top => OhmicContactBC(-0.5),
        :contact_bottom => OhmicContactBC(0.5),
    ),
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

callbacks = default_callbacks_builder()

sol, semi = solve(
    problem,
    model,
    config;
    callbacks=callbacks,
    name="quick_start",
)
```

The physics core can be loaded without `Trixi`; the backend extension activates
when `using Trixi` is present in the environment.

## BLG Reference Convention

For nonlinear straight-channel studies, the repository now standardizes a
dimensionless BLG-oriented reference convention while keeping the solver API
dimensionless:

- `mu0 = 1`
- `vF = 1`
- `mass = 2`
- `gamma_mr = 0`
- `gamma_ee = 0`
- straight-channel total length `L = 1`

These values are available through `ElectronKinetics.blg_reference_setup()`.
The choice `mu0 = 1` and `vF = 1` implies `mass = 2` because the nonlinear
transport model uses `vF = sqrt(2 * mu0 / mass)`.

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
