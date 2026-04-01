# Solve

## Public Entry Points

- `SolverConfig` stores numerical settings such as polynomial degree, CFL,
  final time, residual tolerance, and auto-harmonic bounds.
- `MeshBuildConfig` stores Julia-side Gmsh meshing settings for `.geo -> .inp`
  generation.
- `TrixiProblem` stores either a mesh path or a geometry path plus the
  boundary-condition map passed to the Trixi backend.
- `estimate_max_harmonic(gamma_mr, gamma_ee; min_harmonic, max_harmonic)` is
  the conservative helper used by `HarmonicBasis(:auto)`.
- `solve(problem, model, config; kwargs...)` is provided by the Trixi extension.

## Trixi-Backed Usage

```julia
using ElectronKinetics
using Trixi

config = SolverConfig(;
    polydeg = 3,
    cfl = 0.2,
    tspan_end = 50.0,
    residual_tol = 1e-8,
    log_every = 50,
)

model = KineticModel2D(
    Isotropic2DFermiSurface(; fermi_velocity=1.0, nu=1.0, mass=1.0, charge=-1.0),
    HarmonicBasis(:auto),
    IsotropicHarmonicStreaming(),
    LinearBGKCollision(0.05, TwoRateProfile(0.40)),
)

problem = TrixiProblem(;
    geometry_path = "mesh.geo",
    boundary_conditions = Dict(
        :walls => MaxwellWallBC(1.0),
        :source => OhmicContactBC(0.5),
        :drain => OhmicContactBC(-0.5),
    ),
)

callbacks = default_callbacks_builder()

sol, semi = solve(problem, model, config; callbacks=callbacks, name = "case_001")
```

For nonlinear harmonic runs:

```julia
model = KineticModel2D(
    Isotropic2DFermiSurface(; fermi_velocity=1.0, nu=1.0, mass=2.0, charge=-1.0),
    HarmonicBasis(:auto),
    IsotropicHarmonicStreaming(),
    QuadraticBGKCollision(0.05, OddQuarticRateProfile(0.40); mu0=1.0, mass=2.0),
)

sol, semi = solve(problem, model, config; callbacks=callbacks, name = "quadratic_bgk_case")
```

The explicit angle-grid reference solver stays available too:

```julia
model = KineticModel2D(
    Isotropic2DFermiSurface(; fermi_velocity=1.0, nu=1.0, mass=2.0, charge=-1.0),
    AngleGrid(128),
    IsotropicAngleStreaming(),
    ExactAngleBGKCollision(; gamma_mr=0.05, gamma_ee=0.40, mu0=1.0, mass=2.0),
)

sol, semi = solve(problem, model, config; callbacks=callbacks, name = "exact_bgk_case")
```

## Behavior Notes

- `ElectronKinetics.solve` is provided by the Trixi extension, not by the core-only package load.
- `TrixiProblem` accepts exactly one of `mesh_path` or `geometry_path`.
- Legacy `solve("mesh.geo", boundary_conditions, ...)` is supported and generates a temporary unique `.inp` automatically.
- `callbacks` is required; use `default_callbacks_builder()` or provide your own SciML callback builder.
- The default builder wires logging/live visualization and finalizes the dashboard on completion.
- `SolverConfig` controls discretization order, CFL, end time, residual target, logging cadence, and auto-harmonic bounds.
- `MeshBuildConfig()` uses temporary unique output paths under `SLURM_TMPDIR` / `TMPDIR` / `tempdir()`.
- `.geo` meshing always writes a unique temporary `.inp` (the `.geo` is the persistent input).
- `HarmonicBasis(:auto)` resolves its working harmonic count from `estimate_max_harmonic(gamma_mr, gamma_ee; ...)`.
- Linear harmonic closures use `LinearBGKCollision(gamma_mr, profile)`.
- Nonlinear harmonic closures use `QuadraticBGKCollision(gamma_mr, profile; mu0, mass, ...)`.
- Angle-grid reference runs use `ExactAngleBGKCollision` or `TwoRateAngleBGKCollision`.
- Warm starts remain available through `u0_override` on the extension solve entry point.

## Auto Harmonic Selector

When a model uses `HarmonicBasis(:auto)`, the extension computes

```math
\gamma_{\mathrm{tot}} = \gamma_{\mathrm{mr}} + \gamma_{\mathrm{ee}},
```

then selects `M` using a conservative logarithmic rule:

- `\gamma_tot <= 1`: `M = max_harmonic_auto`
- `\gamma_tot >= 300`: `M = min_harmonic`
- otherwise: logarithmic interpolation between the configured bounds

So with the default v1 settings:

- `(\gamma_mr, \gamma_ee) = (0, 0)` gives `M = 100`
- `(\gamma_mr, \gamma_ee) = (0, 50)` selects an intermediate `M`
- `(\gamma_mr, \gamma_ee) = (0, 300)` gives `M = 4`

You can tune the auto selector through `SolverConfig(; min_harmonic=..., max_harmonic_auto=...)`
without changing the model structure.
