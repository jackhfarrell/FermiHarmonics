# Solve: main entry point

```@docs
FermiHarmonics.SolveParams
FermiHarmonics.estimate_max_harmonic
FermiHarmonics.solve
```

## Script-Style Usage

```julia
params = SolveParams(;
    max_harmonic = 60,   # optional fixed cap if you call solve(...; max_harmonic=params.max_harmonic)
    polydeg = 3,
    cfl = 0.2,
    tspan_end = 50.0,
    residual_tol = 1e-8,
    log_every = 50,
)

boundary_conditions = Dict(
    :walls => MaxwellWallBC(1.0),
    :source => OhmicContactBC(0.5),
    :drain => OhmicContactBC(-0.5),
)

sol, semi = solve(
    "mesh.inp",
    boundary_conditions,
    params,
    0.05,  # gamma_mr
    0.40;  # gamma_mc
    max_harmonic = :auto,
    visualize = false,
    name = "case_001",
)
```

## Behavior Notes

- `solve` requires a `SolveParams` instance for solver settings.
- Harmonic count is fixed per solve (`nvars = 1 + 2M`), but `M` can be selected automatically per case via `max_harmonic=:auto`.
- Auto mode uses `estimate_max_harmonic(gamma_mr, gamma_mc)` with `gamma_total = gamma_mr + gamma_mc`.
- Default logarithmic map is `gamma_total=0 -> M=150` and `gamma_total>=1000 -> M=4`.
- Optional tuning keys in `SolveParams` are
  `min_harmonic`, `max_harmonic_auto`, and `auto_harmonic_decay_power`.
- Warm starts via `u0_override` support harmonic-count changes by truncating or zero-padding higher modes.

## Auto Harmonic Selector

When `solve(...; max_harmonic=:auto)` is used, the code computes

```math
\gamma_{\mathrm{tot}} = \gamma_{\mathrm{mr}} + \gamma_{\mathrm{mc}},
```

then selects `M` using a conservative logarithmic rule:

- `\gamma_tot <= 1`: `M = max_harmonic_auto`
- `\gamma_tot >= 1000`: `M = min_harmonic`
- otherwise: shifted powered-log interpolation
  `s = (log1p(gamma_tot)-log1p(1))/(log1p(1000)-log1p(1))`,
  then `M ~ M_min + (M_max-M_min) * (1-s)^p`,
  rounded up to an integer

So with default settings (`min_harmonic=4`, `max_harmonic_auto=150`):

- `(\gamma_mr, \gamma_mc) = (0, 0)` gives `M = 150`
- `(\gamma_mr, \gamma_mc) = (0, 100)` gives `M = 34`
- `(\gamma_mr, \gamma_mc) = (0, 200)` gives `M = 21`
- `(\gamma_mr, \gamma_mc) = (0, 1000)` gives `M = 4`

You can tune the range/shape through
`SolveParams(min_harmonic=..., max_harmonic_auto=..., auto_harmonic_decay_power=...)`
without changing the solve workflow.
