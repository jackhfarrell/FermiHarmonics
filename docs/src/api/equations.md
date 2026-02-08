# Equations
We have implemented a custom equations type for use with Trixi. This offers a more fine-tuned way to run simulations that does not go through the `FermiHarmonics.solve()` function.

## Type

```@docs
FermiHarmonics2D
```

## Constructor

```@docs
FermiHarmonics.FermiHarmonics2D(nvars::Integer; gamma_mr, gamma_mc, max_harmonic=0)
```

## Script-Style Usage

```julia
params = SolveParams(;
    max_harmonic = 60,
    # ... other solver settings
)

nvars = 1 + 2 * params.max_harmonic
equations = FermiHarmonics2D(
    nvars;
    gamma_mr = 0.05,
    gamma_mc = 0.4,
    max_harmonic = params.max_harmonic,
)
```

Use `:auto` at solve-time, e.g. `solve(...; max_harmonic=:auto)`.
