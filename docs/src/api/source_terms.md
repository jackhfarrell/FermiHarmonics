# Source Terms

## Combined Source Terms

```@docs
FermiHarmonics.source_terms
```

## Physical Source Terms

```@docs
FermiHarmonics.physical_sources
```

For `transport = :parabolic_nonlinear`, `physical_sources` uses the nonlinear BGK analogue of the
two-rate model: a momentum-relaxing piece toward isotropic equilibrium and a momentum-conserving piece
toward exact drifting local equilibrium.
