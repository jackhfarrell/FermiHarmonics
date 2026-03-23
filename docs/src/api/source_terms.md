# Source Terms

The collision side of the kinetic equation is now routed through the typed
`collision_sources!` interface. Trixi backend equations provide specialized
implementations for the current isotropic fast paths.

The key public hooks are:

- `collision_sources!`
- `mode_rate`

For `transport = :parabolic_nonlinear`, the built-in harmonic closure uses the
nonlinear BGK analogue of the two-rate model: a momentum-relaxing piece toward
isotropic equilibrium and a momentum-conserving piece toward drifting local
equilibrium.
