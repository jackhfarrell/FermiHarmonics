# Boundary Conditions

Boundary conditions remain backend-facing in v1. The core package exposes the
physical boundary-condition types used by the Trixi extension, but does not yet
define a backend-neutral execution interface.

## Public Types

- `MaxwellWallBC`
- `OhmicContactBC`
