# SLURM Utilities API

SLURM and sweep-management helpers are not part of the new public `ElectronKinetics`
v1 API. Existing project-specific scripts can still carry local automation, but
the package interface is now focused on:

- backend-agnostic model construction in `src/ElectronKinetics.jl`
- backend execution through Julia package extensions such as `ElectronKineticsTrixiExt`
- analysis and export helpers that operate on solved backend problems

## Behavior Notes

- Project-level SLURM automation should live beside the project that owns it.
- If these helpers are promoted back into the package later, they should follow
  the same backend-extension pattern as the solver interface.
