# Model Construction
FermiFlows v1 builds the public physics API from composable typed objects in the
core package. The Trixi extension lowers these model objects to backend-specific
equation adapters internally; those adapters are no longer the primary public
API.

## Core Types

- `AbstractFermiSurface2D`
- `AbstractAngularDiscretization2D`
- `AbstractStreamingOperator2D`
- `AbstractModeRateProfile`
- `AbstractCollisionModel2D`
- `Isotropic2DFermiSurface`
- `HarmonicBasis`
- `AngleGrid`
- `IsotropicHarmonicStreaming`
- `IsotropicAngleStreaming`
- `TwoRateProfile`
- `OddQuarticRateProfile`
- `ConstantModeRateProfile`
- `CustomModeRateProfile`
- `LinearBGKCollision`
- `QuadraticBGKCollision`
- `ExactAngleBGKCollision`
- `TwoRateAngleBGKCollision`
- `KineticModel2D`
- `BandSpec`

## Example

```julia
surface = Isotropic2DFermiSurface(; vF=1.0, nu=1.0, mass=1.0, charge=-1.0)
discretization = HarmonicBasis(:auto)
streaming = IsotropicHarmonicStreaming()
collision = LinearBGKCollision(0.05, TwoRateProfile(0.40))

model = KineticModel2D(
    surface,
    discretization,
    streaming,
    collision;
    reference = blg_reference_setup(),
)
```

Multiband linear models are constructed by passing a `bands=[BandSpec(...), ...]`
vector to `KineticModel2D`. In v1, multiband support is intentionally limited to
the existing two-band linear harmonic workflow.
