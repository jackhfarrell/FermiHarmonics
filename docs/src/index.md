# ElectronKinetics
ElectronKinetics is a Julia framework for simulating space-resolved electron flow in
2D kinetic theory. The package now exposes a backend-agnostic physics core and
loads the current `Trixi.jl` solver integration through a Julia package
extension.

## Model Details
We solve a linearized Boltzmann equation:
```math
\partial_t \phi + v_F(\cos\theta\,\partial_x+\sin\theta\,\partial_y)\phi
= \mathcal{C}[\phi].
```
The right hand side is the collision integral, featuring physical terms designed to mitigate Gibbs phenomena. We approximate the distribution function `\phi` as a Fourier series with maximum harmonic `M`:
```math
\phi(x,y,\theta,t)=\frac{a_0}{2}+\sum_{m=1}^{M}\left[a_m\cos(m\theta)+b_m\sin(m\theta)\right].
```
For harmonic discretizations, the solver can adaptively pick `M` based on the
strength of the collision rates. High-damping regimes use fewer harmonics while
weakly damped regimes retain higher angular resolution.

For linear transport we adopt a relaxation-time-like (BGK) approximation for the collision integral, so that the Boltzmann equation in harmonic basis reads, for $m=0$
```math
\partial_t a_0 + v_F\left(\partial_x a_1 + \partial_y b_1\right) = -\gamma_0 a_0,
```
and for $m \ge 1$
```math
\partial_t a_m
 + \frac{v_F}{2}\partial_x(a_{m-1} + a_{m+1})
 + \frac{v_F}{2}\partial_y(b_{m+1} - b_{m-1})
 = -\gamma_m a_m,
```
```math
\partial_t b_m
 + \frac{v_F}{2}\partial_x(b_{m-1} + b_{m+1})
 + \frac{v_F}{2}\partial_y(a_{m-1} - a_{m+1})
 = -\gamma_m b_m.
```
For the scattering rates, the default linear model is a two-time BGK closure,
 
```math
\gamma_0 = 0
```
```math
\gamma_1 = \gamma_{\mathrm{mr}}
```
```math
\gamma_n = \gamma_{\mathrm{mr}} + \gamma_{\mathrm{mc}}, (n \ge 2).
```

In v1 the collision closure is factored through typed collision models and
mode-rate profiles. Built-in profiles include a two-rate BGK model, an
odd-quartic nonlinear damping profile, a constant `gamma(m)`, and typed custom
closures.

For `transport = :parabolic_nonlinear`, the streaming term is evaluated from the exact parabolic-band flux
and the collision term relaxes toward local equilibrium on the drifting Fermi-disk manifold. The nonlinear
solver reports physical observables `n`, `jx`, and `jy` in analysis output.

## BLG Reference Convention

For nonlinear straight-channel studies, the canonical dimensionless reference
case is:

- `mu0 = 1`
- `vF = 1`
- `mass = 2`
- `gamma_mr = 0`
- `gamma_mc = 0`
- total channel length `L = 1`

This convention is exposed through `ElectronKinetics.blg_reference_setup()`.
Because the nonlinear parabolic-band transport uses
`vF = sqrt(2 * mu0 / mass)`, setting `mu0 = 1` and `vF = 1` fixes the solver
mass to `2`.

The intended physical interpretation is BLG-inspired, but the solver remains
dimensionless. Physical values like `m = 0.03 m_e` belong to the back-mapping
layer in analysis and notes, not to the low-level solver API. Exact
neutrality is also out of scope for the current nonlinear parabolic-band
transport, which requires `mu0 > 0`.

## Solve Entry Point

The main public interfaces are documented in:

- [Solve API](api/solve.md)
- [Model API](api/equations.md)
- [Mesh Guide](mesh.md)

## Quick Start

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

config = SolverConfig(;
    polydeg = 3,
    cfl = 0.8,
    tspan_end = 50.0,
    residual_tol = 1e-5,
    log_every = 200,
    min_harmonic = 4,
    max_harmonic_auto = 150,
)

problem = TrixiProblem(;
    mesh_path = "projects/square_bells_ucsb/mesh/square_bells.inp",
    boundary_conditions = Dict(
        :walls => MaxwellWallBC(1.0),
        :contact_top => OhmicContactBC(-0.5),
        :contact_bottom => OhmicContactBC(0.5),
    ),
)

sol, semi = solve(problem, model, config; visualize = false, name = "quick_start")
```


## API Reference

```@contents
Pages = [
    "api/equations.md",
    "api/boundary_conditions.md",
    "api/source_terms.md",
    "api/io_utils.md",
    "api/solve.md",
]
Depth = 2
```

## Index

```@index
```
