# FermiHarmonics
This is a Julia code for simulating a toy model of Fermi liquid transport with momentum relaxing and momentum conserving collisions.  I use `Trixi.jl` which is a library for high-accuracy PDE solutions.

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
By default, we adaptively pick `M` based on how strong the damping is from collisions. The minimum `M` is `4` and the maximum `M` is `150`. This way, simulations in high damping regimes use fewer harmonics, while weakly damped cases keep higher angular resolution.

We adopt a relaxation-time-like (BGK) approximation for the collision integral, so that the Boltzmann equation in Harmonic basis reads, for $m=0$
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
For the scattering rates, we adopt a two-time (linearized BGK) model to capture the ballistic-hydrodynamic-diffusive crossover, 
```math
\gamma_0 = 0
```
```math
\gamma_1 = \gamma_{\mathrm{mr}}
```
```math
\gamma_n = \gamma_{\mathrm{mr}} + \gamma_{\mathrm{mc}}, (n \ge 2).
```

In this code release, source terms are purely physical; we do not apply additional numerical tail damping.

## Solve Entry Point

The main solve interface is documented in:

- [Solve API](api/solve.md)
- [Mesh Guide](mesh.md)

## Quick Start

```julia
using FermiHarmonics
using Trixi  # for boundary condition constructors

boundary_conditions = Dict(
    :walls => MaxwellWallBC(1.0),
    :contact_top => OhmicContactBC(-0.5),
    :contact_bottom => OhmicContactBC(0.5),
)

params = SolveParams(;
    polydeg = 3,
    cfl = 0.8,
    tspan_end = 50.0,
    residual_tol = 1e-5,
    log_every = 200,
    min_harmonic = 4,
    max_harmonic_auto = 150,
)

sol, semi = solve(
    "projects/square_bells_ucsb/mesh/square_bells.inp",
    boundary_conditions,
    params,
    0.0,   # gamma_mr
    50.0;  # gamma_mc
    max_harmonic = :auto,
    visualize = false,
    name = "quick_start",
)
```


## API Reference

```@contents
Pages = [
    "api/equations.md",
    "api/boundary_conditions.md",
    "api/source_terms.md",
    "api/io_utils.md",
    "api/solve.md",
    "api/slurm_utils.md",
]
Depth = 2
```

## Index

```@index
```
