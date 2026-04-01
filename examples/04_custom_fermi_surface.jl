# # Custom Fermi Surface
#
# This example demonstrates how to define a custom Fermi surface
# with user-supplied velocity function fermi_velocity(θ).
#
# Custom surfaces are useful for:
# - Systems with known analytic formulas (elliptic, hexagonal, etc.)
# - Arbitrary band structures
# - Comparing theoretical vs. experimental dispersions

using ElectronKinetics
using Trixi

live = "--live" in ARGS
if live
    try
        @eval using GLMakie
    catch err
        error("Live visualization requested, but GLMakie is not available: $(err)")
    end
end

# Example 1: Elliptic Fermi surface (built-in analytic type)
# Anisotropy: fermi_velocity(θ) = fermi_velocity0 / √(cos²θ + sin²θ/a²)
# aspect ratio a > 1 means velocity higher in y-direction
elliptic_surface = EllipticFermiSurface2D(
    fermi_velocity0 = 1.0,           # Fermi velocity (isotropic limit)
    aspect = 2.0,        # Anisotropy: v_y/v_x ratio
    nu = 1.0,            # density of states
    mass = 1.0,
    charge = -1.0
)

# Example 2: Custom user-defined surface with arbitrary function
# Must provide: fermi_velocity_func, max_fermi_velocity (upper bound for CFL), and properties
custom_surface = GeneralFermiSurface2D(
    # fermi_velocity(θ) function: can be any formula
    # Here: modulated with 4-fold symmetry (square-like band)
    θ -> 1.0 * (1.0 + 0.3 * cos(4θ)),
    name = :square_modulated,
    max_fermi_velocity = 1.3,        # Maximum of fermi_velocity(θ) — MUST be tight upper bound
    nu = 1.0,
    mass = 1.0,
    charge = -1.0
)

# Use custom surface in a model
model = KineticModel2D(
    custom_surface,
    HarmonicBasis(:auto),
    IsotropicHarmonicStreaming(),
    LinearBGKCollision(0.1, TwoRateProfile(0.4))
)

println("Custom Fermi surface model assembled")
println("Surface type: $(typeof(custom_surface))")
println("Vf(θ=0): $(custom_surface.fermi_velocity_func(0.0))")
println("Vf(θ=π/4): $(custom_surface.fermi_velocity_func(π/4))")

project_root = normpath(joinpath(@__DIR__, ".."))
mesh_path = joinpath(project_root, "assets", "square_bells.inp")
geometry_path = joinpath(project_root, "assets", "square_bells.geo")
boundary_conditions = Dict(
    :walls => MaxwellWallBC(1.0),
    :contact_top => OhmicContactBC(-0.5),
    :contact_bottom => OhmicContactBC(0.5),
)

problem = if isfile(mesh_path)
    TrixiProblem(; mesh_path=mesh_path, boundary_conditions=boundary_conditions)
elseif isfile(geometry_path)
    TrixiProblem(; geometry_path=geometry_path, boundary_conditions=boundary_conditions)
else
    error("No mesh found. Expected $(mesh_path) or $(geometry_path).")
end

config = SolverConfig(;
    polydeg = 3,
    cfl = 0.8,
    tspan_end = 50.0,
    residual_tol = 1e-5,
    log_every = 200,
)

callbacks = default_callbacks_builder()

sol, semi = solve(
    problem,
    model,
    config;
    callbacks=callbacks,
    visualize=live,
    visualization_mode=:mesh_native,
    name="example_custom_surface",
)
