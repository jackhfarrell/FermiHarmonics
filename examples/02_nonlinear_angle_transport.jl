# # Nonlinear Angle Grid Transport
#
# This example demonstrates nonlinear kinetic transport using an angle grid
# with parabolic band approximation.
#
# Angle grids are required for:
# - Strong (nonlinear) driving fields
# - Detailed angle-dependent physics
# - Problems where harmonic basis is insufficient
#
# Note: Nonlinear transport requires specifying band structure (mu0, mass)
# and cannot use the linear BGK collision models.

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

# For nonlinear transport, we typically use a parabolic band approximation
# Surface can be isotropic; the band details are in the collision model
surface = Isotropic2DFermiSurface(
    fermi_velocity = 1.0,           # Reference Fermi velocity
    nu = 1.0,           # density of states
    mass = 1.0,         # effective mass
    charge = -1.0       # electron charge
)

# Direct angle grid sampling: represents f(θ) at θ_n = 2πn/ntheta
# Resolution: ntheta angles sampled, use 32-64 for typical problems
discretization = AngleGrid(32)  # 32 equally-spaced angles

# Quadratic BGK collision: parabolic band + momentum-relaxing + e-e scattering
collision = QuadraticBGKCollision(
    0.05,                              # gamma_mr: momentum-relaxing rate
    OddQuarticRateProfile(0.3);        # gamma_ee profile (quartic m-dependence)
    mu0 = 0.5,                         # band bottom (chemical potential)
    mass = 1.0,                        # band mass
    electrostatic_coupling = 0.0,      # e-e interaction strength (optional)
    theta_oversample = 2               # Oversampling for FFT (optional)
)

# For angle grids, use IsotropicAngleStreaming (only option)
model = KineticModel2D(
    surface,
    discretization,
    IsotropicAngleStreaming(),     # Required for angle grids
    collision
)

println("Nonlinear angle grid model assembled successfully")
println("Surface: $(typeof(surface))")
println("Discretization: $(typeof(discretization))")
println("Collision: $(typeof(collision))")
println("Number of angles: 32")

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
    name="example_nonlinear_angle",
)
