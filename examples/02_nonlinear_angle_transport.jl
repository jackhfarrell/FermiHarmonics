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

# For nonlinear transport, we typically use a parabolic band approximation
# Surface can be isotropic; the band details are in the collision model
surface = Isotropic2DFermiSurface(
    vF = 1.0,           # Reference Fermi velocity
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
