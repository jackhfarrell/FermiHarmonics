# # Linear Harmonic Basis Transport
#
# This example demonstrates the most common use case: linear Fermi-liquid
# transport using harmonic basis expansion with automatic mode estimation.
#
# The harmonic basis is memory-efficient and suitable for:
# - Linear response regime (small perturbations)
# - Systems with well-separated timescales (weak damping)
# - Routine kinetic transport calculations

using ElectronKinetics

# Define carrier properties (e.g., graphene electrons)
# Surface specifies the Fermi surface shape and quasiparticle properties
surface = Isotropic2DFermiSurface(
    vF = 1.0,           # Fermi velocity (energy units)
    nu = 1.0,           # density of states
    mass = 1.0,         # effective mass
    charge = -1.0       # electron charge (negative for electrons)
)

# Choose angular discretization with automatic mode estimation
# :auto estimates the maximum harmonic from collision parameters
# Alternative: HarmonicBasis(10) for fixed 10 modes
discretization = HarmonicBasis(:auto)

# Define collision model using linear BGK approximation
# Parameters:
#   gamma_mr = 0.1  — momentum-relaxing (drag) scattering rate
#   TwoRateProfile(0.4) — mode-dependent rates: γ=0 for m<2, γ=0.4 for m≥2
collision = LinearBGKCollision(
    0.1,                           # momentum-relaxing rate
    TwoRateProfile(0.4)            # mode-dependent rate profile
)

# Assemble the kinetic model by combining components
model = KineticModel2D(
    surface,
    discretization,
    IsotropicHarmonicStreaming(),  # Only choice for harmonic basis
    collision
)

# The model is now ready to use with a solver backend (e.g., Trixi)
# See the TrixiProblem documentation for mesh-based solving
println("Linear harmonic model assembled successfully")
println("Surface: $(typeof(surface))")
println("Discretization: $(typeof(discretization))")
println("Collision: $(typeof(collision))")
