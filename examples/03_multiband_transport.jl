# # Multiband Linear Transport
#
# This example demonstrates how to model multiple carrier species
# (e.g., electrons and holes) with coupled kinetic transport.
#
# Multiband support:
# - Uses linear harmonic basis (not angle grid)
# - Each band has its own surface (Fermi velocity, mass, etc.)
# - Bands are coupled through inter-band drag (optional)
# - Useful for compensated semimetals, graphene bilayers, etc.

using ElectronKinetics

# Define two carrier bands (electrons and holes)
electron_band = Band(
    :electrons,                          # band name
    Isotropic2DFermiSurface(
        vF = 1.0,
        nu = 0.5,
        mass = 0.5,
        charge = -1.0
    ),
    gamma_mr = 0.1,                      # momentum-relaxing rate
    gamma_ee = 0.05                      # e-e scattering rate
)

hole_band = Band(
    :holes,
    Isotropic2DFermiSurface(
        vF = 0.8,
        nu = 0.5,
        mass = 0.6,
        charge = +1.0
    ),
    gamma_mr = 0.08,
    gamma_ee = 0.04
)

# Create model with multiple bands and optional drag coupling
# Bands are passed as vector; first band determines surface for streaming
# gamma_drag: inter-band momentum transfer rate (0.0 = no coupling)
model = KineticModel2D(
    electron_band.surface,               # Base surface (for streaming matrices)
    HarmonicBasis(:auto),
    IsotropicHarmonicStreaming(),
    LinearBGKCollision(0.1),             # Collision model for main band
    bands = [electron_band, hole_band],  # All bands as vector
    gamma_drag = 0.02                    # Inter-band drag coupling strength
)

println("Multiband linear harmonic model assembled successfully")
println("Number of bands: 2")
println("Band names: $(band.name for band in [electron_band, hole_band])")
println("Drag coupling: gamma_drag = 0.02")
