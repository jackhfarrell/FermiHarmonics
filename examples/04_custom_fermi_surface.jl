# # Custom Fermi Surface
#
# This example demonstrates how to define a custom Fermi surface
# with user-supplied velocity function vF(θ).
#
# Custom surfaces are useful for:
# - Systems with known analytic formulas (elliptic, hexagonal, etc.)
# - Arbitrary band structures
# - Comparing theoretical vs. experimental dispersions

using ElectronKinetics

# Example 1: Elliptic Fermi surface (built-in analytic type)
# Anisotropy: vF(θ) = vF0 / √(cos²θ + sin²θ/a²)
# aspect ratio a > 1 means velocity higher in y-direction
elliptic_surface = EllipticFermiSurface2D(
    vF0 = 1.0,           # Fermi velocity (isotropic limit)
    aspect = 2.0,        # Anisotropy: v_y/v_x ratio
    nu = 1.0,            # density of states
    mass = 1.0,
    charge = -1.0
)

# Example 2: Custom user-defined surface with arbitrary function
# Must provide: vF_func, max_vF (upper bound for CFL), and properties
custom_surface = GeneralFermiSurface2D(
    # vF(θ) function: can be any formula
    # Here: modulated with 4-fold symmetry (square-like band)
    θ -> 1.0 * (1.0 + 0.3 * cos(4θ)),
    name = :square_modulated,
    max_vF = 1.3,        # Maximum of vF(θ) — MUST be tight upper bound
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
println("Vf(θ=0): $(custom_surface.vF_func(0.0))")
println("Vf(θ=π/4): $(custom_surface.vF_func(π/4))")
