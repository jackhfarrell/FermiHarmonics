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
using Trixi

live = "--live" in ARGS
if live
    try
        @eval using GLMakie
    catch err
        error("Live visualization requested, but GLMakie is not available: $(err)")
    end
end

# Define two carrier bands (electrons and holes)
electron_band = Band(
    :electrons,                          # band name
    Isotropic2DFermiSurface(
        fermi_velocity = 1.0,
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
        fermi_velocity = 0.8,
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
    name="example_multiband",
)
