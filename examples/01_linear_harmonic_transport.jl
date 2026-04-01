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
using Trixi

live = "--live" in ARGS
mesh_only = "--mesh-only" in ARGS
if live
    try
        @eval using GLMakie
    catch err
        error("Live visualization requested, but GLMakie is not available: $(err)")
    end
end

# Define carrier properties (e.g., graphene electrons)
# Surface specifies the Fermi surface shape and quasiparticle properties
surface = Isotropic2DFermiSurface(
    fermi_velocity = 1.0,           # Fermi velocity (energy units)
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

project_root = normpath(joinpath(@__DIR__, ".."))
mesh_path = joinpath(project_root, "assets", "square_bells.inp")
geometry_path = joinpath(project_root, "assets", "square_bells.geo")
boundary_conditions = Dict(
    :walls => MaxwellWallBC(1.0),
    :contact_top => OhmicContactBC(-0.5),
    :contact_bottom => OhmicContactBC(0.5),
)

config = SolverConfig(;
    polydeg = 3,
    cfl = 0.8,
    tspan_end = 50.0,
    residual_tol = 1e-5,
    log_every = 200,
)

problem = if isfile(mesh_path)
    TrixiProblem(; mesh_path=mesh_path, boundary_conditions=boundary_conditions)
elseif isfile(geometry_path)
    TrixiProblem(;
        geometry_path=geometry_path,
        boundary_conditions=boundary_conditions,
        mesh_build=MeshBuildConfig(output_mode=:persistent, mesh_scale=3.0),
    )
else
    error("No mesh found. Expected $(mesh_path) or $(geometry_path).")
end

if live
    preview_mesh(
        problem,
        model,
        config;
        visualization_mode=:mesh_native,
        mesh_scale=3.0,
        wait_for_close=true,
        name="example_linear_harmonic_mesh",
    )
end

if mesh_only
    exit()
end

callbacks = default_callbacks_builder()

sol, semi = solve(
    problem,
    model,
    config;
    callbacks=callbacks,
    visualize=live,
    visualization_mode=:mesh_native,
    name="example_linear_harmonic",
)
