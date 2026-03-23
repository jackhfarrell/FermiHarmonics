# T = 0 Fermi liquid with circular Fermi Surface in Trixi.jl

"""
    FermiHarmonics

FermiHarmonics solves the 2D linearized Boltzmann system in harmonic variables
on a fixed mesh (no adaptive mesh refinement):

``\\partial_t \\mathbf{u} + A_x\\,\\partial_x \\mathbf{u} + A_y\\,\\partial_y \\mathbf{u}
= \\mathbf{S}_{\\mathrm{phys}}(\\mathbf{u}).``

State layout for maximum harmonic ``M``:
``\\mathbf{u} = (a_0, a_1, b_1, \\ldots, a_M, b_M)^T``, with ``n_{\\mathrm{vars}} = 1 + 2M``.

Physical scattering source:
- ``a_0`` undamped,
- ``a_1, b_1`` damped at ``\\gamma_{mr}``,
- higher harmonics damped at ``\\gamma_{mr} + \\gamma_{mc}``.
"""
module FermiHarmonics
using SciMLBase
using Trixi
using LinearAlgebra
using SparseArrays
using StaticArrays
using HDF5
using FFTW
using NLsolve
using Plots

# Include core physics and interfaces
include("equations.jl")
include("nonlinear_transport.jl")
include("reference_setup.jl")
include("source_terms.jl")
include("boundary_conditions.jl")
include("trixi_interface.jl")
include("slurm_utils.jl")

include("runner.jl")


# Simulation infrastructure
include("io_utils.jl")

# Export main types and functions
export FermiHarmonics2D,
       FermiAngles2D,
       SolveParams,
       MaxwellWallBC,
       OhmicContactBC,
       blg_reference_setup,
       estimate_max_harmonic,
       save_solution_custom,
       save_for_analysis,
       evaluate_solution,
       evaluate_observables,
       solve_status,
       solve,
       write_sweep_metadata!,
       archive_mesh!,
       copy_mesh_to_scratch,
       select_cases,
       grid_lookup,
       submit_sweep!

end # module FermiHarmonics
