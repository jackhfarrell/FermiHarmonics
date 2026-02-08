using Documenter
using FermiHarmonics

# Build documentation
makedocs(;
    modules = [FermiHarmonics],
    sitename = "FermiHarmonics.jl",
    remotes = nothing,  # Disable remote source links
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://docs.jackhfarrell.com",
        assets = String[],
        mathengine = Documenter.MathJax3(),  # LaTeX rendering via MathJax3
    ),
    pages = [
        "Home" => "index.md",
        "Mesh" => "mesh.md",
        "API Reference" => [
            "Equations" => "api/equations.md",
            "Boundary Conditions" => "api/boundary_conditions.md",
            "Source Terms" => "api/source_terms.md",
            "I/O Utilities" => "api/io_utils.md",
            "Solve" => "api/solve.md",
            "SLURM Utilities" => "api/slurm_utils.md",
        ],
    ],
    checkdocs = :none,
)

deploydocs(;
    repo = "github.com/jackhfarrell/FermiHarmonics.git",
    devbranch = "main",
    cname = "docs.jackhfarrell.com",
)
