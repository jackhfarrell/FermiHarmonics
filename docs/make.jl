using Documenter
using FermiFlows
using Trixi

# Build documentation
makedocs(;
    modules = [FermiFlows],
    sitename = "FermiFlows.jl",
    remotes = nothing,  # Disable remote source links
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://fermiflows.jackhfarrell.com",
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

if get(ENV, "CI", "false") == "true"
    deploydocs(;
        repo = "github.com/jackhfarrell/FermiHarmonics.git",
        devbranch = "main",
        cname = "fermiflows.jackhfarrell.com",
    )
end
