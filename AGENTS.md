# AGENTS.md

Minimal guidance for agents working in this Julia repository.

## References

- Project overview and setup: [README.md](/Users/jfarrell/Desktop/ElectronKinetics.jl/README.md)
- Project documentation: <https://fermiharmonics.jackhfarrell.com>
- Julia docs home: <https://docs.julialang.org/en/v1/>
- Julia manual: <https://docs.julialang.org/en/v1/manual/>
- Julia style guide: <https://docs.julialang.org/en/v1/manual/style-guide/>
- Julia performance tips: <https://docs.julialang.org/en/v1/manual/performance-tips/>
- Makie docs: <https://docs.makie.org/stable/>
- Pkg docs: <https://pkgdocs.julialang.org/v1/>
- Trixi docs: <https://trixi-framework.github.io/TrixiDocumentation/stable/>
- OrdinaryDiffEq docs: <https://docs.sciml.ai/OrdinaryDiffEq/stable/>
- Julia testing docs: <https://docs.julialang.org/en/v1/stdlib/Test/>

## Project Hygiene

- Treat this as a standard Julia package with `Project.toml`, `Manifest.toml`, `src/`, `test/`, and `docs/`.
- Run Julia commands with the project environment active, typically `julia --project=.`.
- Prefer package code in `src/`, tests in `test/`, and small runnable examples in `demo/` or `projects/`.
- Keep dependencies declared in [Project.toml](/Users/jfarrell/Desktop/ElectronKinetics.jl/Project.toml) and avoid introducing unnecessary new packages.
- When changing behavior, add or update tests under [test/](/Users/jfarrell/Desktop/ElectronKinetics.jl/test).
- Prefer clear, type-stable, allocation-aware Julia code and use the official Julia docs above as the first reference for language, package, testing, and performance questions.
