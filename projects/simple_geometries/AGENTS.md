# AGENTS.md

Minimal guidance for agents working in `projects/simple_geometries`.

## References

- Repository overview and setup: [README.md](/Users/jfarrell/Desktop/FermiHarmonics/README.md)
- Julia docs home: <https://docs.julialang.org/en/v1/>
- Julia manual: <https://docs.julialang.org/en/v1/manual/>
- Julia style guide: <https://docs.julialang.org/en/v1/manual/style-guide/>
- Julia performance tips: <https://docs.julialang.org/en/v1/manual/performance-tips/>
- Pkg docs: <https://pkgdocs.julialang.org/v1/>
- Julia testing docs: <https://docs.julialang.org/en/v1/stdlib/Test/>

## Project Hygiene

- Treat this as a Julia-centered subproject inside the main package repository.
- Run Julia commands from the repository root with the project environment active, typically `julia --project=.`.
- Keep reusable Julia code in [src/](/Users/jfarrell/Desktop/FermiHarmonics/src), geometry-specific scripts in [projects/simple_geometries/scripts](/Users/jfarrell/Desktop/FermiHarmonics/projects/simple_geometries/scripts), and analysis artifacts in sibling subfolders only when they are specific to this subproject.
- Keep dependencies declared in [Project.toml](/Users/jfarrell/Desktop/FermiHarmonics/Project.toml) and avoid adding unnecessary packages for one-off scripts.
- When behavior changes in shared solver logic, add or update tests under [test/](/Users/jfarrell/Desktop/FermiHarmonics/test).
- Prefer clear, type-stable, allocation-aware Julia code and use the official Julia docs above as the first reference for language, package, testing, and performance questions.
