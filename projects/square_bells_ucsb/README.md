# square_bells_ucsb

This folder contains scripts, meshes, and notebooks for square-bells UCSB parameter studies.

## What is tracked

- `scripts/`: run and sweep submission scripts.
- `mesh/`: mesh and geometry inputs.
- `analyze_sweep.ipynb`: sweep-level analysis notebook (cleaned for publication).
- `inspect_observables.ipynb`: single-run observable inspection notebook (cleaned for publication).

## What is not tracked

Generated results are intentionally excluded from git:

- `results/`
- large generated `.h5` data files

## Expected local data layout

The notebooks are configured to look for repo-local generated data under:

- `results/example_sweep/data/*.h5`
- `results/example_single_run/data/*.h5`

You can change these paths in the first data-loading cells.

## Typical workflow

1. Run a simulation or sweep via `scripts/`.
2. Save analysis outputs under `results/<run_or_sweep>/data/`.
3. Open notebooks and point to the relevant `results/.../data` path.
