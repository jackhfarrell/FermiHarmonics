# square_bells_ucsb

This folder contains scripts, meshes, and notebooks for square-bells UCSB parameter studies.

## What is tracked

- `scripts/`: run and sweep submission scripts.
- `mesh/`: mesh and geometry inputs.
- `analyze_sweep.ipynb`: sweep-level analysis notebook (cleaned for publication).

## What is not tracked

Generated results are intentionally excluded from git, including:

- `results/`
- large generated `.h5` data files

## Cluster sweep workflow (SLURM)

Before submitting, configure your SLURM account/partition/qos in:

- `projects/square_bells_ucsb/scripts/square_bells_array.jl`

Relevant block:

- `sbatch = Dict(...)`
  - `:account`
  - `:partition`
  - `:qos`
  - and any other site-specific settings you need.

Then submit both wall-scattering sweeps:

```bash
julia --project=. projects/square_bells_ucsb/scripts/submit_two_p_sweeps.jl
```

This submits:

- one sweep at `p_scatter=1.0` with 250 jobs
- one sweep at `p_scatter=0.0` with 250 jobs
- total: **500 SLURM array jobs**

Each sweep writes results under:

- `projects/square_bells_ucsb/results/square_bells_ucsb_sweep_p=<...>_<timestamp>/data/*.h5`

## Analyze a completed sweep

1. Open `projects/square_bells_ucsb/analyze_sweep.ipynb`.
2. Set `SWEEP_DATA_DIR` in the first configuration cell to your `.../data` directory.
3. Choose `(I_IDX, J_IDX)` and plotting options.
4. Run the notebook top-to-bottom.
