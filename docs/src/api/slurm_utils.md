# SLURM Utilities API

```@docs
FermiHarmonics.submit_sweep!
FermiHarmonics.write_sweep_metadata!
FermiHarmonics.archive_mesh!
FermiHarmonics.copy_mesh_to_scratch
FermiHarmonics.select_cases
FermiHarmonics.grid_lookup
```

## Behavior Notes

- `submit_sweep!` returns a `Cmd` only in `dry_run=true`; otherwise it submits and returns `nothing`.
- `write_sweep_metadata!` writes/overwrites `sweep_metadata.toml`.
- `select_cases` uses `SLURM_ARRAY_TASK_ID` from environment and defaults to task `1` if unset.
