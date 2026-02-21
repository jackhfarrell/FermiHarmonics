from pathlib import Path
import re

import h5py
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# Research-style script: edit constants here as needed.
PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS_ROOT = PROJECT_ROOT / "results"
SWEEP_GLOB = "simple_geometries_junction_sweep_*"
X_CUT_LEFT = -0.4
X_CUT_RIGHT = 0.4
BULK_Y_MIN = -0.2
BULK_Y_MAX = 0.2
GAMMA_MR_MIN = 1e-2
GAMMA_MR_MAX = 1e2
OUTFILE = Path(__file__).resolve().parent / "junction_bulk_delta_a0_left_minus_right.png"


sweep_dirs = sorted([p for p in RESULTS_ROOT.glob(SWEEP_GLOB) if p.is_dir()], key=lambda p: p.stat().st_mtime)
if not sweep_dirs:
    raise RuntimeError(f"No sweep directory matching {SWEEP_GLOB!r} found in {RESULTS_ROOT}")

sweep_dir = sweep_dirs[-1]
data_dir = sweep_dir / "data"
files = sorted(data_dir.glob("*.h5"))
if not files:
    raise RuntimeError(f"No .h5 files found in {data_dir}")

print(f"Using sweep directory: {sweep_dir}")
print(f"Found {len(files)} files")

pattern = re.compile(r"gamma_mc=([0-9eE+\-.]+)_gamma_mr=([0-9eE+\-.]+)")
rows = []

for fpath in files:
    m = pattern.search(fpath.name)
    if m is None:
        continue

    gamma_mc = float(m.group(1))
    gamma_mr = float(m.group(2))

    with h5py.File(fpath, "r") as h5:
        a0 = np.asarray(h5["a0"]).T
        x = np.asarray(h5["x"])
        y = np.asarray(h5["y"])
        mask = np.asarray(h5["mask"]).astype(bool).T

    ix_left = int(np.argmin(np.abs(x - X_CUT_LEFT)))
    ix_right = int(np.argmin(np.abs(x - X_CUT_RIGHT)))
    bulk_y = (y >= BULK_Y_MIN) & (y <= BULK_Y_MAX)

    a0_left_line = np.where(mask[ix_left, :] & bulk_y, a0[ix_left, :], np.nan)
    a0_right_line = np.where(mask[ix_right, :] & bulk_y, a0[ix_right, :], np.nan)

    a0_left_bulk = np.nanmean(a0_left_line)
    a0_right_bulk = np.nanmean(a0_right_line)
    delta_a0 = a0_left_bulk - a0_right_bulk

    if not np.isfinite(a0_left_bulk) or not np.isfinite(a0_right_bulk):
        delta_a0_mag = np.nan
    else:
        delta_a0_mag = np.abs(delta_a0)

    rows.append((gamma_mr, gamma_mc, delta_a0_mag))

if not rows:
    raise RuntimeError("No parseable files found (filename regex did not match).")

rows = [r for r in rows if GAMMA_MR_MIN <= r[0] <= GAMMA_MR_MAX]
if not rows:
    raise RuntimeError("No rows remain after gamma_mr truncation.")

gamma_mr_vals = np.array(sorted({r[0] for r in rows}))
gamma_mc_vals = np.array(sorted({r[1] for r in rows}))

delta_a0_mag_grid = np.full((len(gamma_mr_vals), len(gamma_mc_vals)), np.nan)

for gamma_mr, gamma_mc, delta_a0_mag in rows:
    i = np.searchsorted(gamma_mr_vals, gamma_mr)
    j = np.searchsorted(gamma_mc_vals, gamma_mc)
    delta_a0_mag_grid[i, j] = delta_a0_mag

fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)

finite_delta = delta_a0_mag_grid[np.isfinite(delta_a0_mag_grid)]
if finite_delta.size == 0:
    raise RuntimeError("No finite bulk |Δa0| values to plot.")

vmax = float(np.quantile(finite_delta, 0.99))
if not np.isfinite(vmax) or np.isclose(vmax, 0.0):
    vmax = float(np.nanmax(finite_delta))
    if not np.isfinite(vmax) or np.isclose(vmax, 0.0):
        vmax = 1.0

cmap = sns.color_palette("rocket", as_cmap=True).copy()
cmap.set_bad("0.85")
pcm = ax.pcolormesh(
    gamma_mr_vals,
    gamma_mc_vals,
    delta_a0_mag_grid.T,
    shading="auto",
    cmap=cmap,
    vmin=0.0,
    vmax=vmax,
)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(GAMMA_MR_MIN, GAMMA_MR_MAX)
ax.set_xlabel(r"$\gamma_{mr}$")
ax.set_ylabel(r"$\gamma_{mc}$")
ax.set_title(r"Bulk side potential difference magnitude: $|\Delta a_0|$")
cb = fig.colorbar(pcm, ax=ax)
cb.set_label(r"$|\Delta a_0|$")

fig.suptitle(
    "Junction bulk side-to-side potential difference magnitude\n"
    + rf"$a_{{0,\mathrm{{left}}}}=\langle a_0(x={X_CUT_LEFT:.2f},y)\rangle_{{{BULK_Y_MIN:.2f}\leq y\leq {BULK_Y_MAX:.2f}}}$, "
    + rf"$a_{{0,\mathrm{{right}}}}=\langle a_0(x={X_CUT_RIGHT:.2f},y)\rangle_{{{BULK_Y_MIN:.2f}\leq y\leq {BULK_Y_MAX:.2f}}}$, "
    + r"$|\Delta a_0| = |a_{0,\mathrm{left}} - a_{0,\mathrm{right}}|$"
)

fig.savefig(OUTFILE, dpi=220)
print(f"Saved figure: {OUTFILE}")
