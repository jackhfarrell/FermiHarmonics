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
Y_CUT = 0.4
X_CUT_RIGHT = 0.4
GAMMA_MR_MIN = 1e-2
GAMMA_MR_MAX = 1e2
OUTFILE = Path(__file__).resolve().parent / "junction_i2_over_i3_heatmap.png"


def integrate_over_valid_segments(coord: np.ndarray, values: np.ndarray, valid: np.ndarray) -> float:
    total = 0.0
    n = len(coord)
    i = 0
    while i < n:
        if not valid[i]:
            i += 1
            continue
        j = i
        while j < n and valid[j]:
            j += 1
        if j - i >= 2:
            total += np.trapezoid(values[i:j], coord[i:j])
        i = j
    return total


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
        a1 = np.asarray(h5["a1"]).T
        b1 = np.asarray(h5["b1"]).T
        x = np.asarray(h5["x"])
        y = np.asarray(h5["y"])
        mask = np.asarray(h5["mask"]).astype(bool).T

    iy = int(np.argmin(np.abs(y - Y_CUT)))
    ix_right = int(np.argmin(np.abs(x - X_CUT_RIGHT)))

    b1_line = np.where(mask[:, iy], b1[:, iy], np.nan)
    valid_x = np.isfinite(b1_line)
    if np.count_nonzero(valid_x) < 2:
        i2 = np.nan
    else:
        i2 = np.abs(integrate_over_valid_segments(x, b1_line, valid_x))

    a1_right_line = np.where(mask[ix_right, :], a1[ix_right, :], np.nan)
    valid_y_right = np.isfinite(a1_right_line)
    if np.count_nonzero(valid_y_right) < 2:
        i3 = np.nan
    else:
        i3 = np.abs(integrate_over_valid_segments(y, a1_right_line, valid_y_right))

    if not np.isfinite(i2) or not np.isfinite(i3) or np.isclose(i3, 0.0):
        ratio = np.nan
    else:
        ratio = i2 / i3

    rows.append((gamma_mr, gamma_mc, i2, i3, ratio))

if not rows:
    raise RuntimeError("No parseable files found (filename regex did not match).")

rows = [r for r in rows if GAMMA_MR_MIN <= r[0] <= GAMMA_MR_MAX]
if not rows:
    raise RuntimeError("No rows remain after gamma_mr truncation.")

gamma_mr_vals = np.array(sorted({r[0] for r in rows}))
gamma_mc_vals = np.array(sorted({r[1] for r in rows}))

ratio_grid = np.full((len(gamma_mr_vals), len(gamma_mc_vals)), np.nan)
for gamma_mr, gamma_mc, _, _, ratio in rows:
    i = np.searchsorted(gamma_mr_vals, gamma_mr)
    j = np.searchsorted(gamma_mc_vals, gamma_mc)
    ratio_grid[i, j] = ratio

fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)
finite_ratio = ratio_grid[np.isfinite(ratio_grid)]
if finite_ratio.size == 0:
    raise RuntimeError("No finite I2/I3 values to plot.")

vmax = float(np.quantile(finite_ratio, 0.99))
vmin = float(np.quantile(finite_ratio, 0.01))
if not np.isfinite(vmin) or not np.isfinite(vmax) or vmax <= vmin:
    vmin, vmax = float(np.nanmin(finite_ratio)), float(np.nanmax(finite_ratio))
    if not np.isfinite(vmin) or not np.isfinite(vmax) or vmax <= vmin:
        vmin, vmax = 0.0, 1.0

cmap = sns.color_palette("rocket", as_cmap=True).copy()
cmap.set_bad("0.85")

pcm = ax.pcolormesh(gamma_mr_vals, gamma_mc_vals, ratio_grid.T, shading="auto", cmap=cmap, vmin=vmin, vmax=vmax)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(GAMMA_MR_MIN, GAMMA_MR_MAX)
ax.set_xlabel(r"$\gamma_{mr}$")
ax.set_ylabel(r"$\gamma_{mc}$")
ax.set_title(
    r"Heatmap of $I_2 / I_3$"
    + "\n"
    + rf"$I_2=|\int b_1(x,y={Y_CUT:.2f})\,dx|$, $I_3=|\int a_1(x={X_CUT_RIGHT:.2f},y)\,dy|$, $10^{{-2}}\leq\gamma_{{mr}}\leq10^2$"
)
cb = fig.colorbar(pcm, ax=ax)
cb.set_label(r"$I_2 / I_3$")

fig.savefig(OUTFILE, dpi=220)
print(f"Saved figure: {OUTFILE}")
