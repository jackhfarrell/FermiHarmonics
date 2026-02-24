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
X_CUT_WEST = -0.4
X_CUT_EAST = 0.4
Y_CUT_NORTH = 0.4
GAMMA_MR_MIN = 1e-2
GAMMA_MR_MAX = 1e2
GAMMA_MC_MIN = 1e-2
GAMMA_MC_MAX = 1e2
DELTA_V = 1.0
OUTFILE = Path(__file__).resolve().parent / "junction_total_conductance_heatmap.png"


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

    ix_west = int(np.argmin(np.abs(x - X_CUT_WEST)))
    ix_east = int(np.argmin(np.abs(x - X_CUT_EAST)))
    iy_north = int(np.argmin(np.abs(y - Y_CUT_NORTH)))

    a1_west_line = np.where(mask[ix_west, :], a1[ix_west, :], np.nan)
    valid_y_west = np.isfinite(a1_west_line)
    if np.count_nonzero(valid_y_west) < 2:
        i_w = np.nan
    else:
        i_w = np.abs(integrate_over_valid_segments(y, a1_west_line, valid_y_west))

    a1_east_line = np.where(mask[ix_east, :], a1[ix_east, :], np.nan)
    valid_y_east = np.isfinite(a1_east_line)
    if np.count_nonzero(valid_y_east) < 2:
        i_e = np.nan
    else:
        i_e = np.abs(integrate_over_valid_segments(y, a1_east_line, valid_y_east))

    b1_north_line = np.where(mask[:, iy_north], b1[:, iy_north], np.nan)
    valid_x_north = np.isfinite(b1_north_line)
    if np.count_nonzero(valid_x_north) < 2:
        i_n = np.nan
    else:
        i_n = np.abs(integrate_over_valid_segments(x, b1_north_line, valid_x_north))

    if not np.isfinite(i_w) or not np.isfinite(i_e) or not np.isfinite(i_n) or np.isclose(DELTA_V, 0.0):
        g_total = np.nan
    else:
        i_total = i_w + i_e + i_n
        g_total = i_total / DELTA_V

    rows.append((gamma_mr, gamma_mc, g_total))

if not rows:
    raise RuntimeError("No parseable files found (filename regex did not match).")

rows = [r for r in rows if GAMMA_MR_MIN <= r[0] <= GAMMA_MR_MAX and GAMMA_MC_MIN <= r[1] <= GAMMA_MC_MAX]
if not rows:
    raise RuntimeError("No rows remain after gamma_mr/gamma_mc truncation.")

gamma_mr_vals = np.array(sorted({r[0] for r in rows}))
gamma_mc_vals = np.array(sorted({r[1] for r in rows}))

g_total_grid = np.full((len(gamma_mr_vals), len(gamma_mc_vals)), np.nan)

for gamma_mr, gamma_mc, g_total in rows:
    i = np.searchsorted(gamma_mr_vals, gamma_mr)
    j = np.searchsorted(gamma_mc_vals, gamma_mc)
    g_total_grid[i, j] = g_total

fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)
finite_g = g_total_grid[np.isfinite(g_total_grid)]
if finite_g.size == 0:
    raise RuntimeError("No finite total conductance values to plot.")

#vmin = float(np.quantile(finite_g, 0.01))
#vmax = float(np.quantile(finite_g, 0.99))
#if not np.isfinite(vmin) or not np.isfinite(vmax) or vmax <= vmin:
#    vmin, vmax = float(np.nanmin(finite_g)), float(np.nanmax(finite_g))
#    if not np.isfinite(vmin) or not np.isfinite(vmax) or vmax <= vmin:
#        vmin, vmax = 0.0, 1.0

cmap = sns.color_palette("crest", as_cmap=True).copy()
cmap.set_bad("0.85")
pcm = ax.pcolormesh(gamma_mr_vals, gamma_mc_vals, g_total_grid.T, shading="auto", cmap=cmap)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(GAMMA_MR_MIN, GAMMA_MR_MAX)
ax.set_ylim(GAMMA_MC_MIN, GAMMA_MC_MAX)
ax.set_xlabel(r"$\gamma_{mr}$")
ax.set_ylabel(r"$\gamma_{mc}$")
ax.set_title(
    r"Junction total conductance heatmap"
    + "\n"
    + rf"$G_{{\mathrm{{tot}}}}=(I_W+I_E+I_N)/\Delta V$, $\Delta V={DELTA_V:.2f}$"
)
cb = fig.colorbar(pcm, ax=ax)
cb.set_label(r"$G_{\mathrm{tot}}$")

fig.savefig(OUTFILE, dpi=220)
print(f"Saved figure: {OUTFILE}")
