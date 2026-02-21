from pathlib import Path
import re

import h5py
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
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
OUTFILE = Path(__file__).resolve().parent / "junction_terminal_ratios_ie_over_iw_in_over_ie_heatmap.png"
SMOOTHING_PASSES = 0


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


def smooth_nan_grid(grid: np.ndarray, passes: int = 1) -> np.ndarray:
    # Simple NaN-aware 3x3 box smoothing.
    k = np.ones((3, 3), dtype=float)
    out = grid.copy()
    for _ in range(max(0, passes)):
        padded = np.pad(out, 1, mode="edge")
        next_out = out.copy()
        for i in range(out.shape[0]):
            for j in range(out.shape[1]):
                w = padded[i:i + 3, j:j + 3]
                m = np.isfinite(w)
                if np.any(m):
                    next_out[i, j] = np.sum(w[m] * k[m]) / np.sum(k[m])
                else:
                    next_out[i, j] = np.nan
        out = next_out
    return out


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

    if not np.isfinite(i_e) or not np.isfinite(i_w) or np.isclose(i_w, 0.0):
        ie_over_iw = np.nan
    else:
        ie_over_iw = i_e / i_w

    if not np.isfinite(i_n) or not np.isfinite(i_e) or np.isclose(i_e, 0.0):
        in_over_ie = np.nan
    else:
        in_over_ie = i_n / i_e

    rows.append((gamma_mr, gamma_mc, i_e, i_n, i_w, ie_over_iw, in_over_ie))

if not rows:
    raise RuntimeError("No parseable files found (filename regex did not match).")

rows = [r for r in rows if GAMMA_MR_MIN <= r[0] <= GAMMA_MR_MAX and GAMMA_MC_MIN <= r[1] <= GAMMA_MC_MAX]
if not rows:
    raise RuntimeError("No rows remain after gamma_mr/gamma_mc truncation.")

gamma_mr_vals = np.array(sorted({r[0] for r in rows}))
gamma_mc_vals = np.array(sorted({r[1] for r in rows}))

ie_over_iw_grid = np.full((len(gamma_mr_vals), len(gamma_mc_vals)), np.nan)
in_over_ie_grid = np.full((len(gamma_mr_vals), len(gamma_mc_vals)), np.nan)

for gamma_mr, gamma_mc, _, _, _, ie_over_iw, in_over_ie in rows:
    i = np.searchsorted(gamma_mr_vals, gamma_mr)
    j = np.searchsorted(gamma_mc_vals, gamma_mc)
    ie_over_iw_grid[i, j] = ie_over_iw
    in_over_ie_grid[i, j] = in_over_ie

fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.6), constrained_layout=True, sharex=True, sharey=True)
plots = [
    (ie_over_iw_grid, r"$I_E / I_W$", "rocket"),
    (in_over_ie_grid, r"$I_N / I_E$", "mako"),
]

for ax, (grid, title, palette) in zip(axes, plots):
    grid_smooth = smooth_nan_grid(grid, passes=SMOOTHING_PASSES)
    finite = grid_smooth[np.isfinite(grid_smooth)]
    if finite.size == 0:
        raise RuntimeError(f"No finite values to plot for {title}.")
    vmin = float(np.quantile(finite, 0.01))
    vmax = float(np.quantile(finite, 0.99))
    if not np.isfinite(vmin) or not np.isfinite(vmax) or vmax <= vmin:
        vmin, vmax = float(np.nanmin(finite)), float(np.nanmax(finite))
        if not np.isfinite(vmin) or not np.isfinite(vmax) or vmax <= vmin:
            vmin, vmax = 0.0, 1.0

    cmap = sns.color_palette(palette, as_cmap=True).copy()
    cmap.set_bad("0.85")
    z = np.ma.masked_invalid(grid_smooth.T)
    fill_levels = np.linspace(vmin, vmax, 80)
    y_min = float(np.nanmin(gamma_mc_vals))
    y_max = float(np.nanmax(gamma_mc_vals))
    x_min = float(np.nanmin(gamma_mr_vals))
    x_max = float(np.nanmax(gamma_mr_vals))
    domain_aspect = (np.log10(y_max) - np.log10(y_min)) / (np.log10(x_max) - np.log10(x_min))
    pcm = ax.contourf(
        gamma_mr_vals,
        gamma_mc_vals,
        z,
        levels=fill_levels,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        extend="both",
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(GAMMA_MR_MIN, GAMMA_MR_MAX)
    ax.set_ylim(GAMMA_MC_MIN, GAMMA_MC_MAX)
    ax.set_box_aspect(domain_aspect)
    ax.set_xlabel(r"$\gamma_{mr}$")
    ax.set_title(title)
    cb = fig.colorbar(pcm, ax=ax)
    cb.formatter = mticker.FormatStrFormatter("%.2g")
    cb.update_ticks()
    cb.set_label(title)

axes[0].set_ylabel(r"$\gamma_{mc}$")

fig.savefig(OUTFILE, dpi=220)
print(f"Saved figure: {OUTFILE}")
