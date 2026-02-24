from pathlib import Path
import re

import h5py
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.colors import ListedColormap

# Uses the same definitions/cuts as plot_junction_terminal_fractions_esw_heatmap.py
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

OUT_FIG = Path(__file__).resolve().parent / "junction_fraction_jacobian_diagnostics.png"
OUT_CSV = Path(__file__).resolve().parent / "junction_fraction_jacobian_summary.csv"


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


def finite_differences_log(grid: np.ndarray, logx: np.ndarray, logy: np.ndarray):
    """Return partials dgrid/dlogx and dgrid/dlogy with central differences."""
    nx, ny = grid.shape
    dlogx = np.full_like(grid, np.nan, dtype=float)
    dlogy = np.full_like(grid, np.nan, dtype=float)

    for i in range(nx):
        for j in range(ny):
            if not np.isfinite(grid[i, j]):
                continue

            # d/dlogx (gamma_mr axis)
            if i > 0 and i < nx - 1:
                f0, f1, f2 = grid[i - 1, j], grid[i, j], grid[i + 1, j]
                if np.isfinite(f0) and np.isfinite(f1) and np.isfinite(f2):
                    dlogx[i, j] = (f2 - f0) / (logx[i + 1] - logx[i - 1])
            elif i == 0 and nx > 1:
                f1, f2 = grid[i, j], grid[i + 1, j]
                if np.isfinite(f1) and np.isfinite(f2):
                    dlogx[i, j] = (f2 - f1) / (logx[i + 1] - logx[i])
            elif i == nx - 1 and nx > 1:
                f0, f1 = grid[i - 1, j], grid[i, j]
                if np.isfinite(f0) and np.isfinite(f1):
                    dlogx[i, j] = (f1 - f0) / (logx[i] - logx[i - 1])

            # d/dlogy (gamma_mc axis)
            if j > 0 and j < ny - 1:
                f0, f1, f2 = grid[i, j - 1], grid[i, j], grid[i, j + 1]
                if np.isfinite(f0) and np.isfinite(f1) and np.isfinite(f2):
                    dlogy[i, j] = (f2 - f0) / (logy[j + 1] - logy[j - 1])
            elif j == 0 and ny > 1:
                f1, f2 = grid[i, j], grid[i, j + 1]
                if np.isfinite(f1) and np.isfinite(f2):
                    dlogy[i, j] = (f2 - f1) / (logy[j + 1] - logy[j])
            elif j == ny - 1 and ny > 1:
                f0, f1 = grid[i, j - 1], grid[i, j]
                if np.isfinite(f0) and np.isfinite(f1):
                    dlogy[i, j] = (f1 - f0) / (logy[j] - logy[j - 1])

    return dlogx, dlogy


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
    i_w = np.nan if np.count_nonzero(valid_y_west) < 2 else np.abs(integrate_over_valid_segments(y, a1_west_line, valid_y_west))

    a1_east_line = np.where(mask[ix_east, :], a1[ix_east, :], np.nan)
    valid_y_east = np.isfinite(a1_east_line)
    i_e = np.nan if np.count_nonzero(valid_y_east) < 2 else np.abs(integrate_over_valid_segments(y, a1_east_line, valid_y_east))

    b1_north_line = np.where(mask[:, iy_north], b1[:, iy_north], np.nan)
    valid_x_north = np.isfinite(b1_north_line)
    i_n = np.nan if np.count_nonzero(valid_x_north) < 2 else np.abs(integrate_over_valid_segments(x, b1_north_line, valid_x_north))

    ie_over_iw = np.nan if (not np.isfinite(i_e) or not np.isfinite(i_w) or np.isclose(i_w, 0.0)) else i_e / i_w
    in_over_ie = np.nan if (not np.isfinite(i_n) or not np.isfinite(i_e) or np.isclose(i_e, 0.0)) else i_n / i_e
    rows.append((gamma_mr, gamma_mc, ie_over_iw, in_over_ie))

if not rows:
    raise RuntimeError("No parseable files found (filename regex did not match).")

rows = [r for r in rows if GAMMA_MR_MIN <= r[0] <= GAMMA_MR_MAX and GAMMA_MC_MIN <= r[1] <= GAMMA_MC_MAX]
if not rows:
    raise RuntimeError("No rows remain after gamma_mr/gamma_mc truncation.")

gamma_mr_vals = np.array(sorted({r[0] for r in rows}))
gamma_mc_vals = np.array(sorted({r[1] for r in rows}))
log_gamma_mr = np.log10(gamma_mr_vals)
log_gamma_mc = np.log10(gamma_mc_vals)

ie_over_iw_grid = np.full((len(gamma_mr_vals), len(gamma_mc_vals)), np.nan)
in_over_ie_grid = np.full((len(gamma_mr_vals), len(gamma_mc_vals)), np.nan)
for gamma_mr, gamma_mc, ie_over_iw, in_over_ie in rows:
    i = np.searchsorted(gamma_mr_vals, gamma_mr)
    j = np.searchsorted(gamma_mc_vals, gamma_mc)
    ie_over_iw_grid[i, j] = ie_over_iw
    in_over_ie_grid[i, j] = in_over_ie

# J = [[d(IE/IW)/dlog(gamma_mr), d(IE/IW)/dlog(gamma_mc)],
#      [d(IN/IE)/dlog(gamma_mr), d(IN/IE)/dlog(gamma_mc)]]
d_ie_iw_dlog_mr, d_ie_iw_dlog_mc = finite_differences_log(ie_over_iw_grid, log_gamma_mr, log_gamma_mc)
d_in_ie_dlog_mr, d_in_ie_dlog_mc = finite_differences_log(in_over_ie_grid, log_gamma_mr, log_gamma_mc)

det_j = d_ie_iw_dlog_mr * d_in_ie_dlog_mc - d_ie_iw_dlog_mc * d_in_ie_dlog_mr

# 2x2 condition number via singular values.
cond_j = np.full_like(det_j, np.nan, dtype=float)
for i in range(det_j.shape[0]):
    for j in range(det_j.shape[1]):
        J = np.array(
            [
                [d_ie_iw_dlog_mr[i, j], d_ie_iw_dlog_mc[i, j]],
                [d_in_ie_dlog_mr[i, j], d_in_ie_dlog_mc[i, j]],
            ],
            dtype=float,
        )
        if not np.all(np.isfinite(J)):
            continue
        svals = np.linalg.svd(J, compute_uv=False)
        if svals[-1] <= 0.0:
            cond_j[i, j] = np.inf
        else:
            cond_j[i, j] = svals[0] / svals[-1]

finite_det = det_j[np.isfinite(det_j)]
finite_cond = cond_j[np.isfinite(cond_j)]
if finite_det.size == 0 or finite_cond.size == 0:
    raise RuntimeError("No finite Jacobians to summarize.")

abs_det = np.abs(det_j)
eps_det = 1e-6
good_invertibility = np.isfinite(abs_det) & (abs_det > eps_det)
good_condition = np.isfinite(cond_j) & (cond_j < 10.0)
good_both = good_invertibility & good_condition
finite_mask = np.isfinite(det_j) & np.isfinite(cond_j)

print(f"Finite Jacobian cells: {np.count_nonzero(finite_mask)} / {finite_mask.size}")
print(f"|det(J)| > {eps_det:g}: {np.count_nonzero(good_invertibility)} cells")
print(f"cond(J) < 10: {np.count_nonzero(good_condition)} cells")
print(f"Both criteria: {np.count_nonzero(good_both)} cells")

print("det(J) summary:")
print(f"  min={np.nanmin(finite_det):.6g}, max={np.nanmax(finite_det):.6g}, median={np.nanmedian(finite_det):.6g}")
print("|det(J)| summary:")
print(f"  p10={np.nanpercentile(np.abs(finite_det), 10):.6g}, p50={np.nanpercentile(np.abs(finite_det), 50):.6g}, p90={np.nanpercentile(np.abs(finite_det), 90):.6g}")
print("cond(J) summary:")
print(f"  p10={np.nanpercentile(finite_cond, 10):.6g}, p50={np.nanpercentile(finite_cond, 50):.6g}, p90={np.nanpercentile(finite_cond, 90):.6g}, max={np.nanmax(finite_cond):.6g}")

# Save CSV summary for downstream use.
with open(OUT_CSV, "w", encoding="utf-8") as f:
    f.write("gamma_mr,gamma_mc,ie_over_iw,in_over_ie,d_ie_iw_dlog_mr,d_ie_iw_dlog_mc,d_in_ie_dlog_mr,d_in_ie_dlog_mc,det_j,abs_det_j,cond_j,good_invertibility,good_condition,good_both\n")
    for i, gamma_mr in enumerate(gamma_mr_vals):
        for j, gamma_mc in enumerate(gamma_mc_vals):
            f.write(
                f"{gamma_mr},{gamma_mc},{ie_over_iw_grid[i,j]},{in_over_ie_grid[i,j]},{d_ie_iw_dlog_mr[i,j]},"
                f"{d_ie_iw_dlog_mc[i,j]},{d_in_ie_dlog_mr[i,j]},{d_in_ie_dlog_mc[i,j]},{det_j[i,j]},"
                f"{abs_det[i,j]},{cond_j[i,j]},{int(good_invertibility[i,j])},{int(good_condition[i,j])},{int(good_both[i,j])}\n"
            )
print(f"Saved CSV: {OUT_CSV}")

fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.8), constrained_layout=True, sharex=True, sharey=True)

z1 = np.log10(np.abs(det_j))
cmap1 = sns.color_palette("mako", as_cmap=True).copy()
cmap1.set_bad("0.85")
v1 = z1[np.isfinite(z1)]
vmin1, vmax1 = (np.nanpercentile(v1, 2), np.nanpercentile(v1, 98)) if v1.size else (-6, 0)
pcm1 = axes[0].pcolormesh(gamma_mr_vals, gamma_mc_vals, z1.T, shading="auto", cmap=cmap1, vmin=vmin1, vmax=vmax1)
axes[0].set_title(r"$\log_{10}| \det J |$")
cb1 = fig.colorbar(pcm1, ax=axes[0])
cb1.set_label(r"$\log_{10}| \det J |$")

z2 = np.log10(cond_j)
cmap2 = sns.color_palette("rocket", as_cmap=True).copy()
cmap2.set_bad("0.85")
v2 = z2[np.isfinite(z2)]
vmin2, vmax2 = (np.nanpercentile(v2, 2), np.nanpercentile(v2, 98)) if v2.size else (0, 4)
pcm2 = axes[1].pcolormesh(gamma_mr_vals, gamma_mc_vals, z2.T, shading="auto", cmap=cmap2, vmin=vmin2, vmax=vmax2)
axes[1].set_title(r"$\log_{10}\kappa(J)$")
cb2 = fig.colorbar(pcm2, ax=axes[1])
cb2.set_label(r"$\log_{10}\kappa(J)$")

mask_int = np.where(good_both, 1.0, np.where(finite_mask, 0.0, np.nan))
cmap3 = ListedColormap(["#D9D9D9", "#1b9e77"])
cmap3.set_bad("0.85")
pcm3 = axes[2].pcolormesh(gamma_mr_vals, gamma_mc_vals, mask_int.T, shading="auto", cmap=cmap3, vmin=0.0, vmax=1.0)
axes[2].set_title(r"Good Region ($|\det J|>10^{-6}$, $\kappa<10$)")
cb3 = fig.colorbar(pcm3, ax=axes[2], ticks=[0, 1])
cb3.ax.set_yticklabels(["no", "yes"])

for ax in axes:
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(GAMMA_MR_MIN, GAMMA_MR_MAX)
    ax.set_ylim(GAMMA_MC_MIN, GAMMA_MC_MAX)
    ax.set_xlabel(r"$\gamma_{mr}$")
axes[0].set_ylabel(r"$\gamma_{mc}$")

fig.suptitle(r"Jacobian diagnostics for $(I_E/I_W, I_N/I_E)$ in $\log_{10}(\gamma)$ space")
fig.savefig(OUT_FIG, dpi=220)
print(f"Saved figure: {OUT_FIG}")
