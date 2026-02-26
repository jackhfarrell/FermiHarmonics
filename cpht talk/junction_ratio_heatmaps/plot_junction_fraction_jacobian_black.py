from pathlib import Path
import re

import h5py
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import seaborn as sns

plt.rcParams.update(
    {
        "text.usetex": True,
        "font.family": "serif",
        "font.serif": ["Computer Modern Roman", "Times New Roman", "DejaVu Serif"],
        "font.size": 24,
        "axes.titlesize": 24,
        "axes.labelsize": 24,
        "xtick.labelsize": 20,
        "ytick.labelsize": 20,
    }
)

# Data/source configuration
REPO_ROOT = Path(__file__).resolve().parents[2]
RESULTS_ROOT = REPO_ROOT / "projects" / "simple_geometries" / "results" / "sweeps"
SWEEP_GLOB = "simple_geometries_junction_sweep_*"
X_CUT_WEST = -0.4
X_CUT_EAST = 0.4
Y_CUT_NORTH = 0.4
GAMMA_MR_MIN = 1e-2
GAMMA_MR_MAX = 1e1
GAMMA_MC_MIN = 1e-2
GAMMA_MC_MAX = 1e2
SMOOTHING_PASSES = 0
JACOBIAN_SMOOTHING_PASSES = 1

# Output configuration
OUT_DIR = Path(__file__).resolve().parent
OUT_JAC = OUT_DIR / "junction_fraction_jacobian_det_black.png"


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
    kernel = np.ones((3, 3), dtype=float)
    out = grid.copy()
    for _ in range(max(0, passes)):
        padded = np.pad(out, 1, mode="edge")
        next_out = out.copy()
        for i in range(out.shape[0]):
            for j in range(out.shape[1]):
                win = padded[i : i + 3, j : j + 3]
                mask = np.isfinite(win)
                if np.any(mask):
                    next_out[i, j] = np.sum(win[mask] * kernel[mask]) / np.sum(kernel[mask])
                else:
                    next_out[i, j] = np.nan
        out = next_out
    return out


def safe_ratio(num: float, den: float) -> float:
    if not np.isfinite(num) or not np.isfinite(den) or np.isclose(den, 0.0):
        return np.nan
    return num / den


def finite_differences_log(grid: np.ndarray, logx: np.ndarray, logy: np.ndarray):
    nx, ny = grid.shape
    dlogx = np.full_like(grid, np.nan, dtype=float)
    dlogy = np.full_like(grid, np.nan, dtype=float)

    for i in range(nx):
        for j in range(ny):
            if not np.isfinite(grid[i, j]):
                continue

            # d/dlogx
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

            # d/dlogy
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


def style_dark():
    plt.rcParams.update(
        {
            "figure.facecolor": "black",
            "axes.facecolor": "black",
            "savefig.facecolor": "black",
            "axes.edgecolor": "white",
            "axes.labelcolor": "white",
            "xtick.color": "white",
            "ytick.color": "white",
            "text.color": "white",
            "axes.titlecolor": "white",
            "font.size": 24,
            "axes.titlesize": 24,
            "axes.labelsize": 24,
            "xtick.labelsize": 20,
            "ytick.labelsize": 20,
        }
    )


def build_rows(files):
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
        i_w = np.abs(integrate_over_valid_segments(y, a1_west_line, valid_y_west)) if np.count_nonzero(valid_y_west) >= 2 else np.nan

        a1_east_line = np.where(mask[ix_east, :], a1[ix_east, :], np.nan)
        valid_y_east = np.isfinite(a1_east_line)
        i_e = np.abs(integrate_over_valid_segments(y, a1_east_line, valid_y_east)) if np.count_nonzero(valid_y_east) >= 2 else np.nan

        b1_north_line = np.where(mask[:, iy_north], b1[:, iy_north], np.nan)
        valid_x_north = np.isfinite(b1_north_line)
        i_n = np.abs(integrate_over_valid_segments(x, b1_north_line, valid_x_north)) if np.count_nonzero(valid_x_north) >= 2 else np.nan

        rows.append((gamma_mr, gamma_mc, safe_ratio(i_e, i_w), safe_ratio(i_n, i_e)))
    return rows


def main():
    style_dark()

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

    rows = build_rows(files)
    rows = [r for r in rows if GAMMA_MR_MIN <= r[0] <= GAMMA_MR_MAX and GAMMA_MC_MIN <= r[1] <= GAMMA_MC_MAX]
    if not rows:
        raise RuntimeError("No rows remain after gamma truncation.")

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

    # Apply the same smoothing style as ratio heatmaps before derivatives.
    ie_over_iw_grid = smooth_nan_grid(ie_over_iw_grid, passes=SMOOTHING_PASSES)
    in_over_ie_grid = smooth_nan_grid(in_over_ie_grid, passes=SMOOTHING_PASSES)

    d_ie_iw_dlog_mr, d_ie_iw_dlog_mc = finite_differences_log(ie_over_iw_grid, log_gamma_mr, log_gamma_mc)
    d_in_ie_dlog_mr, d_in_ie_dlog_mc = finite_differences_log(in_over_ie_grid, log_gamma_mr, log_gamma_mc)

    det_j = d_ie_iw_dlog_mr * d_in_ie_dlog_mc - d_ie_iw_dlog_mc * d_in_ie_dlog_mr
    abs_det_j = np.abs(det_j)
    abs_det_j = smooth_nan_grid(abs_det_j, passes=JACOBIAN_SMOOTHING_PASSES)
    finite = abs_det_j[np.isfinite(abs_det_j)]
    if finite.size == 0:
        raise RuntimeError("No finite values for det(J)")

    vmax = float(np.nanpercentile(finite, 98))
    if not np.isfinite(vmax) or vmax <= 0:
        vmax = float(np.nanmax(finite))
        if not np.isfinite(vmax) or vmax <= 0:
            vmax = 1.0
    vmin = 0.0

    cmap = sns.color_palette("mako", as_cmap=True).copy()
    cmap.set_bad("#202020")

    fig, ax = plt.subplots(figsize=(5.2, 5.2), constrained_layout=True)
    pcm = ax.pcolormesh(
        gamma_mr_vals,
        gamma_mc_vals,
        np.ma.masked_invalid(abs_det_j.T),
        shading="auto",
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(GAMMA_MR_MIN, GAMMA_MR_MAX)
    ax.set_ylim(GAMMA_MC_MIN, GAMMA_MC_MAX)
    ax.set_box_aspect(1)
    ax.set_xlabel(r"$\gamma_{mr}$")
    ax.set_ylabel(r"$\gamma_{mc}$")
    for spine in ax.spines.values():
        spine.set_color("white")
        spine.set_linewidth(1.2)

    cbar = fig.colorbar(pcm, ax=ax, shrink=0.78, fraction=0.055, pad=0.02)
    cbar.outline.set_edgecolor("white")
    cbar.ax.yaxis.set_tick_params(color="white")
    plt.setp(cbar.ax.get_yticklabels(), color="white")
    cbar.formatter = mticker.FormatStrFormatter("%.2g")
    cbar.update_ticks()
    cbar.set_label(r"$|\det J|$", color="white")

    fig.savefig(OUT_JAC, dpi=1200, facecolor=fig.get_facecolor())
    print(f"Saved: {OUT_JAC}")


if __name__ == "__main__":
    main()
