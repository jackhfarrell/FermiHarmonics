from pathlib import Path
import re

import h5py
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

# Data/source configuration
REPO_ROOT = Path(__file__).resolve().parents[2]
RESULTS_ROOT = REPO_ROOT / "projects" / "simple_geometries" / "results" / "sweeps"
SWEEP_GLOB = "simple_geometries_junction_sweep_*"
X_CUT_WEST = -0.4
X_CUT_EAST = 0.4
Y_CUT_NORTH = 0.4
GAMMA_MR_MIN = 1e-2
GAMMA_MR_MAX = 1e2
GAMMA_MC_MIN = 1e-2
GAMMA_MC_MAX = 1e2
SMOOTHING_PASSES = 0

# Output configuration
OUT_DIR = Path(__file__).resolve().parent
OUT_IN_OVER_IE = OUT_DIR / "junction_in_over_ie_heatmap_black.png"
OUT_IW_OVER_IE = OUT_DIR / "junction_iw_over_ie_heatmap_black.png"


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

        rows.append((gamma_mr, gamma_mc, safe_ratio(i_n, i_e), safe_ratio(i_w, i_e)))
    return rows


def build_grid(rows, col_idx):
    rows = [r for r in rows if GAMMA_MR_MIN <= r[0] <= GAMMA_MR_MAX and GAMMA_MC_MIN <= r[1] <= GAMMA_MC_MAX]
    if not rows:
        raise RuntimeError("No rows remain after gamma_mr/gamma_mc truncation.")

    gamma_mr_vals = np.array(sorted({r[0] for r in rows}))
    gamma_mc_vals = np.array(sorted({r[1] for r in rows}))
    grid = np.full((len(gamma_mr_vals), len(gamma_mc_vals)), np.nan)

    for row in rows:
        gamma_mr, gamma_mc = row[0], row[1]
        i = np.searchsorted(gamma_mr_vals, gamma_mr)
        j = np.searchsorted(gamma_mc_vals, gamma_mc)
        grid[i, j] = row[col_idx]
    return gamma_mr_vals, gamma_mc_vals, smooth_nan_grid(grid, passes=SMOOTHING_PASSES)


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
            "font.size": 16,
            "axes.titlesize": 20,
            "axes.labelsize": 18,
            "xtick.labelsize": 14,
            "ytick.labelsize": 14,
        }
    )


def plot_single(gamma_mr_vals, gamma_mc_vals, grid, title, outfile):
    finite = grid[np.isfinite(grid)]
    if finite.size == 0:
        raise RuntimeError(f"No finite values to plot for {title}.")

    vmin = float(np.quantile(finite, 0.01))
    vmax = float(np.quantile(finite, 0.99))
    if not np.isfinite(vmin) or not np.isfinite(vmax) or vmax <= vmin:
        vmin, vmax = float(np.nanmin(finite)), float(np.nanmax(finite))
        if vmax <= vmin:
            vmin, vmax = 0.0, 1.0

    cmap = plt.get_cmap("magma").copy()
    cmap.set_bad("#202020")

    fig, ax = plt.subplots(figsize=(7.2, 6.4), constrained_layout=True)
    z = np.ma.masked_invalid(grid.T)
    levels = np.linspace(vmin, vmax, 80)
    pcm = ax.contourf(
        gamma_mr_vals,
        gamma_mc_vals,
        z,
        levels=levels,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        extend="both",
    )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(GAMMA_MR_MIN, GAMMA_MR_MAX)
    ax.set_ylim(GAMMA_MC_MIN, GAMMA_MC_MAX)
    ax.set_xlabel(r"$\gamma_{mr}$")
    ax.set_ylabel(r"$\gamma_{mc}$")
    ax.set_title(title, pad=10)
    for spine in ax.spines.values():
        spine.set_color("white")
        spine.set_linewidth(1.2)

    cbar = fig.colorbar(pcm, ax=ax)
    cbar.outline.set_edgecolor("white")
    cbar.ax.yaxis.set_tick_params(color="white")
    plt.setp(cbar.ax.get_yticklabels(), color="white")
    cbar.formatter = mticker.FormatStrFormatter("%.2g")
    cbar.update_ticks()
    cbar.set_label(title, color="white")

    fig.savefig(outfile, dpi=260, facecolor=fig.get_facecolor())
    print(f"Saved: {outfile}")


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
    if not rows:
        raise RuntimeError("No parseable rows found (filename regex did not match).")

    gamma_mr_vals, gamma_mc_vals, in_over_ie_grid = build_grid(rows, col_idx=2)
    _, _, iw_over_ie_grid = build_grid(rows, col_idx=3)

    plot_single(gamma_mr_vals, gamma_mc_vals, in_over_ie_grid, r"$I_N / I_E$", OUT_IN_OVER_IE)
    plot_single(gamma_mr_vals, gamma_mc_vals, iw_over_ie_grid, r"$I_W / I_E$", OUT_IW_OVER_IE)


if __name__ == "__main__":
    main()
