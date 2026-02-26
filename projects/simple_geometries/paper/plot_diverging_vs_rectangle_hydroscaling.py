from pathlib import Path
import csv

import h5py
import matplotlib.pyplot as plt
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS_ROOT = PROJECT_ROOT / "results"

RECT_SWEEP = RESULTS_ROOT / "simple_geometries_rectangle_gamma_mc_sweep_sweep_2026-02-25_142706"
DIV_SWEEP = RESULTS_ROOT / "simple_geometries_diverging_nozzle_gamma_mc_sweep_sweep_2026-02-25_142715"

RECT_CSV = Path(__file__).resolve().parent / "rectangle_sweep" / "rectangle_conductance_vs_gamma_mc.csv"
DIV_CSV = Path(__file__).resolve().parent / "diverging_sweep" / "diverging_conductance_vs_gamma_mc.csv"

OUT_PNG = Path(__file__).resolve().parent / "diverging_vs_rectangle_hydroscaling.png"


def load_curve(csv_path: Path):
    gamma_mc = []
    conductance = []
    with csv_path.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            gmc = float(row["gamma_mc"])
            g = float(row["conductance"])
            if np.isfinite(g):
                gamma_mc.append(gmc)
                conductance.append(g)
    gamma_mc = np.asarray(gamma_mc, dtype=float)
    conductance = np.asarray(conductance, dtype=float)
    order = np.argsort(gamma_mc)
    return gamma_mc[order], conductance[order]


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


def geometric_resistance_metric(sweep_dir: Path) -> float:
    files = sorted((sweep_dir / "data").glob("*.h5"))
    if not files:
        raise RuntimeError(f"No data files in {sweep_dir / 'data'}")
    fpath = files[0]
    with h5py.File(fpath, "r") as h5:
        x = np.asarray(h5["x"])
        y = np.asarray(h5["y"])
        mask = np.asarray(h5["mask"]).astype(bool).T

    widths = np.full(len(y), np.nan, dtype=float)
    for j in range(len(y)):
        valid = mask[:, j]
        if np.count_nonzero(valid) >= 2:
            xv = x[valid]
            widths[j] = float(np.max(xv) - np.min(xv))

    invw = np.where(np.isfinite(widths) & (widths > 0), 1.0 / widths, np.nan)
    valid = np.isfinite(invw)
    if np.count_nonzero(valid) < 2:
        raise RuntimeError(f"Could not compute width profile integral from {fpath}")
    return integrate_over_valid_segments(y, invw, valid)


def main():
    gamma_rect, g_rect = load_curve(RECT_CSV)
    gamma_div, g_div = load_curve(DIV_CSV)

    r_rect = geometric_resistance_metric(RECT_SWEEP)
    r_div = geometric_resistance_metric(DIV_SWEEP)

    g_rect_scaled = g_rect * r_rect
    g_div_scaled = g_div * r_div

    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.4), constrained_layout=True)

    ax = axes[0]
    ax.plot(gamma_rect, g_rect, "-o", lw=1.8, ms=4.0, label="Rectangle")
    ax.plot(gamma_div, g_div, "-o", lw=1.8, ms=4.0, label="Diverging")
    ax.set_xscale("log")
    ax.set_xlabel(r"$\gamma_{mc}$")
    ax.set_ylabel(r"$G$")
    ax.set_title("Raw Conductance")
    ax.grid(True, which="both", alpha=0.28)
    ax.legend(frameon=False)

    ax = axes[1]
    ax.plot(gamma_rect, g_rect_scaled, "-o", lw=1.8, ms=4.0, label="Rectangle")
    ax.plot(gamma_div, g_div_scaled, "-o", lw=1.8, ms=4.0, label="Diverging")
    ax.set_xscale("log")
    ax.set_xlabel(r"$\gamma_{mc}$")
    ax.set_ylabel(r"$G \times \int dy/w(y)$")
    ax.set_title(r"Hydrodynamic Geometric Scaling")
    ax.grid(True, which="both", alpha=0.28)
    ax.legend(frameon=False)

    fig.savefig(OUT_PNG, dpi=240)

    print(f"Rectangle integral ∫dy/w(y) = {r_rect:.8g}")
    print(f"Diverging integral ∫dy/w(y) = {r_div:.8g}")
    print(f"Saved: {OUT_PNG}")


if __name__ == "__main__":
    main()
