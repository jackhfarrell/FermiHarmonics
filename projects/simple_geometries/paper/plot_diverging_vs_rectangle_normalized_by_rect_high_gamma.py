from pathlib import Path
import csv

import matplotlib.pyplot as plt
import numpy as np

PAPER_DIR = Path(__file__).resolve().parent
RECT_CSV = PAPER_DIR / "rectangle_sweep" / "rectangle_conductance_vs_gamma_mc.csv"
DIV_CSV = PAPER_DIR / "diverging_sweep" / "diverging_conductance_vs_gamma_mc.csv"
OUTFILE = PAPER_DIR / "diverging_vs_rectangle_normalized_by_rect_high_gamma.png"


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


def main():
    gamma_rect, g_rect = load_curve(RECT_CSV)
    gamma_div, g_div = load_curve(DIV_CSV)

    g_rect_high = float(g_rect[-1])  # high-gamma reference (max gamma_mc in rectangle sweep)
    if not np.isfinite(g_rect_high) or g_rect_high <= 0:
        raise RuntimeError("Rectangle high-gamma conductance is invalid for normalization.")

    g_rect_norm = g_rect / g_rect_high
    g_div_norm = g_div / g_rect_high

    r_rect_norm = g_rect_high / g_rect
    r_div_norm = g_rect_high / g_div

    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.4), constrained_layout=True)

    ax = axes[0]
    ax.plot(gamma_rect, g_rect_norm, "-o", lw=1.8, ms=4.0, label="Rectangle")
    ax.plot(gamma_div, g_div_norm, "-o", lw=1.8, ms=4.0, label="Diverging")
    ax.set_xscale("log")
    ax.set_xlabel(r"$\gamma_{mc}$")
    ax.set_ylabel(r"$G / G_{\mathrm{rect}}(\gamma_{mc,\max})$")
    ax.set_title("Conductance Normalized By Rectangle High-$\\gamma_{mc}$")
    ax.grid(True, which="both", alpha=0.28)
    ax.legend(frameon=False)

    ax = axes[1]
    ax.plot(gamma_rect, r_rect_norm, "-o", lw=1.8, ms=4.0, label="Rectangle")
    ax.plot(gamma_div, r_div_norm, "-o", lw=1.8, ms=4.0, label="Diverging")
    ax.set_xscale("log")
    ax.set_xlabel(r"$\gamma_{mc}$")
    ax.set_ylabel(r"$R / R_{\mathrm{rect}}(\gamma_{mc,\max})$")
    ax.set_title("Resistance View (Same Normalization)")
    ax.grid(True, which="both", alpha=0.28)
    ax.legend(frameon=False)

    fig.savefig(OUTFILE, dpi=240)
    print(f"Rectangle high-gamma reference G = {g_rect_high:.8g}")
    print(f"Saved: {OUTFILE}")


if __name__ == "__main__":
    main()
