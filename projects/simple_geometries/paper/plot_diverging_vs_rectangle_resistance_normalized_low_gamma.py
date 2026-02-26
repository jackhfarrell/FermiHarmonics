from pathlib import Path
import csv

import matplotlib.pyplot as plt
import numpy as np

PAPER_DIR = Path(__file__).resolve().parent
RECT_CSV = PAPER_DIR / "rectangle_sweep" / "rectangle_conductance_vs_gamma_mc.csv"
DIV_CSV = PAPER_DIR / "diverging_sweep" / "diverging_conductance_vs_gamma_mc.csv"
OUTFILE = PAPER_DIR / "diverging_vs_rectangle_resistance_normalized_low_gamma.png"


def load_curve(csv_path: Path):
    gamma_mc = []
    conductance = []
    with csv_path.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            gmc = float(row["gamma_mc"])
            g = float(row["conductance"])
            if np.isfinite(g) and g > 0:
                gamma_mc.append(gmc)
                conductance.append(g)
    gamma_mc = np.asarray(gamma_mc, dtype=float)
    conductance = np.asarray(conductance, dtype=float)
    order = np.argsort(gamma_mc)
    return gamma_mc[order], conductance[order]


def main():
    gamma_rect, g_rect = load_curve(RECT_CSV)
    gamma_div, g_div = load_curve(DIV_CSV)

    r_rect = 1.0 / g_rect
    r_div = 1.0 / g_div

    r_rect_norm = r_rect / r_rect[0]  # normalize by each geometry's low-gamma resistance
    r_div_norm = r_div / r_div[0]

    fig, ax = plt.subplots(figsize=(7.6, 4.6), constrained_layout=True)
    ax.plot(gamma_rect, r_rect_norm, "-o", lw=1.8, ms=4.0, label="Rectangle")
    ax.plot(gamma_div, r_div_norm, "-o", lw=1.8, ms=4.0, label="Diverging")
    ax.set_xscale("log")
    ax.set_xlabel(r"$\gamma_{mc}$")
    ax.set_ylabel(r"$R(\gamma_{mc}) / R(\gamma_{mc,\min})$")
    ax.set_title("Resistance Normalized By Low-$\\gamma_{mc}$ Value")
    ax.grid(True, which="both", alpha=0.28)
    ax.legend(frameon=False)

    fig.savefig(OUTFILE, dpi=240)
    print(f"Saved: {OUTFILE}")


if __name__ == "__main__":
    main()
