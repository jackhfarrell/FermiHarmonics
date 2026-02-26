from pathlib import Path
import csv

import matplotlib.pyplot as plt
import numpy as np

PAPER_DIR = Path(__file__).resolve().parent
RECT_CSV = PAPER_DIR / "rectangle_sweep" / "rectangle_conductance_vs_gamma_mc.csv"
DIV_CSV = PAPER_DIR / "diverging_sweep" / "diverging_conductance_vs_gamma_mc.csv"
DOG_CSV = PAPER_DIR / "dogleg_sweep" / "dogleg_conductance_vs_gamma_mc.csv"
OUTFILE = PAPER_DIR / "all_geometries_conductance_normalized_high_gamma.png"


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
    gamma_dog, g_dog = load_curve(DOG_CSV)

    # Normalize by each geometry's highest-viscosity (max gamma_mc) value.
    g_rect_norm = g_rect / g_rect[-1]
    g_div_norm = g_div / g_div[-1]
    g_dog_norm = g_dog / g_dog[-1]

    fig, ax = plt.subplots(figsize=(7.8, 4.8), constrained_layout=True)
    ax.plot(gamma_rect, g_rect_norm, "-o", lw=1.8, ms=4.0, label="Rectangle")
    ax.plot(gamma_div, g_div_norm, "-o", lw=1.8, ms=4.0, label="Diverging")
    ax.plot(gamma_dog, g_dog_norm, "-o", lw=1.8, ms=4.0, label="Dog-leg")
    ax.set_xscale("log")
    ax.set_xlabel(r"$\gamma_{mc}$")
    ax.set_ylabel(r"$G(\gamma_{mc}) / G(\gamma_{mc,\max})$")
    ax.set_title("Conductance Normalized at Highest Viscosity")
    ax.grid(True, which="both", alpha=0.28)
    ax.legend(frameon=False)

    fig.savefig(OUTFILE, dpi=240)
    print(f"Saved: {OUTFILE}")


if __name__ == "__main__":
    main()
