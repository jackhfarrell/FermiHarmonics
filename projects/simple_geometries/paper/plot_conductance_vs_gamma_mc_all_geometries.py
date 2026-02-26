from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

PAPER_DIR = Path(__file__).resolve().parent

SOURCES = [
    ("Diverging nozzle", PAPER_DIR / "diverging_sweep" / "diverging_conductance_vs_gamma_mc.csv"),
    ("Rectangle", PAPER_DIR / "rectangle_sweep" / "rectangle_conductance_vs_gamma_mc.csv"),
    ("Dog-leg", PAPER_DIR / "dogleg_sweep" / "dogleg_conductance_vs_gamma_mc.csv"),
]

OUTFILE = PAPER_DIR / "conductance_vs_gamma_mc_all_geometries.png"


def load_csv(path: Path):
    data = np.genfromtxt(path, delimiter=",", names=True)
    gamma_mc = np.asarray(data["gamma_mc"], dtype=float)
    conductance = np.asarray(data["conductance"], dtype=float)
    order = np.argsort(gamma_mc)
    return gamma_mc[order], conductance[order]


def main():
    fig, ax = plt.subplots(figsize=(7.6, 4.8), constrained_layout=True)

    for label, csv_path in SOURCES:
        if not csv_path.is_file():
            raise RuntimeError(f"Missing CSV: {csv_path}")
        gamma_mc, conductance = load_csv(csv_path)
        ax.plot(gamma_mc, conductance, "-o", lw=1.8, ms=4.0, label=label)

    ax.set_xscale("log")
    ax.set_xlabel(r"$\gamma_{mc}$")
    ax.set_ylabel(r"$G = I_{\mathrm{norm}}/|\Delta a_0|$")
    ax.set_title(r"Conductance vs $\gamma_{mc}$ (all geometries)")
    ax.grid(True, which="both", alpha=0.28)
    ax.legend(frameon=False)

    fig.savefig(OUTFILE, dpi=240)
    print(f"Saved: {OUTFILE}")


if __name__ == "__main__":
    main()
