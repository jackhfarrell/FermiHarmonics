from pathlib import Path
import csv
import re

import h5py
import matplotlib.pyplot as plt
import numpy as np

PAPER_DIR = Path(__file__).resolve().parent
RESULTS_ROOT = PAPER_DIR.parent / "results"

GEOMS = {
    "rectangle": {
        "csv": PAPER_DIR / "rectangle_sweep" / "rectangle_conductance_vs_gamma_mc.csv",
        "sweep_dir": RESULTS_ROOT / "simple_geometries_rectangle_gamma_mc_sweep_sweep_2026-02-25_142706",
    },
    "diverging": {
        "csv": PAPER_DIR / "diverging_sweep" / "diverging_conductance_vs_gamma_mc.csv",
        "sweep_dir": RESULTS_ROOT / "simple_geometries_diverging_nozzle_gamma_mc_sweep_sweep_2026-02-25_142715",
    },
    "dogleg": {
        "csv": PAPER_DIR / "dogleg_sweep" / "dogleg_conductance_vs_gamma_mc.csv",
        "sweep_dir": RESULTS_ROOT / "simple_geometries_dogleg_gamma_mc_sweep_sweep_2026-02-25_145852",
    },
}

OUT_PNG = PAPER_DIR / "ballistic_diffuse_perimeter_law_check.png"
OUT_CSV = PAPER_DIR / "ballistic_diffuse_perimeter_law_summary.csv"


def low_gamma_conductance(csv_path: Path):
    rows = []
    with csv_path.open("r", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            gmc = float(row["gamma_mc"])
            g = float(row["conductance"])
            if np.isfinite(g) and g > 0:
                rows.append((gmc, g))
    if not rows:
        raise RuntimeError(f"No finite conductance values in {csv_path}")
    rows.sort(key=lambda t: t[0])
    return rows[0]


def first_h5(sweep_dir: Path):
    files = sorted((sweep_dir / "data").glob("*.h5"))
    if not files:
        raise RuntimeError(f"No h5 files in {sweep_dir / 'data'}")
    # prefer smallest gamma_mc file if available
    pat = re.compile(r"gamma_mc=([0-9eE+\-.]+)")
    tagged = []
    for f in files:
        m = pat.search(f.name)
        if m:
            tagged.append((float(m.group(1)), f))
    if tagged:
        tagged.sort(key=lambda t: t[0])
        return tagged[0][1]
    return files[0]


def area_perimeter_from_mask(h5_path: Path):
    with h5py.File(h5_path, "r") as h5:
        x = np.asarray(h5["x"], dtype=float)
        y = np.asarray(h5["y"], dtype=float)
        mask = np.asarray(h5["mask"]).astype(bool).T

    dx = float(np.mean(np.diff(x)))
    dy = float(np.mean(np.diff(y)))
    area = float(np.count_nonzero(mask) * dx * dy)

    # Grid-edge perimeter estimate: count domain/non-domain transitions,
    # including outer boundary by padding with False.
    padded = np.pad(mask, ((1, 1), (1, 1)), mode="constant", constant_values=False)
    p = 0.0
    # vertical grid edges (normal x)
    trans_x = padded[:-1, :] != padded[1:, :]
    p += float(np.count_nonzero(trans_x) * dy)
    # horizontal grid edges (normal y)
    trans_y = padded[:, :-1] != padded[:, 1:]
    p += float(np.count_nonzero(trans_y) * dx)

    return area, p


def main():
    rows = []
    for name, cfg in GEOMS.items():
        gamma_mc_min, g0 = low_gamma_conductance(cfg["csv"])
        h5_path = first_h5(cfg["sweep_dir"])
        area, perim = area_perimeter_from_mask(h5_path)
        ap = area / perim if perim > 0 else np.nan
        rows.append((name, gamma_mc_min, g0, area, perim, ap))

    rows.sort(key=lambda t: t[0])

    # Fit G ~ c*(A/P) through origin as quick perimeter-law check.
    ap_vec = np.array([r[5] for r in rows], dtype=float)
    g_vec = np.array([r[2] for r in rows], dtype=float)
    c = float(np.dot(ap_vec, g_vec) / np.dot(ap_vec, ap_vec))
    g_pred = c * ap_vec

    ref_idx = [i for i, r in enumerate(rows) if r[0] == "rectangle"][0]
    g_ref = g_vec[ref_idx]
    gpred_ref = g_pred[ref_idx]

    with OUT_CSV.open("w", encoding="utf-8") as f:
        f.write("geometry,gamma_mc_min,G_low,A,P,A_over_P,G_pred,ratio_to_rect,pred_ratio_to_rect\n")
        for i, r in enumerate(rows):
            f.write(
                f"{r[0]},{r[1]:.16g},{r[2]:.16g},{r[3]:.16g},{r[4]:.16g},{r[5]:.16g},"
                f"{g_pred[i]:.16g},{(r[2]/g_ref):.16g},{(g_pred[i]/gpred_ref):.16g}\n"
            )

    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.2), constrained_layout=True)

    # Panel 1: normalized observed vs perimeter-law prediction.
    names = [r[0] for r in rows]
    x = np.arange(len(rows))
    obs_norm = g_vec / g_ref
    pred_norm = g_pred / gpred_ref
    axes[0].plot(x, obs_norm, "o-", lw=1.8, label="Observed low-$\\gamma_{mc}$ G")
    axes[0].plot(x, pred_norm, "s--", lw=1.6, label=r"Perimeter-law fit $\propto A/P$")
    axes[0].set_xticks(x, names)
    axes[0].set_ylabel("Normalized to rectangle")
    axes[0].set_title("Ballistic-Diffuse Scaling Check")
    axes[0].grid(True, alpha=0.28)
    axes[0].legend(frameon=False, fontsize=9)

    # Panel 2: raw G vs A/P with fit line.
    xx = np.linspace(0.0, 1.05 * np.max(ap_vec), 200)
    axes[1].plot(ap_vec, g_vec, "o", ms=6)
    for i, n in enumerate(names):
        axes[1].annotate(n, (ap_vec[i], g_vec[i]), textcoords="offset points", xytext=(4, 4), fontsize=9)
    axes[1].plot(xx, c * xx, "-", lw=1.6, label=rf"fit: $G \approx {c:.3g}(A/P)$")
    axes[1].set_xlabel(r"$A/P$")
    axes[1].set_ylabel(r"$G$ at lowest $\gamma_{mc}$")
    axes[1].set_title(r"Low-$\gamma_{mc}$ Conductance vs $A/P$")
    axes[1].grid(True, alpha=0.28)
    axes[1].legend(frameon=False, fontsize=9)

    fig.savefig(OUT_PNG, dpi=240)
    print(f"Saved: {OUT_PNG}")
    print(f"Saved: {OUT_CSV}")


if __name__ == "__main__":
    main()
