from pathlib import Path
import re

import h5py
import matplotlib.pyplot as plt
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
RESULTS_ROOT = PROJECT_ROOT / "results"
RESULTS_ROOT_SWEEPS = PROJECT_ROOT / "results" / "sweeps"
SWEEP_GLOBS = [
    "simple_geometries_rectangle_gamma_mc_sweep_sweep_*",
    "simple_geometries_rectangle_gamma_mc_sweep_*",
    "simple_geometries_rectangle_sweep_*",
]
HARDCODED_SWEEP_DIR = RESULTS_ROOT / "simple_geometries_rectangle_gamma_mc_sweep_sweep_2026-02-25_142706"

GAMMA_MR_TARGET = 1e-2
GAMMA_MR_TOL = 1e-12
CONTACT_SAMPLE_DEPTH = 2
A0_BULK_LOW_FRAC = 0.30
A0_BULK_HIGH_FRAC = 0.70

OUT_CSV = Path(__file__).resolve().parent / "rectangle_conductance_vs_gamma_mc.csv"
OUT_PNG = Path(__file__).resolve().parent / "rectangle_conductance_vs_gamma_mc.png"
OUT_PNG_ZOOM = Path(__file__).resolve().parent / "rectangle_conductance_vs_gamma_mc_zoom.png"


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


def line_integral_abs_current(x: np.ndarray, field_line: np.ndarray, valid: np.ndarray) -> float:
    return np.abs(integrate_over_valid_segments(x, field_line, valid))


def line_span(x: np.ndarray, valid: np.ndarray) -> float:
    if np.count_nonzero(valid) < 2:
        return np.nan
    xv = x[valid]
    return float(np.max(xv) - np.min(xv))


def pick_latest_sweep_dir():
    if HARDCODED_SWEEP_DIR.is_dir():
        return HARDCODED_SWEEP_DIR
    candidates = []
    for root in (RESULTS_ROOT, RESULTS_ROOT_SWEEPS):
        for pattern in SWEEP_GLOBS:
            candidates.extend([p for p in root.glob(pattern) if p.is_dir()])
    if not candidates:
        raise RuntimeError(f"No sweep directories matching {SWEEP_GLOBS} found in {RESULTS_ROOT} or {RESULTS_ROOT_SWEEPS}")
    return sorted(candidates, key=lambda p: p.stat().st_mtime)[-1]


def parse_rows(files):
    pat = re.compile(r"gamma_mc=([0-9eE+\-.]+)_gamma_mr=([0-9eE+\-.]+)")
    rows = []
    for fpath in files:
        m = pat.search(fpath.name)
        if m is None:
            continue
        gamma_mc = float(m.group(1))
        gamma_mr = float(m.group(2))
        if not np.isclose(gamma_mr, GAMMA_MR_TARGET, atol=GAMMA_MR_TOL, rtol=0.0):
            continue

        with h5py.File(fpath, "r") as h5:
            a0 = np.asarray(h5["a0"]).T
            b1 = np.asarray(h5["b1"]).T
            x = np.asarray(h5["x"])
            mask = np.asarray(h5["mask"]).astype(bool).T

        y_valid = np.where(np.any(mask, axis=0))[0]
        if len(y_valid) <= 2 * CONTACT_SAMPLE_DEPTH:
            continue
        iy_bottom = int(y_valid[CONTACT_SAMPLE_DEPTH])
        iy_top = int(y_valid[-1 - CONTACT_SAMPLE_DEPTH])

        valid_bottom = mask[:, iy_bottom] & np.isfinite(b1[:, iy_bottom])
        valid_top = mask[:, iy_top] & np.isfinite(b1[:, iy_top])
        if np.count_nonzero(valid_bottom) < 2 or np.count_nonzero(valid_top) < 2:
            continue

        i_bottom = line_integral_abs_current(x, b1[:, iy_bottom], valid_bottom)
        i_top = line_integral_abs_current(x, b1[:, iy_top], valid_top)
        i_total = 0.5 * (i_bottom + i_top)

        w_bottom = line_span(x, valid_bottom)
        w_top = line_span(x, valid_top)
        width_eff = 0.5 * (w_bottom + w_top)
        i_norm = i_total / width_eff if np.isfinite(width_eff) and width_eff > 0 else np.nan

        iy_bulk_lo = int(y_valid[int(round(A0_BULK_LOW_FRAC * (len(y_valid) - 1)))])
        iy_bulk_hi = int(y_valid[int(round(A0_BULK_HIGH_FRAC * (len(y_valid) - 1)))])
        a0_bottom = np.nanmean(np.where(mask[:, iy_bulk_lo], a0[:, iy_bulk_lo], np.nan))
        a0_top = np.nanmean(np.where(mask[:, iy_bulk_hi], a0[:, iy_bulk_hi], np.nan))
        delta_a0 = np.abs(a0_bottom - a0_top)

        conductance = i_norm / delta_a0 if np.isfinite(i_norm) and np.isfinite(delta_a0) and delta_a0 > 0 else np.nan
        rows.append((gamma_mc, gamma_mr, i_total, i_norm, delta_a0, conductance))

    return rows


def write_csv(rows):
    with OUT_CSV.open("w", encoding="utf-8") as f:
        f.write("gamma_mc,gamma_mr,i_total,i_norm,delta_a0_measured,conductance\n")
        for r in rows:
            f.write(",".join(f"{v:.16g}" for v in r) + "\n")


def main():
    sweep_dir = pick_latest_sweep_dir()
    files = sorted((sweep_dir / "data").glob("*.h5"))
    if not files:
        raise RuntimeError(f"No .h5 files found in {sweep_dir / 'data'}")
    rows = parse_rows(files)
    if not rows:
        raise RuntimeError(f"No usable rows found for gamma_mr={GAMMA_MR_TARGET}")
    rows = sorted(rows, key=lambda t: t[0])
    write_csv(rows)

    gamma_mc = np.array([r[0] for r in rows])
    conductance = np.array([r[5] for r in rows])

    fig, ax = plt.subplots(figsize=(7.2, 4.6), constrained_layout=True)
    ax.plot(gamma_mc, conductance, "-o", lw=1.8, ms=4.5)
    ax.set_xscale("log")
    ax.set_xlabel(r"$\gamma_{mc}$")
    ax.set_ylabel(r"$G = I_{\mathrm{norm}}/|\Delta a_0|$")
    ax.set_title(r"Rectangle Conductance vs $\gamma_{mc}$")
    ax.grid(True, which="both", alpha=0.28)
    fig.savefig(OUT_PNG, dpi=220)

    # Zoomed view on low-to-mid gamma_mc region for easier visual inspection.
    fig_zoom, ax_zoom = plt.subplots(figsize=(7.2, 4.6), constrained_layout=True)
    ax_zoom.plot(gamma_mc, conductance, "-o", lw=1.8, ms=4.5)
    ax_zoom.set_xlabel(r"$\gamma_{mc}$")
    ax_zoom.set_ylabel(r"$G = I_{\mathrm{norm}}/|\Delta a_0|$")
    ax_zoom.set_title(r"Rectangle Conductance vs $\gamma_{mc}$ (Zoom)")
    ax_zoom.set_xlim(0.0, 20.0)
    finite = conductance[np.isfinite(conductance)]
    if finite.size:
        y0 = float(np.min(finite))
        y1 = float(np.max(finite[gamma_mc <= 20.0])) if np.any(gamma_mc <= 20.0) else float(np.max(finite))
        pad = 0.08 * max(y1 - y0, 1e-12)
        ax_zoom.set_ylim(y0 - pad, y1 + pad)
    ax_zoom.grid(True, which="both", alpha=0.28)
    fig_zoom.savefig(OUT_PNG_ZOOM, dpi=220)

    print(f"Using sweep directory: {sweep_dir}")
    print(f"Rows used: {len(rows)}")
    print(f"Saved CSV: {OUT_CSV}")
    print(f"Saved plot: {OUT_PNG}")
    print(f"Saved zoom plot: {OUT_PNG_ZOOM}")


if __name__ == "__main__":
    main()
