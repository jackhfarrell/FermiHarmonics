from pathlib import Path
import re

import h5py
import matplotlib.pyplot as plt
import numpy as np

# Simple analysis script for diverging-nozzle gamma_mc sweeps.
# Conductance definition used here:
#   G = I_norm / |Δa0_measured|
# where I_norm is the contact-averaged total current normalized by contact width,
# and Δa0_measured is the measured a0 drop between near-contact slices.

PROJECT_ROOT = Path(__file__).resolve().parents[2]
RESULTS_ROOT = PROJECT_ROOT / "results"
RESULTS_ROOT_SWEEPS = PROJECT_ROOT / "results" / "sweeps"
SWEEP_GLOBS = [
    "simple_geometries_diverging_nozzle_gamma_mc_sweep_sweep_*",
    "simple_geometries_diverging_nozzle_gamma_mc_sweep_*",
    "simple_geometries_diverging_nozzle_sweep_*",
]
HARDCODED_SWEEP_DIR = RESULTS_ROOT / "simple_geometries_diverging_nozzle_gamma_mc_sweep_sweep_2026-02-25_142715"

# Pick one gamma_mr slice (for gamma_mc-only sweep this is just the fixed value).
GAMMA_MR_TARGET = 1e-2
GAMMA_MR_TOL = 1e-12

# Sample lines this many y-indices in from the first/last in-domain y-lines.
CONTACT_SAMPLE_DEPTH = 2
A0_BULK_LOW_FRAC = 0.30
A0_BULK_HIGH_FRAC = 0.70

OUT_CSV = Path(__file__).resolve().parent / "diverging_conductance_vs_gamma_mc.csv"
OUT_PNG = Path(__file__).resolve().parent / "diverging_conductance_vs_gamma_mc.png"


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
        raise RuntimeError(
            f"No sweep directories matching {SWEEP_GLOBS} found in {RESULTS_ROOT} or {RESULTS_ROOT_SWEEPS}"
        )
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

        # Identify in-domain y-lines and sample slightly inboard from contacts.
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
    header = "gamma_mc,gamma_mr,i_total,i_norm,delta_a0_measured,conductance\n"
    with OUT_CSV.open("w", encoding="utf-8") as f:
        f.write(header)
        for r in rows:
            f.write(",".join(f"{v:.16g}" for v in r) + "\n")


def main():
    sweep_dir = pick_latest_sweep_dir()
    data_dir = sweep_dir / "data"
    files = sorted(data_dir.glob("*.h5"))
    if not files:
        raise RuntimeError(f"No .h5 files found in {data_dir}")

    rows = parse_rows(files)
    if not rows:
        raise RuntimeError(
            f"No usable rows found for gamma_mr={GAMMA_MR_TARGET} in {data_dir}. "
            "Check GAMMA_MR_TARGET or available sweep files."
        )

    rows = sorted(rows, key=lambda t: t[0])
    write_csv(rows)

    gamma_mc = np.array([r[0] for r in rows])
    conductance = np.array([r[5] for r in rows])

    fig, ax = plt.subplots(figsize=(7.2, 4.6), constrained_layout=True)
    ax.plot(gamma_mc, conductance, "-o", lw=1.8, ms=4.5)
    ax.set_xscale("log")
    ax.set_xlabel(r"$\gamma_{mc}$")
    ax.set_ylabel(r"$G = I_{\mathrm{norm}}/|\Delta a_0|$")
    ax.set_title(r"Diverging Nozzle Conductance vs $\gamma_{mc}$")
    ax.grid(True, which="both", alpha=0.28)
    fig.savefig(OUT_PNG, dpi=220)

    print(f"Using sweep directory: {sweep_dir}")
    print(f"Rows used: {len(rows)}")
    print(f"Saved CSV: {OUT_CSV}")
    print(f"Saved plot: {OUT_PNG}")


if __name__ == "__main__":
    main()
