from pathlib import Path
import re

import h5py
import matplotlib.pyplot as plt
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
RESULTS_ROOT = PROJECT_ROOT / "results"
RESULTS_ROOT_SWEEPS = PROJECT_ROOT / "results" / "sweeps"
SWEEP_GLOBS = [
    "simple_geometries_dogleg_gamma_mc_sweep_sweep_*",
    "simple_geometries_dogleg_gamma_mc_sweep_*",
    "simple_geometries_dogleg_sweep_*",
]
HARDCODED_SWEEP_DIR = RESULTS_ROOT / "simple_geometries_dogleg_gamma_mc_sweep_sweep_2026-02-25_145852"

GAMMA_MR_TARGET = 1e-2
GAMMA_MR_TOL = 1e-12
CONTACT_SAMPLE_DEPTH = 2

OUT_CSV = Path(__file__).resolve().parent / "dogleg_conductance_vs_gamma_mc.csv"
OUT_PNG = Path(__file__).resolve().parent / "dogleg_conductance_vs_gamma_mc.png"


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


def line_span(coord: np.ndarray, valid: np.ndarray) -> float:
    if np.count_nonzero(valid) < 2:
        return np.nan
    c = coord[valid]
    return float(np.max(c) - np.min(c))


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
            a1 = np.asarray(h5["a1"]).T
            b1 = np.asarray(h5["b1"]).T
            x = np.asarray(h5["x"])
            y = np.asarray(h5["y"])
            mask = np.asarray(h5["mask"]).astype(bool).T

        # Source is near lower-y boundary (horizontal contact): use b1(x) at low y.
        y_valid = np.where(np.any(mask, axis=0))[0]
        x_valid = np.where(np.any(mask, axis=1))[0]
        if len(y_valid) <= CONTACT_SAMPLE_DEPTH or len(x_valid) <= CONTACT_SAMPLE_DEPTH:
            continue
        iy_source = int(y_valid[CONTACT_SAMPLE_DEPTH])
        ix_drain = int(x_valid[-1 - CONTACT_SAMPLE_DEPTH])  # Drain is near max-x boundary (vertical contact).

        valid_source = mask[:, iy_source] & np.isfinite(b1[:, iy_source])
        valid_drain = mask[ix_drain, :] & np.isfinite(a1[ix_drain, :])
        if np.count_nonzero(valid_source) < 2 or np.count_nonzero(valid_drain) < 2:
            continue

        i_source = np.abs(integrate_over_valid_segments(x, b1[:, iy_source], valid_source))
        i_drain = np.abs(integrate_over_valid_segments(y, a1[ix_drain, :], valid_drain))
        i_total = 0.5 * (i_source + i_drain)

        w_source = line_span(x, valid_source)
        w_drain = line_span(y, valid_drain)
        width_eff = 0.5 * (w_source + w_drain)
        i_norm = i_total / width_eff if np.isfinite(width_eff) and width_eff > 0 else np.nan

        # Bulk a0 slices: source-side in vertical leg interior, drain-side in horizontal leg interior.
        x_min = float(np.min(x[x_valid]))
        x_max = float(np.max(x[x_valid]))
        y_min = float(np.min(y[y_valid]))
        y_max = float(np.max(y[y_valid]))
        w_est = line_span(x, valid_source)
        if not np.isfinite(w_est) or w_est <= 0:
            continue

        x_corner = x_min + w_est
        y_corner = y_max - w_est
        y_source_bulk_target = y_min + 0.5 * max(y_corner - y_min, 0.0)
        x_drain_bulk_target = x_corner + 0.5 * max(x_max - x_corner, 0.0)

        iy_source_bulk = int(y_valid[np.argmin(np.abs(y[y_valid] - y_source_bulk_target))])
        ix_drain_bulk = int(x_valid[np.argmin(np.abs(x[x_valid] - x_drain_bulk_target))])

        a0_source = np.nanmean(np.where(mask[:, iy_source_bulk], a0[:, iy_source_bulk], np.nan))
        a0_drain = np.nanmean(np.where(mask[ix_drain_bulk, :], a0[ix_drain_bulk, :], np.nan))
        delta_a0 = np.abs(a0_source - a0_drain)

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
    ax.set_title(r"Dog-Leg Conductance vs $\gamma_{mc}$")
    ax.grid(True, which="both", alpha=0.28)
    fig.savefig(OUT_PNG, dpi=220)

    print(f"Using sweep directory: {sweep_dir}")
    print(f"Rows used: {len(rows)}")
    print(f"Saved CSV: {OUT_CSV}")
    print(f"Saved plot: {OUT_PNG}")


if __name__ == "__main__":
    main()
