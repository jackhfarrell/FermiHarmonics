#!/usr/bin/env python3
"""
Plot resistance vs temperature for a rectangle T sweep.

Resistance metrics:
- R_bulk_A0 = |A0_bulk_bottom - A0_bulk_top| / |I_total|
- R_total = |bias| / |I_total|

where I_total is the cross-section-integrated b1 current, averaged over a central
bulk y-band.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
import tomllib

import h5py
import matplotlib.pyplot as plt
import numpy as np


def band_mean(y: np.ndarray, values: np.ndarray, y_lo: float, y_hi: float) -> float:
    in_band = (y >= y_lo) & (y <= y_hi) & np.isfinite(values)
    if not np.any(in_band):
        return float("nan")
    return float(np.mean(values[in_band]))


def mean_profile_over_x(field: np.ndarray, mask: np.ndarray) -> np.ndarray:
    nx, ny = field.shape
    profile = np.full(ny, np.nan, dtype=float)
    for j in range(ny):
        valid = mask[:, j] & np.isfinite(field[:, j])
        if np.any(valid):
            profile[j] = float(np.mean(field[valid, j]))
    return profile


def total_current_profile_over_x(x: np.ndarray, b1: np.ndarray, mask: np.ndarray) -> np.ndarray:
    nx, ny = b1.shape
    profile = np.full(ny, np.nan, dtype=float)
    for j in range(ny):
        valid = mask[:, j] & np.isfinite(b1[:, j])
        if np.count_nonzero(valid) < 2:
            continue
        x_row = x[valid]
        j_row = b1[valid, j]
        profile[j] = float(np.trapezoid(j_row, x_row))
    return profile


def parse_summary_rows(sweep_dir: Path) -> list[dict[str, str]]:
    summary_dir = sweep_dir / "summary"
    csv_files = sorted(summary_dir.glob("summary_task_*.csv"))
    if not csv_files:
        all_csv = summary_dir / "summary_all.csv"
        if all_csv.exists():
            csv_files = [all_csv]
    if not csv_files:
        raise FileNotFoundError(f"No summary CSV found in {summary_dir}")

    rows: list[dict[str, str]] = []
    for csv_path in csv_files:
        with csv_path.open("r", newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                if row.get("status", "").strip() == "ok":
                    rows.append(row)
    if not rows:
        raise RuntimeError(f"No successful rows (status=ok) found in {summary_dir}")
    return rows


def read_bias(sweep_dir: Path) -> float:
    metadata_path = sweep_dir / "sweep_metadata.toml"
    if not metadata_path.exists():
        return 1.0
    with metadata_path.open("rb") as f:
        data = tomllib.load(f)
    return float(data.get("sweep_metadata", {}).get("bias", 1.0))


def compute_metrics(sweep_dir: Path, rows: list[dict[str, str]], bias: float) -> list[dict[str, float | str]]:
    metrics: list[dict[str, float | str]] = []

    for row in rows:
        t_val = float(row["T"])
        gamma_mc = float(row["gamma_mc"])
        gamma_mr = float(row["gamma_mr"])
        rel_output = row["output_file"]
        h5_path = sweep_dir / rel_output

        with h5py.File(h5_path, "r") as h5:
            x = np.asarray(h5["x"], dtype=float)
            y = np.asarray(h5["y"], dtype=float)
            a0 = np.asarray(h5["a0"], dtype=float)
            b1 = np.asarray(h5["b1"], dtype=float)
            mask = np.asarray(h5["mask"]).astype(bool)

        a0_y = mean_profile_over_x(a0, mask)
        i_y = total_current_profile_over_x(x, b1, mask)

        y_min = float(np.min(y))
        y_max = float(np.max(y))
        length = y_max - y_min

        # Match the bulk-band convention used in rectangle contact-resistance analysis.
        a0_bottom_inner = band_mean(y, a0_y, y_min + 0.30 * length, y_min + 0.40 * length)
        a0_top_inner = band_mean(y, a0_y, y_max - 0.40 * length, y_max - 0.30 * length)
        i_total_bulk = band_mean(y, i_y, y_min + 0.45 * length, y_min + 0.55 * length)

        delta_a0_bulk = abs(a0_bottom_inner - a0_top_inner)
        i_abs = abs(i_total_bulk) if math.isfinite(i_total_bulk) else float("nan")

        if i_abs > 0:
            r_bulk_a0 = delta_a0_bulk / i_abs
            r_total = abs(bias) / i_abs
        else:
            r_bulk_a0 = float("nan")
            r_total = float("nan")

        metrics.append(
            {
                "T": t_val,
                "gamma_mc": gamma_mc,
                "gamma_mr": gamma_mr,
                "I_total_bulk": i_total_bulk,
                "A0_bulk_bottom": a0_bottom_inner,
                "A0_bulk_top": a0_top_inner,
                "DeltaA0_bulk": delta_a0_bulk,
                "R_bulk_A0": r_bulk_a0,
                "R_total": r_total,
                "output_file": rel_output,
            }
        )

    metrics.sort(key=lambda r: float(r["T"]))
    return metrics


def write_metrics_csv(metrics: list[dict[str, float | str]], out_csv: Path) -> None:
    fieldnames = [
        "T",
        "gamma_mc",
        "gamma_mr",
        "I_total_bulk",
        "A0_bulk_bottom",
        "A0_bulk_top",
        "DeltaA0_bulk",
        "R_bulk_A0",
        "R_total",
        "output_file",
    ]
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(metrics)


def make_plot(metrics: list[dict[str, float | str]], bias: float, out_png: Path) -> None:
    t = np.asarray([float(r["T"]) for r in metrics], dtype=float)
    r_bulk = np.asarray([float(r["R_bulk_A0"]) for r in metrics], dtype=float)
    r_total = np.asarray([float(r["R_total"]) for r in metrics], dtype=float)

    valid_bulk = np.isfinite(t) & np.isfinite(r_bulk) & (t > 0) & (r_bulk > 0)
    valid_total = np.isfinite(t) & np.isfinite(r_total) & (t > 0) & (r_total > 0)

    fig, ax = plt.subplots(figsize=(8.2, 5.2))
    if np.any(valid_bulk):
        ax.loglog(t[valid_bulk], r_bulk[valid_bulk], "o-", lw=1.8, ms=4, label=r"$R_{\mathrm{bulk},A_0}=|\Delta A_0|/|I|$")
    if np.any(valid_total):
        ax.loglog(t[valid_total], r_total[valid_total], "s-", lw=1.8, ms=4, label=r"$R_{\mathrm{total}}=|\mathrm{bias}|/|I|$")

    ax.set_xlabel("Temperature T")
    ax.set_ylabel("Resistance")
    ax.set_title(f"Rectangle T Sweep Resistance (bias={bias:g})")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=180)
    plt.close(fig)


def default_sweep_dir(results_root: Path) -> Path:
    candidates = sorted(
        [p for p in results_root.glob("*rectangle_T_sweep*") if p.is_dir()],
        key=lambda p: p.stat().st_mtime,
    )
    if not candidates:
        raise FileNotFoundError(f"No rectangle T sweep directories found under {results_root}")
    return candidates[-1]


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot resistance vs T for rectangle T sweep output.")
    parser.add_argument(
        "--sweep-dir",
        type=Path,
        default=None,
        help="Sweep directory path. Default: latest *rectangle_T_sweep* in projects/simple_geometries/results",
    )
    parser.add_argument("--out-png", type=Path, default=None, help="Output plot path (.png)")
    parser.add_argument("--out-csv", type=Path, default=None, help="Output metrics CSV path")
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[4]
    results_root = repo_root / "projects" / "simple_geometries" / "results"
    sweep_dir = args.sweep_dir or default_sweep_dir(results_root)

    out_png = args.out_png or (sweep_dir / "summary" / "resistance_vs_T_bulkA0_totalCurrent.png")
    out_csv = args.out_csv or (sweep_dir / "summary" / "resistance_vs_T_bulkA0_totalCurrent.csv")

    rows = parse_summary_rows(sweep_dir)
    bias = read_bias(sweep_dir)
    metrics = compute_metrics(sweep_dir, rows, bias)
    write_metrics_csv(metrics, out_csv)
    make_plot(metrics, bias, out_png)

    print(f"sweep_dir={sweep_dir}")
    print(f"points={len(metrics)}")
    print(f"csv={out_csv}")
    print(f"png={out_png}")


if __name__ == "__main__":
    main()

