from __future__ import annotations

import csv
import os
import re
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parent / ".mplconfig"))
os.environ.setdefault("XDG_CACHE_HOME", str(Path(__file__).resolve().parent / ".cache"))

import h5py
import matplotlib.pyplot as plt
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[3]
SWEEP_NAME = "simple_geometries_junction_gamma_mc_sweep_sweep_2026-02-11_144931"
DATA_DIR = PROJECT_ROOT / "projects" / "simple_geometries" / "results" / SWEEP_NAME / "data"

X_CUT = 0.4
Y_CUT = 0.4
DENOM_EPS = 1e-14

OUT_DIR = Path(__file__).resolve().parent / "figures"
OUT_FIG = OUT_DIR / "junction_i1_over_i2_gamma_mc_sweep.png"
OUT_CSV = OUT_DIR / "junction_i1_i2_gamma_mc_table.csv"


def parse_gamma_mc(path: Path) -> float:
    match = re.search(r"gamma_mc=([^_]+)", path.stem)
    if match is None:
        raise ValueError(f"Could not parse gamma_mc from {path.name}")
    return float(match.group(1))


def integrate_simpson(values: np.ndarray, h: float) -> float:
    vals = np.array(values, dtype=float, copy=True)
    vals[~np.isfinite(vals)] = 0.0
    n = vals.size
    if n < 2:
        return 0.0
    if n == 2:
        return 0.5 * h * (vals[0] + vals[1])
    if n % 2 == 1:
        s = vals[0] + vals[-1]
        s += 4.0 * np.sum(vals[1:-1:2])
        s += 2.0 * np.sum(vals[2:-2:2])
        return (h / 3.0) * s
    s = vals[0] + vals[-2]
    s += 4.0 * np.sum(vals[1:-2:2])
    s += 2.0 * np.sum(vals[2:-3:2])
    simpson_part = (h / 3.0) * s
    trap_part = 0.5 * h * (vals[-2] + vals[-1])
    return simpson_part + trap_part


if not DATA_DIR.is_dir():
    alt_dir = PROJECT_ROOT / "projects" / "simple_geometries" / "results" / "results" / SWEEP_NAME / "data"
    if alt_dir.is_dir():
        DATA_DIR = alt_dir
    else:
        raise FileNotFoundError(f"No data directory found for sweep {SWEEP_NAME}")

h5_files = sorted(DATA_DIR.glob("*.h5"))
if not h5_files:
    raise FileNotFoundError(f"No h5 files in {DATA_DIR}")

rows: list[tuple[float, float, float, float]] = []

for h5_path in h5_files:
    gamma_mc = parse_gamma_mc(h5_path)
    with h5py.File(h5_path, "r") as f:
        x = np.asarray(f["x"], dtype=float)
        y = np.asarray(f["y"], dtype=float)
        a1 = np.asarray(f["a1"], dtype=float).T
        b1 = np.asarray(f["b1"], dtype=float).T

    ix_left = int(np.argmin(np.abs(x + X_CUT)))
    iy = int(np.argmin(np.abs(y - Y_CUT)))
    dx = float(x[1] - x[0]) if x.size > 1 else 1.0
    dy = float(y[1] - y[0]) if y.size > 1 else 1.0

    a1_line = a1[ix_left, :]
    b1_line = b1[:, iy]

    i1 = integrate_simpson(a1_line, dy)
    i2 = integrate_simpson(b1_line, dx)

    ratio = np.nan
    if np.isfinite(i1) and np.isfinite(i2) and abs(i2) > DENOM_EPS:
        ratio = i1 / i2

    rows.append((gamma_mc, i1, i2, ratio))

rows_sorted = sorted(rows, key=lambda r: r[0])
gamma_mc_vals = np.array([r[0] for r in rows_sorted], dtype=float)
ratios = np.array([r[3] for r in rows_sorted], dtype=float)

finite = np.isfinite(ratios)
if np.count_nonzero(finite) == 0:
    raise RuntimeError("No finite I1/I2 values.")

fig, ax = plt.subplots(figsize=(6.8, 4.6))
ax.plot(gamma_mc_vals[finite], ratios[finite], marker="o", lw=1.5, ms=4)
ax.set_xscale("log")
ax.set_xlabel(r"$\gamma_{\mathrm{mc}}$")
ax.set_ylabel(r"$I_1/I_2$")
ax.set_title(r"$I_1/I_2$ vs $\gamma_{\mathrm{mc}}$")
ax.grid(True, alpha=0.3)

OUT_DIR.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT_FIG, dpi=300, bbox_inches="tight")
plt.close(fig)

with OUT_CSV.open("w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["gamma_mc", "I1", "I2", "I1_over_I2"])
    writer.writerows(rows_sorted)

print(f"sweep: {SWEEP_NAME}")
print(f"files: {len(h5_files)}")
print(f"figure: {OUT_FIG}")
print(f"table: {OUT_CSV}")
