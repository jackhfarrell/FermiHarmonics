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

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman", "Times New Roman", "DejaVu Serif"],
})

PROJECT_ROOT = Path(__file__).resolve().parents[3]
SWEEP_NAME = "simple_geometries_junction_sweep_2026-02-08_231524"
DATA_DIR = PROJECT_ROOT / "projects" / "simple_geometries" / "results" / "results" / SWEEP_NAME / "data"

X_CUT = 0.4
Y_CUT = 0.4
DENOM_EPS = 1e-14

OUT_DIR = Path(__file__).resolve().parent / "figures"
OUT_FIG = OUT_DIR / "junction_i2_over_i3_heatmap.png"
OUT_CSV = OUT_DIR / "junction_i2_i3_ratio_table.csv"


def parse_gamma(path: Path) -> tuple[float, float]:
    stem = path.stem
    gamma_mr = float(re.search(r"gamma_mr=([^_]+)", stem).group(1))
    gamma_mc = float(re.search(r"gamma_mc=([^_]+)", stem).group(1))
    return gamma_mr, gamma_mc


def log_edges(vals: np.ndarray) -> np.ndarray:
    lv = np.log10(vals)
    mids = 0.5 * (lv[:-1] + lv[1:])
    first = lv[0] - (mids[0] - lv[0])
    last = lv[-1] + (lv[-1] - mids[-1])
    return 10.0 ** np.concatenate(([first], mids, [last]))


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


h5_files = sorted(DATA_DIR.glob("*.h5"))
if not h5_files:
    raise FileNotFoundError(f"No h5 files in {DATA_DIR}")

rows: list[tuple[float, float, float, float, float]] = []

for h5_path in h5_files:
    gamma_mr, gamma_mc = parse_gamma(h5_path)
    with h5py.File(h5_path, "r") as f:
        x = np.asarray(f["x"], dtype=float)
        y = np.asarray(f["y"], dtype=float)
        a1 = np.asarray(f["a1"], dtype=float).T
        b1 = np.asarray(f["b1"], dtype=float).T

    iy = int(np.argmin(np.abs(y - Y_CUT)))
    ix = int(np.argmin(np.abs(x - X_CUT)))
    dx = float(x[1] - x[0]) if x.size > 1 else 1.0
    dy = float(y[1] - y[0]) if y.size > 1 else 1.0

    b1_line = b1[:, iy]
    i2 = integrate_simpson(b1_line, dx)

    a1_line = a1[ix, :]
    i3 = integrate_simpson(a1_line, dy)

    ratio = np.nan
    if np.isfinite(i2) and np.isfinite(i3) and abs(i3) > DENOM_EPS:
        ratio = i2 / i3

    rows.append((gamma_mr, gamma_mc, i2, i3, ratio))

gamma_mr_vals = np.array(sorted({r[0] for r in rows}), dtype=float)
gamma_mc_vals = np.array(sorted({r[1] for r in rows}), dtype=float)
mr_to_i = {v: i for i, v in enumerate(gamma_mr_vals)}
mc_to_j = {v: j for j, v in enumerate(gamma_mc_vals)}

ratio_grid = np.full((len(gamma_mc_vals), len(gamma_mr_vals)), np.nan, dtype=float)
for gamma_mr, gamma_mc, _, _, ratio in rows:
    ratio_grid[mc_to_j[gamma_mc], mr_to_i[gamma_mr]] = ratio

finite = np.isfinite(ratio_grid)
if np.count_nonzero(finite) == 0:
    raise RuntimeError("No finite I2/I3 values.")

vmin = float(np.nanpercentile(ratio_grid[finite], 2))
vmax = float(np.nanpercentile(ratio_grid[finite], 98))
if not np.isfinite(vmin):
    vmin = float(np.nanmin(ratio_grid[finite]))
if not np.isfinite(vmax):
    vmax = float(np.nanmax(ratio_grid[finite]))
if vmax <= vmin:
    vmax = vmin + 1e-12

fig, ax = plt.subplots(figsize=(7.2, 5.6))
pc = ax.pcolormesh(
    log_edges(gamma_mr_vals),
    log_edges(gamma_mc_vals),
    np.ma.masked_invalid(ratio_grid),
    shading="auto",
    cmap="viridis",
    vmin=vmin,
    vmax=vmax,
)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel(r"$\gamma_{\mathrm{mr}}$")
ax.set_ylabel(r"$\gamma_{\mathrm{mc}}$")
ax.set_title(r"$I_2/I_3,\ I_2=\int b_1(x,y=0.4)\,dx,\ I_3=\int a_1(x=0.4,y)\,dy$")
fig.colorbar(pc, ax=ax, label=r"$I_2/I_3$")

OUT_DIR.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT_FIG, dpi=300, bbox_inches="tight")
plt.close(fig)

with OUT_CSV.open("w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["gamma_mr", "gamma_mc", "I2", "I3", "I2_over_I3"])
    writer.writerows(sorted(rows, key=lambda r: (r[0], r[1])))

print(f"sweep: {SWEEP_NAME}")
print(f"files: {len(h5_files)}")
print(f"figure: {OUT_FIG}")
print(f"table: {OUT_CSV}")
