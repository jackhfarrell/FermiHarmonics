from __future__ import annotations

import os
import re
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parent / ".mplconfig"))
os.environ.setdefault("XDG_CACHE_HOME", str(Path(__file__).resolve().parent / ".cache"))

import h5py
import matplotlib.pyplot as plt
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[3]
SWEEP_NAME = "simple_geometries_junction_sweep_2026-02-08_231524"
DATA_DIR = PROJECT_ROOT / "projects" / "simple_geometries" / "results" / "results" / SWEEP_NAME / "data"

# Select run by gamma index:
# I_IDX -> gamma_mr index
# J_IDX -> gamma_mc index
I_IDX = 0
J_IDX = 0

OUT_PATH = Path(__file__).resolve().parent / "figures" / "junction_run_by_index.png"


def parse_gamma(path: Path) -> tuple[float, float]:
    stem = path.stem
    m_mr = re.search(r"gamma_mr=([^_]+)", stem)
    m_mc = re.search(r"gamma_mc=([^_]+)", stem)
    if m_mr is None or m_mc is None:
        raise ValueError(f"Could not parse gamma values from {path.name}")
    return float(m_mr.group(1)), float(m_mc.group(1))


h5_files = sorted(DATA_DIR.glob("*.h5"))
if not h5_files:
    raise FileNotFoundError(f"No .h5 files found in {DATA_DIR}")

records: list[tuple[float, float, Path]] = []
for path in h5_files:
    gamma_mr, gamma_mc = parse_gamma(path)
    records.append((gamma_mr, gamma_mc, path))

gamma_mr_vals = sorted({r[0] for r in records})
gamma_mc_vals = sorted({r[1] for r in records})

if not (0 <= I_IDX < len(gamma_mr_vals)):
    raise IndexError(f"I_IDX={I_IDX} out of bounds for gamma_mr grid size {len(gamma_mr_vals)}")
if not (0 <= J_IDX < len(gamma_mc_vals)):
    raise IndexError(f"J_IDX={J_IDX} out of bounds for gamma_mc grid size {len(gamma_mc_vals)}")

gamma_mr_sel = gamma_mr_vals[I_IDX]
gamma_mc_sel = gamma_mc_vals[J_IDX]

path_map = {(gamma_mr, gamma_mc): path for gamma_mr, gamma_mc, path in records}
if (gamma_mr_sel, gamma_mc_sel) not in path_map:
    raise RuntimeError(
        f"No file found for selected indices (I_IDX={I_IDX}, J_IDX={J_IDX}) "
        f"-> (gamma_mr={gamma_mr_sel}, gamma_mc={gamma_mc_sel})"
    )

run_path = path_map[(gamma_mr_sel, gamma_mc_sel)]
with h5py.File(run_path, "r") as f:
    x = np.asarray(f["x"], dtype=float)
    y = np.asarray(f["y"], dtype=float)
    a0 = np.asarray(f["a0"], dtype=float)
    a1 = np.asarray(f["a1"], dtype=float)
    b1 = np.asarray(f["b1"], dtype=float)
    mask = np.asarray(f["mask"]).astype(bool)

a0_plot = np.where(mask, a0, np.nan)
a1_plot = np.where(mask, a1, np.nan)
b1_plot = np.where(mask, b1, np.nan)

fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.2), constrained_layout=True)
fields = [
    (a0_plot, "a0", "viridis"),
    (a1_plot, "a1", "RdBu_r"),
    (b1_plot, "b1", "RdBu_r"),
]

for ax, (field, label, cmap) in zip(axes, fields):
    if label == "a0":
        vmin, vmax = np.nanpercentile(field, [2, 98])
    else:
        vmax = np.nanpercentile(np.abs(field), 98)
        vmin = -vmax
    pcm = ax.pcolormesh(x, y, field, shading="auto", cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_aspect("equal")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(label)
    fig.colorbar(pcm, ax=ax, fraction=0.05, pad=0.02)

fig.suptitle(
    f"junction run by index: I_IDX={I_IDX}, J_IDX={J_IDX}  "
    f"(gamma_mr={gamma_mr_sel:.4g}, gamma_mc={gamma_mc_sel:.4g})"
)

OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT_PATH, dpi=300)
plt.close(fig)

print(f"selected file: {run_path}")
print(f"gamma_mr[{I_IDX}] = {gamma_mr_sel}")
print(f"gamma_mc[{J_IDX}] = {gamma_mc_sel}")
print(f"saved figure: {OUT_PATH}")
