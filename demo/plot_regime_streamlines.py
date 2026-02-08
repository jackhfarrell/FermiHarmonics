#!/usr/bin/env python3
"""Plot current magnitude and streamlines for ballistic, hydrodynamic, and diffusive demo outputs."""

from __future__ import annotations

from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


def load_fields(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    with h5py.File(path, "r") as f:
        x = np.asarray(f["x"])
        y = np.asarray(f["y"])
        a1 = np.asarray(f["a1"])  # stored as (nx, ny)
        b1 = np.asarray(f["b1"])  # stored as (nx, ny)
        mask = np.asarray(f["mask"]).astype(bool)  # stored as (nx, ny)

    u = np.where(mask, a1, np.nan)
    v = np.where(mask, b1, np.nan)
    return x, y, u, v


def integrated_current_at_x(x: np.ndarray, y: np.ndarray, v: np.ndarray, x_cut: float = 0.4) -> float:
    # Use the y-current (b1) integrated along y at fixed x = x_cut.
    ix = int(np.argmin(np.abs(x - x_cut)))
    v_line = v[ix, :]
    valid = np.isfinite(v_line)
    if np.count_nonzero(valid) < 2:
        raise RuntimeError(f"Not enough valid points to integrate current at x={x[ix]:.6g}")
    return float(np.trapz(v_line[valid], y[valid]))


data_dir = Path(__file__).resolve().parent / "data"
out_path = Path(__file__).resolve().parent / "streamlines.png"
dpi = 300
regimes = ["diffusive", "hydrodynamic", "ballistic"]
x_cut = 0.4

fig, axes = plt.subplots(1, 3, figsize=(15, 4.6), constrained_layout=True)
cmap_icefire = sns.color_palette("flare_r", as_cmap=True)

# Determine common axis limits across all regimes
all_xlims, all_ylims = [], []
for regime in regimes:
    file_path = data_dir / f"{regime}.h5"
    if not file_path.is_file():
        raise FileNotFoundError(f"Missing file: {file_path}")
    x, y, _, _ = load_fields(file_path)
    all_xlims.extend([x.min(), x.max()])
    all_ylims.extend([y.min(), y.max()])

common_xlim = (min(all_xlims), max(all_xlims))
common_ylim = (min(all_ylims), max(all_ylims))

for ax, regime in zip(axes, regimes):
    file_path = data_dir / f"{regime}.h5"
    x, y, u, v = load_fields(file_path)
    line_current = integrated_current_at_x(x, y, v, x_cut=x_cut)
    line_scale = max(abs(line_current), 1e-14)
    u_norm = u / line_scale
    v_norm = v / line_scale
    speed = np.hypot(u_norm, v_norm)
    speed_masked = np.ma.masked_invalid(speed)
    u_masked = np.ma.masked_invalid(u_norm)
    v_masked = np.ma.masked_invalid(v_norm)

    pcm = ax.pcolormesh(x, y, speed_masked, cmap=cmap_icefire, shading='auto')
    ax.streamplot(
        x,
        y,
        u_masked,
        v_masked,
        color='white',
        density=1.0,
        linewidth=1.0,
        arrowsize=1.0,
        minlength=0.2,
        integration_direction='both',
    )
    fig.colorbar(pcm, ax=ax, pad=0.02, shrink=0.9, label=r"$|\mathbf{j}| / |I_{x=0.4}|$")

    ax.set_xlim(common_xlim)
    ax.set_ylim(common_ylim)
    ax.set_title(regime)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_aspect("equal")

fig.suptitle("FermiHarmonics Demo Streamlines (a1, b1)")
out_path.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(out_path, dpi=dpi)
print(f"saved: {out_path}")
