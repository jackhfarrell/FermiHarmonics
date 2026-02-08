"""
Plot current magnitude and streamlines for ballistic, hydrodynamic, and diffusive demo outputs.
"""

from __future__ import annotations

from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman", "Times New Roman", "DejaVu Serif"],
    "font.size": 16,
    "axes.labelsize": 18,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
})

def integrated_current_at_x(x: np.ndarray, y: np.ndarray, v: np.ndarray, x_cut: float = 0.4) -> float:
    # Use y-current (b1) integrated along y at fixed x = x_cut.
    ix = int(np.argmin(np.abs(x - x_cut)))
    v_line = v[ix, :]
    valid = np.isfinite(v_line)
    if np.count_nonzero(valid) < 2:
        raise RuntimeError(f"Not enough valid points to integrate current at x={x[ix]:.6g}")
    return float(np.trapz(v_line[valid], y[valid]))


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


data_dir = Path(__file__).resolve().parent / "data"
out_path = Path(__file__).resolve().parent / "streamlines.png"
dpi = 300
regimes = ["diffusive", "hydrodynamic", "ballistic"]
x_cut = 0.4

fig = plt.figure(figsize=(18, 5.6), constrained_layout=False)
gs = fig.add_gridspec(
    1,
    4,
    width_ratios=[1.0, 1.0, 1.0, 0.05],
    wspace=0.18,
    left=0.06,
    right=0.94,
    bottom=0.08,
    top=0.92,
)
axes = [fig.add_subplot(gs[0, i]) for i in range(3)]
cax = fig.add_subplot(gs[0, 3])
cmap_rocket = sns.color_palette("rocket", as_cmap=True)
last_pcm = None

# Determine common axis limits across all regimes
all_xlims, all_ylims = [], []
field_data: list[tuple[str, np.ndarray, np.ndarray, np.ndarray, np.ndarray]] = []
speed_norm_max = 0.0
for regime in regimes:
    file_path = data_dir / f"{regime}.h5"
    if not file_path.is_file():
        raise FileNotFoundError(f"Missing file: {file_path}")
    x, y, u, v = load_fields(file_path)
    all_xlims.extend([x.min(), x.max()])
    all_ylims.extend([y.min(), y.max()])
    line_current = integrated_current_at_x(x, y, v, x_cut=x_cut)
    line_scale = max(abs(line_current), 1e-14)
    u_norm = u / line_scale
    v_norm = v / line_scale
    speed_norm = np.hypot(u_norm, v_norm)
    local_max = np.nanmax(speed_norm)
    if np.isfinite(local_max):
        speed_norm_max = max(speed_norm_max, float(local_max))
    field_data.append((regime, x, y, u_norm, v_norm))

common_xlim = (min(all_xlims), max(all_xlims))
common_ylim = (min(all_ylims), max(all_ylims))
if speed_norm_max <= 0:
    speed_norm_max = 1.0

for ax, (regime, x, y, u_norm, v_norm) in zip(axes, field_data):
    speed = np.hypot(u_norm, v_norm)
    speed_masked = np.ma.masked_invalid(speed)
    u_masked = np.ma.masked_invalid(u_norm)
    v_masked = np.ma.masked_invalid(v_norm)

    pcm = ax.pcolormesh(
        x, y, speed_masked, cmap=cmap_rocket, shading='auto', vmin=0.0, vmax=speed_norm_max
    )
    last_pcm = pcm
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
    ax.set_xlim(common_xlim)
    ax.set_ylim(common_ylim)
    ax.text(
        0.5,
        1.02,
        regime.capitalize(),
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=18,
    )
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_aspect("equal")

out_path.parent.mkdir(parents=True, exist_ok=True)
if last_pcm is not None:
    fig.colorbar(
        last_pcm,
        cax=cax,
        label=r"$|\mathbf{j}|$ (arbitrary units)",
    )
fig.savefig(out_path, dpi=dpi)
print(f"saved: {out_path}")
