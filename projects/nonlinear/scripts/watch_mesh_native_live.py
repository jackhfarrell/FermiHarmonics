from __future__ import annotations

import argparse
import time
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np


def load_fields(path: Path) -> dict[str, np.ndarray] | None:
    if not path.is_file():
        return None
    with h5py.File(path, "r") as f:
        x = np.asarray(f["x"])
        y = np.asarray(f["y"])
        triangles = np.asarray(f["triangles"], dtype=np.int32)
        if triangles.ndim == 2 and triangles.shape[0] == 3 and triangles.shape[1] != 3:
            triangles = triangles.T
        jx = np.asarray(f["jx"])
        jy = np.asarray(f["jy"])
        t = float(f.attrs.get("time", np.nan))
    return {
        "x": x,
        "y": y,
        "triangles": triangles,
        "jx": jx,
        "jy": jy,
        "speed": np.hypot(jx, jy),
        "time": t,
    }


def build_streamline_grid(fields: dict[str, np.ndarray], grid_size: int):
    triang = mtri.Triangulation(fields["x"], fields["y"], fields["triangles"])
    trifinder = triang.get_trifinder()
    jx_interp = mtri.LinearTriInterpolator(triang, fields["jx"])
    jy_interp = mtri.LinearTriInterpolator(triang, fields["jy"])

    x_min, x_max = np.min(fields["x"]), np.max(fields["x"])
    y_min, y_max = np.min(fields["y"]), np.max(fields["y"])
    x_grid = np.linspace(x_min, x_max, grid_size)
    y_grid = np.linspace(y_min, y_max, grid_size)
    xx, yy = np.meshgrid(x_grid, y_grid)
    inside = trifinder(xx, yy) >= 0
    jx_grid = np.where(inside, np.asarray(jx_interp(xx, yy)), np.nan)
    jy_grid = np.where(inside, np.asarray(jy_interp(xx, yy)), np.nan)
    return triang, x_grid, y_grid, jx_grid, jy_grid


def main() -> None:
    parser = argparse.ArgumentParser(description="Watch a live mesh-native HDF5 file and refresh a plot window.")
    parser.add_argument("input", type=Path)
    parser.add_argument("--poll-seconds", type=float, default=1.0)
    parser.add_argument("--stream-grid", type=int, default=180)
    parser.add_argument("--density", type=float, default=0.7)
    args = parser.parse_args()

    plt.ion()
    fig, ax = plt.subplots(figsize=(7.6, 4.8))
    ax.set_facecolor("#111111")
    ax.text(
        0.5,
        0.5,
        f"Waiting for live data:\n{args.input.name}",
        ha="center",
        va="center",
        color="white",
        transform=ax.transAxes,
    )
    ax.set_axis_off()
    plt.show(block=False)
    last_mtime = None

    while plt.fignum_exists(fig.number):
        if args.input.is_file():
            mtime = args.input.stat().st_mtime
            if last_mtime is None or mtime > last_mtime:
                last_mtime = mtime
                fields = load_fields(args.input)
                if fields is not None:
                    triang, x_grid, y_grid, jx_grid, jy_grid = build_streamline_grid(fields, args.stream_grid)
                    ax.clear()
                    ax.set_axis_on()
                    pcm = ax.tripcolor(triang, fields["speed"], shading="gouraud", cmap="mako")
                    ax.streamplot(
                        x_grid,
                        y_grid,
                        jx_grid,
                        jy_grid,
                        color="white",
                        density=args.density,
                        linewidth=0.8,
                        arrowsize=0.9,
                        minlength=0.15,
                        maxlength=5.0,
                        integration_direction="both",
                    )
                    ax.triplot(triang, color=(1.0, 1.0, 1.0, 0.05), linewidth=0.2)
                    ax.set_aspect("equal")
                    ax.set_xlabel("x")
                    ax.set_ylabel("y")
                    ax.set_title(f"Tesla valve live current magnitude (t={fields['time']:.3f})")
                    if len(fig.axes) == 1:
                        fig.colorbar(pcm, ax=ax, pad=0.02, label="|j|")
                    else:
                        fig.axes[1].cla()
                        fig.colorbar(pcm, cax=fig.axes[1], label="|j|")
                    fig.canvas.draw_idle()
                    fig.canvas.flush_events()
        plt.pause(args.poll_seconds)


if __name__ == "__main__":
    main()
