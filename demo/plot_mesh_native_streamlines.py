"""
Render a smooth curved-geometry current plot from a mesh-native HDF5 analysis file.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
import seaborn as sns


plt.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": ["Computer Modern Roman", "Times New Roman", "DejaVu Serif"],
        "font.size": 14,
        "axes.labelsize": 16,
        "axes.titlesize": 18,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "figure.dpi": 150,
        "savefig.dpi": 300,
    }
)


def load_fields(path: Path) -> dict[str, np.ndarray]:
    with h5py.File(path, "r") as f:
        x = np.asarray(f["x"])
        y = np.asarray(f["y"])
        triangles = np.asarray(f["triangles"], dtype=np.int32)
        if triangles.ndim == 2 and triangles.shape[0] == 3 and triangles.shape[1] != 3:
            triangles = triangles.T
        jx = np.asarray(f["jx"])
        jy = np.asarray(f["jy"])
        density = np.asarray(f["n"]) if "n" in f else np.asarray(f["a0"])
        time = float(f.attrs.get("time", np.nan))

    speed = np.hypot(jx, jy)
    return {
        "x": x,
        "y": y,
        "triangles": triangles,
        "jx": jx,
        "jy": jy,
        "speed": speed,
        "density": density,
        "time": time,
    }


def build_streamline_grid(fields: dict[str, np.ndarray], grid_size: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray] | None:
    triang = mtri.Triangulation(fields["x"], fields["y"], fields["triangles"])
    try:
        trifinder = triang.get_trifinder()
        jx_interp = mtri.LinearTriInterpolator(triang, fields["jx"])
        jy_interp = mtri.LinearTriInterpolator(triang, fields["jy"])
    except RuntimeError:
        return None

    x_min, x_max = np.min(fields["x"]), np.max(fields["x"])
    y_min, y_max = np.min(fields["y"]), np.max(fields["y"])
    x_grid = np.linspace(x_min, x_max, grid_size)
    y_grid = np.linspace(y_min, y_max, grid_size)
    xx, yy = np.meshgrid(x_grid, y_grid)

    inside = trifinder(xx, yy) >= 0
    jx_grid = np.where(inside, np.asarray(jx_interp(xx, yy)), np.nan)
    jy_grid = np.where(inside, np.asarray(jy_interp(xx, yy)), np.nan)
    return x_grid, y_grid, jx_grid, jy_grid


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="Input mesh-native HDF5 analysis file")
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output image path (defaults next to the input file)",
    )
    parser.add_argument(
        "--density",
        type=float,
        default=1.1,
        help="Streamline density passed to Matplotlib",
    )
    parser.add_argument(
        "--stream-grid",
        type=int,
        default=250,
        help="Auxiliary regular-grid resolution used only for streamline tracing",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    fields = load_fields(args.input)

    output = args.output
    if output is None:
        output = args.input.with_suffix("").with_name(args.input.stem + "_mesh_native_streamlines.png")

    triang = mtri.Triangulation(fields["x"], fields["y"], fields["triangles"])
    streamline_grid = build_streamline_grid(fields, args.stream_grid)

    fig, ax = plt.subplots(figsize=(7.4, 4.4), constrained_layout=True)
    cmap = sns.color_palette("mako", as_cmap=True)
    pcm = ax.tripcolor(
        triang,
        fields["speed"],
        shading="gouraud",
        cmap=cmap,
    )

    if streamline_grid is not None:
        x_grid, y_grid, jx_grid, jy_grid = streamline_grid
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
    else:
        stride = max(1, len(fields["x"]) // 1200)
        ax.quiver(
            fields["x"][::stride],
            fields["y"][::stride],
            fields["jx"][::stride],
            fields["jy"][::stride],
            color="white",
            alpha=0.8,
            scale_units="xy",
            scale=None,
            width=0.0025,
        )

    ax.triplot(triang, color=(1.0, 1.0, 1.0, 0.05), linewidth=0.2)
    ax.set_aspect("equal")
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$y$")
    ax.set_title(rf"Mesh-native current magnitude and streamlines ($t={fields['time']:.3f}$)")

    cbar = fig.colorbar(pcm, ax=ax, pad=0.02)
    cbar.set_label(r"$|\mathbf{j}|$")

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)
    print(f"saved: {output}")


if __name__ == "__main__":
    main()
