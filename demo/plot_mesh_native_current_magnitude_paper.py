"""
Render a paper-ready current-magnitude plot from a mesh-native HDF5 analysis file.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
import seaborn as sns


FIGURE_WIDTH_IN = 3.375
FIGURE_HEIGHT_IN = FIGURE_WIDTH_IN * 3.0 / 4.0


plt.rcParams.update(
    {
        "text.usetex": True,
        "font.family": "serif",
        "font.serif": ["Computer Modern Roman", "Times New Roman", "DejaVu Serif"],
        "font.size": 11,
        "axes.labelsize": 11,
        "xtick.labelsize": 11,
        "ytick.labelsize": 11,
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "axes.linewidth": 0.7,
        "xtick.major.width": 0.7,
        "ytick.major.width": 0.7,
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

    return {
        "x": x,
        "y": y,
        "triangles": triangles,
        "speed": np.hypot(jx, jy),
    }


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
        "--cmap",
        type=str,
        default="mako",
        help="Matplotlib/seaborn colormap name",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    fields = load_fields(args.input)

    output = args.output
    if output is None:
        output = args.input.with_suffix("").with_name(args.input.stem + "_paper_current_magnitude.png")

    triang = mtri.Triangulation(fields["x"], fields["y"], fields["triangles"])
    cmap = sns.color_palette(args.cmap, as_cmap=True)

    fig, ax = plt.subplots(figsize=(FIGURE_WIDTH_IN, FIGURE_HEIGHT_IN), constrained_layout=True)
    pcm = ax.tripcolor(
        triang,
        fields["speed"],
        shading="gouraud",
        cmap=cmap,
    )

    ax.set_aspect("equal")
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$y$")

    cbar = fig.colorbar(pcm, ax=ax, pad=0.02, fraction=0.06)
    cbar.set_label(r"$|\mathbf{j}|$")
    cbar.outline.set_linewidth(0.7)

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)
    print(f"saved: {output}")


if __name__ == "__main__":
    main()
