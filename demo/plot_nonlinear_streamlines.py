"""
Render a publication-style nonlinear current plot from an HDF5 analysis file.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


plt.rcParams.update(
    {
        "text.usetex": True,
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
        jx = np.asarray(f["jx"])
        jy = np.asarray(f["jy"])
        mask = np.asarray(f["mask"]).astype(bool)
        density = np.asarray(f["n"])
        time = float(f.attrs.get("time", np.nan))

    jx = np.where(mask, jx, np.nan)
    jy = np.where(mask, jy, np.nan)
    density = np.where(mask, density, np.nan)
    speed = np.hypot(jx, jy)

    return {
        "x": x,
        "y": y,
        "jx": jx,
        "jy": jy,
        "speed": speed,
        "density": density,
        "mask": mask,
        "time": time,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="Input HDF5 analysis file")
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output image path (defaults next to the input file)",
    )
    parser.add_argument(
        "--density",
        type=float,
        default=1.25,
        help="Streamline density passed to Matplotlib",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    fields = load_fields(args.input)

    output = args.output
    if output is None:
        output = args.input.with_suffix("").with_name(args.input.stem + "_streamlines.png")

    x = fields["x"]
    y = fields["y"]
    jx = np.ma.masked_invalid(fields["jx"])
    jy = np.ma.masked_invalid(fields["jy"])
    speed = np.ma.masked_invalid(fields["speed"])

    fig, ax = plt.subplots(figsize=(7.4, 4.4), constrained_layout=True)
    cmap = sns.color_palette("mako", as_cmap=True)

    pcm = ax.pcolormesh(
        x,
        y,
        speed,
        shading="auto",
        cmap=cmap,
    )
    ax.streamplot(
        x,
        y,
        jx,
        jy,
        color="white",
        density=args.density,
        linewidth=0.8,
        arrowsize=0.9,
        minlength=0.15,
        maxlength=5.0,
        integration_direction="both",
    )

    ax.set_aspect("equal")
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$y$")
    ax.set_title(rf"Nonlinear current magnitude and streamlines ($t={fields['time']:.3f}$)")

    cbar = fig.colorbar(pcm, ax=ax, pad=0.02)
    cbar.set_label(r"$|\mathbf{j}|$")

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)
    print(f"saved: {output}")


if __name__ == "__main__":
    main()
