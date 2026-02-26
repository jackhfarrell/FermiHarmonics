from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Arc, Circle
import numpy as np


def main():
    fig, ax = plt.subplots(figsize=(6, 6), dpi=220)
    fig.patch.set_facecolor("black")
    ax.set_facecolor("black")

    # Circular Fermi surface
    kf = 1.0
    circle = Circle((0.0, 0.0), kf, fill=False, edgecolor="white", linewidth=2.2)
    ax.add_patch(circle)

    # Momentum vector p at angle theta
    theta = np.deg2rad(35.0)
    px, py = kf * np.cos(theta), kf * np.sin(theta)
    ax.arrow(
        0.0,
        0.0,
        px,
        py,
        color="white",
        width=0.01,
        head_width=0.08,
        head_length=0.1,
        length_includes_head=True,
        zorder=5,
    )

    # Reference x-axis
    ax.plot([0, 1.2], [0, 0], color="white", linewidth=1.2, alpha=0.9)

    # Angle arc and labels
    arc = Arc((0.0, 0.0), 0.55, 0.55, angle=0, theta1=0, theta2=np.rad2deg(theta), color="white", linewidth=1.6)
    ax.add_patch(arc)
    ax.text(0.33, 0.13, r"$\theta$", color="white", fontsize=24)
    ax.text(px * 0.72, py * 0.72 + 0.08, r"$\mathbf{p}$", color="white", fontsize=24)
    ax.text(0.0, -1.15, r"Fermi surface", color="white", fontsize=20, ha="center")

    ax.set_aspect("equal")
    ax.set_xlim(-1.35, 1.35)
    ax.set_ylim(-1.35, 1.35)
    ax.axis("off")

    out = Path(__file__).with_name("fermi_surface_cartoon.png")
    fig.savefig(out, dpi=300, facecolor=fig.get_facecolor(), bbox_inches="tight", pad_inches=0.04)
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
