from pathlib import Path
import re

import h5py
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


POINT_RE = re.compile(r"Point\((\d+)\)\s*=\s*\{([^,]+),\s*([^,]+),")
LINE_RE = re.compile(r"Line\((\d+)\)\s*=\s*\{(\d+),\s*(\d+)\};")
LOOP_RE = re.compile(r"Curve Loop\(\d+\)\s*=\s*\{([^}]*)\};")


def parse_geo_boundary(geo_path: Path):
    points = {}
    lines = {}
    loop = []

    for raw in geo_path.read_text(encoding="utf-8").splitlines():
        line = raw.split("//", 1)[0].strip()
        if not line:
            continue

        m = POINT_RE.search(line)
        if m:
            idx = int(m.group(1))
            points[idx] = (float(m.group(2)), float(m.group(3)))
            continue

        m = LINE_RE.search(line)
        if m:
            idx = int(m.group(1))
            lines[idx] = (int(m.group(2)), int(m.group(3)))
            continue

        m = LOOP_RE.search(line)
        if m:
            loop = [int(x.strip()) for x in m.group(1).split(",") if x.strip()]

    if not loop:
        raise RuntimeError(f"No curve loop found in {geo_path}")

    boundary = []
    current = None
    for line_id in loop:
        reverse = line_id < 0
        a, b = lines[abs(line_id)]
        if reverse:
            a, b = b, a
        if current is None:
            boundary.append(points[a])
        elif a != current:
            raise RuntimeError(f"Boundary loop discontinuity at line {line_id}")
        boundary.append(points[b])
        current = b

    return boundary


def integrated_current_at_x(x: np.ndarray, y: np.ndarray, v: np.ndarray, x_cut: float = 0.4) -> float:
    ix = int(np.argmin(np.abs(x - x_cut)))
    v_line = v[ix, :]
    valid = np.isfinite(v_line)
    if np.count_nonzero(valid) < 2:
        raise RuntimeError(f"Not enough valid points to integrate current at x={x[ix]:.6g}")
    return float(np.trapezoid(v_line[valid], y[valid]))


def main():
    repo_root = Path(__file__).resolve().parents[2]
    h5_path = repo_root / "demo" / "data" / "intermediate.h5"
    geo_path = repo_root / "projects" / "square_bells_ucsb" / "mesh" / "square_bells.geo"
    out_path = Path(__file__).with_name("square_bells_intermediate_cpht.png")

    with h5py.File(h5_path, "r") as f:
        x = np.asarray(f["x"])
        y = np.asarray(f["y"])
        a1 = np.asarray(f["a1"])
        b1 = np.asarray(f["b1"])
        mask = np.asarray(f["mask"]).astype(bool)

    u = np.where(mask, a1, np.nan)
    v = np.where(mask, b1, np.nan)

    line_current = integrated_current_at_x(x, y, v, x_cut=0.4)
    line_scale = max(abs(line_current), 1e-14)
    u_norm = u / line_scale
    v_norm = v / line_scale
    speed = np.hypot(u_norm, v_norm)
    speed_masked = np.ma.masked_invalid(speed)

    boundary = parse_geo_boundary(geo_path)
    bx = [p[0] for p in boundary]
    by = [p[1] for p in boundary]

    cmap = sns.color_palette("rocket", as_cmap=True)

    fig, ax = plt.subplots(figsize=(6.5, 5.0), constrained_layout=False)
    fig.patch.set_facecolor("black")
    ax.set_facecolor("black")

    vmax = float(np.nanmax(speed))
    if not np.isfinite(vmax) or vmax <= 0:
        vmax = 1.0
    ax.pcolormesh(x, y, speed_masked, cmap=cmap, shading="auto", vmin=0.0, vmax=vmax)

    ax.streamplot(
        x,
        y,
        np.ma.masked_invalid(u_norm),
        np.ma.masked_invalid(v_norm),
        color="white",
        density=1.0,
        linewidth=1.0,
        arrowsize=1.0,
        minlength=0.2,
        integration_direction="both",
    )

    ax.plot(
        bx,
        by,
        color="white",
        lw=2.2,
        solid_capstyle="round",
        solid_joinstyle="round",
        zorder=10,
    )

    ax.set_xlim(float(np.min(x)), float(np.max(x)))
    ax.set_ylim(float(np.min(y)), float(np.max(y)))
    ax.set_aspect("equal")
    ax.axis("off")

    fig.savefig(out_path, dpi=600, facecolor=fig.get_facecolor(), edgecolor="none", bbox_inches="tight", pad_inches=0.05)
    print(f"saved: {out_path}")


if __name__ == "__main__":
    main()
