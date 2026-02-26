"""
Create square-bells illustration assets:
1) geometry (CPHT style, white contacts) from .geo
2) mesh plot (CPHT style) from .inp
3) hydrodynamic streamline plot (from demo/data/hydrodynamic.h5)
"""

from __future__ import annotations

import re
from collections import Counter
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import Polygon
import numpy as np
import seaborn as sns

FIG_BG = "black"
BULK_FACE = "#5A5A5A"
WALL_COLOR = "white"
CONTACT_COLOR = "white"
WALL_LW = 2.0
CONTACT_LW = 15.6
CONTACT_LENGTH_SCALE = 1.35


def parse_geo_geometry(
    geo_path: Path,
) -> tuple[dict[int, tuple[float, float]], dict[int, tuple[int, int]], dict[str, list[int]], list[int]]:
    points: dict[int, tuple[float, float]] = {}
    lines: dict[int, tuple[int, int]] = {}
    physical_curves: dict[str, list[int]] = {}
    curve_loop: list[int] = []

    point_re = re.compile(r"Point\((\d+)\)\s*=\s*\{([^,]+),\s*([^,]+),")
    line_re = re.compile(r"Line\((\d+)\)\s*=\s*\{(\d+),\s*(\d+)\};")
    phys_re = re.compile(r'Physical Curve\("([^"]+)"\)\s*=\s*\{([^}]*)\};')
    loop_re = re.compile(r"Curve Loop\(\d+\)\s*=\s*\{([^}]*)\};")

    for raw in geo_path.read_text(encoding="utf-8").splitlines():
        line = raw.split("//", maxsplit=1)[0].strip()
        if not line:
            continue

        point_match = point_re.search(line)
        if point_match:
            pid = int(point_match.group(1))
            x = float(point_match.group(2))
            y = float(point_match.group(3))
            points[pid] = (x, y)
            continue

        line_match = line_re.search(line)
        if line_match:
            lid = int(line_match.group(1))
            a = int(line_match.group(2))
            b = int(line_match.group(3))
            lines[lid] = (a, b)
            continue

        phys_match = phys_re.search(line)
        if phys_match:
            name = phys_match.group(1)
            ids = [int(tok.strip()) for tok in phys_match.group(2).split(",") if tok.strip()]
            physical_curves[name] = ids
            continue

        loop_match = loop_re.search(line)
        if loop_match:
            curve_loop = [int(tok.strip()) for tok in loop_match.group(1).split(",") if tok.strip()]

    return points, lines, physical_curves, curve_loop


def parse_inp_mesh(inp_path: Path) -> tuple[dict[int, tuple[float, float]], list[tuple[int, int, int, int]], dict[str, set[int]]]:
    nodes: dict[int, tuple[float, float]] = {}
    quads: list[tuple[int, int, int, int]] = []
    nsets: dict[str, set[int]] = {}

    state = ""
    active_nset = ""

    for raw in inp_path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line:
            continue

        if line.startswith("*"):
            upper = line.upper()
            if upper.startswith("*NODE"):
                state = "nodes"
                active_nset = ""
            elif upper.startswith("*ELEMENT") and "TYPE=CPS4" in upper:
                state = "quads"
                active_nset = ""
            elif upper.startswith("*NSET"):
                match = re.search(r"NSET\s*=\s*([^,\s]+)", line, re.IGNORECASE)
                if match:
                    active_nset = match.group(1)
                    nsets.setdefault(active_nset, set())
                    state = "nset"
                else:
                    state = ""
                    active_nset = ""
            else:
                state = ""
                active_nset = ""
            continue

        if state == "nodes":
            parts = [p.strip() for p in line.split(",")]
            if len(parts) >= 3:
                nid = int(parts[0])
                x = float(parts[1])
                y = float(parts[2])
                nodes[nid] = (x, y)
        elif state == "quads":
            parts = [p.strip() for p in line.split(",")]
            if len(parts) >= 5:
                n1 = int(parts[1])
                n2 = int(parts[2])
                n3 = int(parts[3])
                n4 = int(parts[4])
                quads.append((n1, n2, n3, n4))
        elif state == "nset" and active_nset:
            for token in line.split(","):
                token = token.strip()
                if token:
                    nsets[active_nset].add(int(token))

    return nodes, quads, nsets


def edge_segments_from_ids(
    node_coords: dict[int, tuple[float, float]], edges: list[tuple[int, int]]
) -> list[tuple[tuple[float, float], tuple[float, float]]]:
    segments = []
    for a, b in edges:
        if a in node_coords and b in node_coords:
            segments.append((node_coords[a], node_coords[b]))
    return segments


def draw_segment(ax: plt.Axes, p0: tuple[float, float], p1: tuple[float, float], color: str, lw: float) -> None:
    ax.plot(
        [p0[0], p1[0]],
        [p0[1], p1[1]],
        color=color,
        lw=lw,
        solid_capstyle="projecting",
        zorder=5,
    )


def draw_contact_segment(
    ax: plt.Axes,
    p0: tuple[float, float],
    p1: tuple[float, float],
    color: str,
    lw: float,
    length_scale: float = 1.0,
) -> None:
    x0, y0 = p0
    x1, y1 = p1
    dx = x1 - x0
    dy = y1 - y0
    seg_len = (dx * dx + dy * dy) ** 0.5
    if seg_len == 0.0:
        return
    ux = dx / seg_len
    uy = dy / seg_len
    extra = 0.5 * (length_scale - 1.0) * seg_len
    p0e = (x0 - ux * extra, y0 - uy * extra)
    p1e = (x1 + ux * extra, y1 + uy * extra)
    draw_segment(ax, p0e, p1e, color=color, lw=lw)


def polygon_from_loop(
    points: dict[int, tuple[float, float]], lines: dict[int, tuple[int, int]], loop: list[int]
) -> list[tuple[float, float]]:
    verts = []
    current = None
    for line_id in loop:
        reverse = line_id < 0
        a, b = lines[abs(line_id)]
        if reverse:
            a, b = b, a
        if current is None:
            verts.append(points[a])
        elif a != current:
            raise ValueError(f"Curve loop discontinuity at line {line_id}")
        verts.append(points[b])
        current = b

    if verts and verts[0] == verts[-1]:
        verts.pop()
    return verts


def plot_geometry_with_contacts(
    points: dict[int, tuple[float, float]],
    lines: dict[int, tuple[int, int]],
    physical_curves: dict[str, list[int]],
    curve_loop: list[int],
    out_path: Path,
) -> None:
    fig, ax = plt.subplots(figsize=(6, 6), dpi=200)
    fig.patch.set_facecolor(FIG_BG)
    ax.set_facecolor(FIG_BG)

    bulk_xy = polygon_from_loop(points, lines, curve_loop)
    bulk = Polygon(bulk_xy, closed=True, facecolor=BULK_FACE, edgecolor="none", zorder=1)
    ax.add_patch(bulk)

    boundary_x = [x for x, _ in bulk_xy] + [bulk_xy[0][0]]
    boundary_y = [y for _, y in bulk_xy] + [bulk_xy[0][1]]
    ax.plot(
        boundary_x,
        boundary_y,
        color=WALL_COLOR,
        lw=WALL_LW,
        solid_capstyle="round",
        solid_joinstyle="round",
        zorder=4,
    )

    for contact_name in ("contact_bottom", "contact_top"):
        for line_id in physical_curves.get(contact_name, []):
            a, b = lines[line_id]
            draw_contact_segment(
                ax,
                points[a],
                points[b],
                color=CONTACT_COLOR,
                lw=CONTACT_LW,
                length_scale=CONTACT_LENGTH_SCALE,
            )

    xs = [x for x, _ in bulk_xy]
    ys = [y for _, y in bulk_xy]
    pad = 0.08
    ax.set_xlim(min(xs) - pad, max(xs) + pad)
    ax.set_ylim(min(ys) - pad, max(ys) + pad)
    ax.set_aspect("equal")
    ax.axis("off")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        out_path,
        dpi=300,
        facecolor=fig.get_facecolor(),
        bbox_inches="tight",
        pad_inches=0.06,
    )
    plt.close(fig)


def extract_mesh_edges(quads: list[tuple[int, int, int, int]]) -> tuple[list[tuple[int, int]], list[tuple[int, int]]]:
    all_edges: set[tuple[int, int]] = set()
    counts: Counter[tuple[int, int]] = Counter()

    for n1, n2, n3, n4 in quads:
        for a, b in ((n1, n2), (n2, n3), (n3, n4), (n4, n1)):
            edge = (a, b) if a < b else (b, a)
            all_edges.add(edge)
            counts[edge] += 1

    boundary_edges = [edge for edge, count in counts.items() if count == 1]
    return sorted(all_edges), sorted(boundary_edges)


def plot_mesh_with_contacts(
    nodes: dict[int, tuple[float, float]],
    quads: list[tuple[int, int, int, int]],
    nsets: dict[str, set[int]],
    out_path: Path,
) -> None:
    all_edges, boundary_edges = extract_mesh_edges(quads)
    fig, ax = plt.subplots(figsize=(6, 6), dpi=220)
    fig.patch.set_facecolor(FIG_BG)
    ax.set_facecolor(FIG_BG)

    mesh_segments = edge_segments_from_ids(nodes, all_edges)
    ax.add_collection(LineCollection(mesh_segments, colors="white", linewidths=0.45, alpha=0.72, zorder=2))

    for name in ("contact_bottom", "contact_top"):
        if name not in nsets:
            continue
        node_ids = nsets.get(name, set())
        selected = [edge for edge in boundary_edges if edge[0] in node_ids and edge[1] in node_ids]
        if not selected:
            continue
        contact_segments = edge_segments_from_ids(nodes, selected)
        ax.add_collection(LineCollection(contact_segments, colors=CONTACT_COLOR, linewidths=2.6, zorder=4))

    xs = [xy[0] for xy in nodes.values()]
    ys = [xy[1] for xy in nodes.values()]
    xpad = 0.05 * (max(xs) - min(xs))
    ypad = 0.05 * (max(ys) - min(ys))
    ax.set_xlim(min(xs) - xpad, max(xs) + xpad)
    ax.set_ylim(min(ys) - ypad, max(ys) + ypad)
    ax.set_aspect("equal")
    ax.axis("off")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        out_path,
        dpi=300,
        facecolor=fig.get_facecolor(),
        bbox_inches="tight",
        pad_inches=0.05,
    )
    plt.close(fig)


def plot_hydrodynamic_streamlines(h5_path: Path, out_path: Path) -> None:
    with h5py.File(h5_path, "r") as f:
        x = np.asarray(f["x"])
        y = np.asarray(f["y"])
        a1 = np.asarray(f["a1"])
        b1 = np.asarray(f["b1"])
        mask = np.asarray(f["mask"]).astype(bool)

    u = np.where(mask, a1, np.nan)
    v = np.where(mask, b1, np.nan)
    speed = np.hypot(u, v)
    speed_max = float(np.nanmax(speed))
    if not np.isfinite(speed_max) or speed_max <= 0:
        speed_max = 1.0

    cmap = sns.color_palette("rocket", as_cmap=True)
    fig, ax = plt.subplots(figsize=(8.0, 5.2), constrained_layout=True)
    ax.set_facecolor("black")
    fig.patch.set_facecolor("black")

    ax.pcolormesh(
        x,
        y,
        np.ma.masked_invalid(speed),
        shading="auto",
        cmap=cmap,
        vmin=0.0,
        vmax=speed_max,
    )
    ax.streamplot(
        x,
        y,
        np.ma.masked_invalid(u),
        np.ma.masked_invalid(v),
        color="white",
        density=1.0,
        linewidth=1.0,
        arrowsize=1.0,
        minlength=0.2,
        integration_direction="both",
    )
    ax.set_xlim(float(np.min(x)), float(np.max(x)))
    ax.set_ylim(float(np.min(y)), float(np.max(y)))
    ax.set_aspect("equal")
    ax.axis("off")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        out_path,
        dpi=450,
        facecolor=fig.get_facecolor(),
        edgecolor="none",
        bbox_inches="tight",
        pad_inches=0.05,
    )
    plt.close(fig)


def main() -> None:
    script_dir = Path(__file__).resolve().parent
    project_dir = script_dir.parent
    repo_root = project_dir.parent.parent
    mesh_dir = project_dir / "mesh"
    out_dir = project_dir / "figures"

    geo_path = mesh_dir / "square_bells.geo"
    inp_path = mesh_dir / "square_bells.inp"
    hydro_path = repo_root / "demo" / "data" / "hydrodynamic.h5"

    geometry_out = out_dir / "square_bells_geometry_contacts.png"
    mesh_out = out_dir / "square_bells_mesh.png"
    hydro_out = out_dir / "square_bells_hydrodynamic.png"

    points, lines, physical_curves, curve_loop = parse_geo_geometry(geo_path)
    plot_geometry_with_contacts(points, lines, physical_curves, curve_loop, geometry_out)

    nodes, quads, nsets = parse_inp_mesh(inp_path)
    plot_mesh_with_contacts(nodes, quads, nsets, mesh_out)

    plot_hydrodynamic_streamlines(hydro_path, hydro_out)

    print(f"saved: {geometry_out}")
    print(f"saved: {mesh_out}")
    print(f"saved: {hydro_out}")


if __name__ == "__main__":
    main()
