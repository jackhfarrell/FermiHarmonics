from pathlib import Path
import re

import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

from palette import CONTACT_COLORS


POINT_RE = re.compile(r"Point\((\d+)\)\s*=\s*\{([^,]+),\s*([^,]+),")
LINE_RE = re.compile(r"Line\((\d+)\)\s*=\s*\{(\d+),\s*(\d+)\}")
PHYSICAL_RE = re.compile(r'Physical Curve\("([^"]+)"\)\s*=\s*\{([^}]*)\}')
LOOP_RE = re.compile(r"Curve Loop\(\d+\)\s*=\s*\{([^}]*)\}")
ASSIGN_RE = re.compile(r"^([A-Za-z_]\w*)\s*=\s*([^;]+);")


def draw_segment(ax, p0, p1, color, lw):
    ax.plot(
        [p0[0], p1[0]],
        [p0[1], p1[1]],
        color=color,
        lw=lw,
        solid_capstyle="butt",
        zorder=5,
    )


def draw_contact_segment(ax, p0, p1, color, lw, length_scale=1.0):
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


def parse_geo(path: Path):
    points = {}
    lines = {}
    physical = {}
    curve_loop = []
    vars_map = {}

    def eval_expr(expr: str) -> float:
        # Restrict evaluation to arithmetic expressions and previously defined variables.
        return float(eval(expr, {"__builtins__": {}}, vars_map))

    for raw in path.read_text().splitlines():
        line = raw.split("//", 1)[0].strip()
        if not line or line.startswith("//"):
            continue

        m = ASSIGN_RE.search(line)
        if m:
            vars_map[m.group(1)] = eval_expr(m.group(2).strip())
            continue

        m = POINT_RE.search(line)
        if m:
            idx = int(m.group(1))
            points[idx] = (eval_expr(m.group(2).strip()), eval_expr(m.group(3).strip()))
            continue

        m = LINE_RE.search(line)
        if m:
            idx = int(m.group(1))
            lines[idx] = (int(m.group(2)), int(m.group(3)))
            continue

        m = PHYSICAL_RE.search(line)
        if m:
            name = m.group(1)
            ids = [int(x.strip()) for x in m.group(2).split(",") if x.strip()]
            physical[name] = ids
            continue

        m = LOOP_RE.search(line)
        if m:
            curve_loop = [int(x.strip()) for x in m.group(1).split(",") if x.strip()]

    return points, lines, physical, curve_loop


def polygon_from_loop(points, lines, loop):
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


def main():
    repo_root = Path(__file__).resolve().parents[2]
    geo_path = repo_root / "projects/simple_geometries/meshes/intersection/intersection.geo"
    points, lines, physical, curve_loop = parse_geo(geo_path)
    bulk_xy = polygon_from_loop(points, lines, curve_loop)

    fig, ax = plt.subplots(figsize=(6, 6), dpi=200)
    fig.patch.set_facecolor("black")
    ax.set_facecolor("black")

    bulk = Polygon(bulk_xy, closed=True, facecolor="#5A5A5A", edgecolor="none", zorder=1)
    ax.add_patch(bulk)

    wall_lw = 2.0
    contact_lw = 15.6
    contact_length_scale = 1.35

    for line_id in physical["walls"]:
        a, b = lines[line_id]
        draw_segment(ax, points[a], points[b], color="white", lw=wall_lw)

    contact_colors = {
        "contact_side": CONTACT_COLORS["east_drain"],
        "contact_top": CONTACT_COLORS["north_drain"],
        "contact_bottom": "white",
    }
    for contact_name, color in contact_colors.items():
        for line_id in physical[contact_name]:
            a, b = lines[line_id]
            draw_contact_segment(
                ax,
                points[a],
                points[b],
                color=color,
                lw=contact_lw,
                length_scale=contact_length_scale,
            )

    ax.set_aspect("equal")
    xs = [x for x, _ in bulk_xy]
    ys = [y for _, y in bulk_xy]
    pad = 0.08
    ax.set_xlim(min(xs) - pad, max(xs) + pad)
    ax.set_ylim(min(ys) - pad, max(ys) + pad)
    ax.axis("off")

    out_path = Path(__file__).with_name("intersection_illustration.png")
    fig.savefig(out_path, dpi=300, facecolor=fig.get_facecolor(), bbox_inches="tight", pad_inches=0.06)
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
