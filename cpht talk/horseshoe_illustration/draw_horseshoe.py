from pathlib import Path
import math
import re

import matplotlib.pyplot as plt
from matplotlib.patches import Polygon


POINT_RE = re.compile(r"Point\((\d+)\)\s*=\s*\{([^,]+),\s*([^,]+),")
LINE_RE = re.compile(r"Line\((\d+)\)\s*=\s*\{(\d+),\s*(\d+)\}")
CIRCLE_RE = re.compile(r"Circle\((\d+)\)\s*=\s*\{(\d+),\s*(\d+),\s*(\d+)\}")
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
    curves = {}
    physical = {}
    curve_loop = []
    vars_map = {"Pi": math.pi}

    def eval_expr(expr: str) -> float:
        safe_globals = {"__builtins__": {}, "sqrt": math.sqrt}
        return float(eval(expr, safe_globals, vars_map))

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
            curves[idx] = {"kind": "line", "start": int(m.group(2)), "end": int(m.group(3))}
            continue

        m = CIRCLE_RE.search(line)
        if m:
            idx = int(m.group(1))
            curves[idx] = {
                "kind": "circle",
                "start": int(m.group(2)),
                "center": int(m.group(3)),
                "end": int(m.group(4)),
            }
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

    return points, curves, physical, curve_loop


def oriented_node_ids(curve, reverse=False):
    start = curve["start"]
    end = curve["end"]
    if reverse:
        start, end = end, start
    return start, end


def sample_curve(points, curve, reverse=False, arc_samples=120):
    start_id, end_id = oriented_node_ids(curve, reverse=reverse)

    if curve["kind"] == "line":
        return [points[start_id], points[end_id]]

    center = points[curve["center"]]
    x0, y0 = points[start_id]
    x1, y1 = points[end_id]
    cx, cy = center

    r = math.hypot(x0 - cx, y0 - cy)
    t0 = math.atan2(y0 - cy, x0 - cx)
    t1 = math.atan2(y1 - cy, x1 - cx)

    dt = t1 - t0
    while dt < -math.pi:
        dt += 2 * math.pi
    while dt > math.pi:
        dt -= 2 * math.pi

    n = max(2, int(arc_samples * abs(dt) / math.pi) + 1)
    out = []
    for i in range(n):
        f = i / (n - 1)
        t = t0 + f * dt
        out.append((cx + r * math.cos(t), cy + r * math.sin(t)))
    return out


def draw_curve(ax, points, curve, reverse, color, lw):
    pts = sample_curve(points, curve, reverse=reverse)
    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]
    ax.plot(xs, ys, color=color, lw=lw, solid_capstyle="butt", zorder=5)


def draw_contact_curve(ax, points, curve, reverse, color, lw, length_scale=1.0):
    if curve["kind"] == "line":
        start_id, end_id = oriented_node_ids(curve, reverse=reverse)
        draw_contact_segment(ax, points[start_id], points[end_id], color, lw, length_scale)
        return

    draw_curve(ax, points, curve, reverse=reverse, color=color, lw=lw)


def polygon_from_loop(points, curves, loop):
    verts = []
    current_node = None

    for curve_id in loop:
        reverse = curve_id < 0
        curve = curves[abs(curve_id)]
        start_id, end_id = oriented_node_ids(curve, reverse=reverse)

        if current_node is None:
            pass
        elif start_id != current_node:
            raise ValueError(f"Curve loop discontinuity at curve {curve_id}")

        sampled = sample_curve(points, curve, reverse=reverse)
        if not verts:
            verts.extend(sampled)
        else:
            verts.extend(sampled[1:])

        current_node = end_id

    if verts and verts[0] == verts[-1]:
        verts.pop()
    return verts


def main():
    repo_root = Path(__file__).resolve().parents[2]
    geo_path = repo_root / "projects/simple_geometries/meshes/snake/horseshoe.geo"
    points, curves, physical, curve_loop = parse_geo(geo_path)
    bulk_xy = polygon_from_loop(points, curves, curve_loop)

    fig, ax = plt.subplots(figsize=(6, 6), dpi=200)
    fig.patch.set_facecolor("black")
    ax.set_facecolor("black")

    bulk = Polygon(bulk_xy, closed=True, facecolor="#5A5A5A", edgecolor="none", zorder=1)
    ax.add_patch(bulk)

    wall_lw = 2.0
    contact_lw = 15.6
    contact_length_scale = 1.35

    for curve_id in physical["walls"]:
        reverse = curve_id < 0
        curve = curves[abs(curve_id)]
        draw_curve(ax, points, curve, reverse=reverse, color="white", lw=wall_lw)

    contact_colors = {
        "contact_in": "white",
        "contact_out": "white",
    }
    for contact_name, color in contact_colors.items():
        for curve_id in physical[contact_name]:
            reverse = curve_id < 0
            curve = curves[abs(curve_id)]
            draw_contact_curve(
                ax,
                points,
                curve,
                reverse=reverse,
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

    out_path = Path(__file__).with_name("horseshoe_illustration.png")
    fig.savefig(out_path, dpi=300, facecolor=fig.get_facecolor(), bbox_inches="tight", pad_inches=0.06)
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
