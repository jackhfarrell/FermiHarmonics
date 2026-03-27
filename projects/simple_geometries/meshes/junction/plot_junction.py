import re
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.offsetbox import TextArea, HPacker, AnnotationBbox

# Matplotlib + LaTeX styling (Computer Modern)
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman"],
    "axes.linewidth": 0.0,
})

POINT_RE = re.compile(r"Point\((\d+)\)\s*=\s*\{([^,]+),\s*([^,]+),")
LINE_RE = re.compile(r"Line\((\d+)\)\s*=\s*\{(\d+),\s*(\d+)\}")
PHYSICAL_RE = re.compile(r'Physical Curve\("([^"]+)"\)\s*=\s*\{([^}]*)\}')
LOOP_RE = re.compile(r"Curve Loop\(\d+\)\s*=\s*\{([^}]*)\}")


def parse_geo(path: Path):
    points = {}
    lines = {}
    physical = {}
    curve_loop = []

    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("//"):
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


def draw_contact(ax, p0, p1, color, contact_lw, outline_lw, length_scale=1.15):
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

    ax.add_line(Line2D([p0e[0], p1e[0]], [p0e[1], p1e[1]], color="black", linewidth=outline_lw, solid_capstyle="round", zorder=4))
    ax.add_line(Line2D([p0e[0], p1e[0]], [p0e[1], p1e[1]], color=color, linewidth=contact_lw, solid_capstyle="round", zorder=5))


def darken(color, factor=0.65):
    if isinstance(color, str):
        color = mpl.colors.to_rgb(color)
    return tuple(factor * c for c in color)


label_sep = 2


def add_label(ax, x, y, left_text, right_text, left_color, right_color, box_alignment, fontsize):
    left = TextArea(left_text, textprops=dict(color=left_color, fontsize=fontsize))
    right = TextArea(right_text, textprops=dict(color=right_color, fontsize=fontsize))
    box = HPacker(children=[left, right], align="center", pad=0, sep=label_sep)
    ab = AnnotationBbox(box, (x, y), xycoords="data", frameon=False, box_alignment=box_alignment)
    ax.add_artist(ab)


def main():
    geo_path = Path(__file__).resolve().with_name("junction.geo")
    points, lines, physical, curve_loop = parse_geo(geo_path)
    poly = polygon_from_loop(points, lines, curve_loop)

    fig, ax = plt.subplots(figsize=(4.8, 4.8), dpi=300)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    xs = [x for x, _ in poly]
    ys = [y for _, y in poly]

    # Light grey interior fill
    ax.fill(xs, ys, color="#e0e0e0", zorder=0)

    # Draw boundary (thin black line)
    channel_lw = 1.0
    ax.plot(xs + [xs[0]], ys + [ys[0]], color="black", linewidth=channel_lw, zorder=2)

    # Pastel palette (matching intersection) + extra green for west
    color_s = "#a3c9f4"  # pastel blue (swapped)
    color_n = "#f4a3a3"  # pastel red (swapped)
    color_e = "#f6e39a"  # pastel yellow
    color_w = "#b7e4b7"  # pastel green

    # Darker text colors for readability
    text_s = darken(color_s)
    text_n = darken(color_n)
    text_e = darken(color_e)
    text_w = darken(color_w)

    # Contacts: thicker with black outline matching channel width
    contact_lw = 3.2
    outline_lw = contact_lw + 2 * channel_lw

    contact_colors = {
        "bottom": color_s,
        "middle": color_n,
        "right": color_e,
        "left": color_w,
    }
    for name, color in contact_colors.items():
        for line_id in physical[name]:
            a, b = lines[line_id]
            draw_contact(ax, points[a], points[b], color, contact_lw, outline_lw)

    # Labels
    label_fs = 13
    add_label(ax, 0.0, min(ys) - 0.08, r"$V_S$", r"$= V$", text_s, "black", box_alignment=(0.5, 1.0), fontsize=label_fs)
    add_label(ax, 0.0, max(ys) + 0.08, r"$V_N$", r"$= 0$", text_n, "black", box_alignment=(0.5, 0.0), fontsize=label_fs)
    add_label(ax, max(xs) + 0.08, 0.0, r"$V_E$", r"$= 0$", text_e, "black", box_alignment=(0.0, 0.5), fontsize=label_fs)
    add_label(ax, min(xs) - 0.08, 0.0, r"$V_W$", r"$= 0$", text_w, "black", box_alignment=(1.0, 0.5), fontsize=label_fs)

    pad = 0.12
    ax.set_xlim(min(xs) - pad, max(xs) + pad)
    ax.set_ylim(min(ys) - pad, max(ys) + pad)
    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")

    out_pdf = Path(__file__).with_name("junction_channel.pdf")
    out_png = Path(__file__).with_name("junction_channel.png")
    fig.savefig(out_pdf, bbox_inches="tight", pad_inches=0.02)
    fig.savefig(out_png, bbox_inches="tight", pad_inches=0.02)


if __name__ == "__main__":
    main()
