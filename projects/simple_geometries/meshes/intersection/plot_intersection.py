import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.offsetbox import TextArea, HPacker, AnnotationBbox

# Geometry parameters (match intersection.geo)
Wmain = 0.4
Warm = Wmain
L = 1.0

xm = Wmain / 2.0
ya = Warm / 2.0
xside = xm + L

# Polygon points in order (counterclockwise)
pts = [
    (-xm, -L),
    ( xm, -L),
    ( xm, -ya),
    (xside, -ya),
    (xside,  ya),
    ( xm,  ya),
    ( xm,   L),
    (-xm,   L),
    (-xm, -L),  # close
]

# Matplotlib + LaTeX styling (Computer Modern)
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman"],
    "axes.linewidth": 0.0,
})

fig, ax = plt.subplots(figsize=(4.6, 4.6), dpi=300)
fig.patch.set_facecolor("white")
ax.set_facecolor("white")

xs, ys = zip(*pts)

# Light grey interior fill (slightly darker)
ax.fill(xs, ys, color="#e0e0e0", zorder=0)

# Draw boundary (thinner black line)
channel_lw = 1.0
ax.plot(xs, ys, color="black", linewidth=channel_lw, zorder=2)

# Manual pastel primary colors (initial palette)
color_s = "#a3c9f4"  # pastel blue (swapped)
color_n = "#f4a3a3"  # pastel red (swapped)
color_e = "#f6e39a"  # pastel yellow

# Darken label colors for readability while keeping pastel contacts
def darken(color, factor=0.65):
    if isinstance(color, str):
        color = mpl.colors.to_rgb(color)
    return tuple(factor * c for c in color)

text_s = darken(color_s)
text_n = darken(color_n)
text_e = darken(color_e)

label_fs = 13
label_sep = 2


def add_label(ax, x, y, left_text, right_text, left_color, right_color, box_alignment):
    left = TextArea(left_text, textprops=dict(color=left_color, fontsize=label_fs))
    right = TextArea(right_text, textprops=dict(color=right_color, fontsize=label_fs))
    box = HPacker(children=[left, right], align="center", pad=0, sep=label_sep)
    ab = AnnotationBbox(box, (x, y), xycoords="data", frameon=False, box_alignment=box_alignment)
    ax.add_artist(ab)


# Positions for labels near the three contacts
add_label(ax, 0.0, -L - 0.10, r"$V_S$", r"$= V$", text_s, "black", box_alignment=(0.5, 1.0))
add_label(ax, 0.0,  L + 0.10, r"$V_N$", r"$= 0$", text_n, "black", box_alignment=(0.5, 0.0))
add_label(ax, xside + 0.10, 0.0, r"$V_E$", r"$= 0$", text_e, "black", box_alignment=(0.0, 0.5))

# Contacts: thicker, wider, with black outline matching channel width
contact_lw = 3.2
outline_lw = contact_lw + 2 * channel_lw

# Bottom contact
ax.add_line(Line2D([ -xm,  xm], [-L, -L], color="black", linewidth=outline_lw, solid_capstyle="round"))
ax.add_line(Line2D([ -xm,  xm], [-L, -L], color=color_s, linewidth=contact_lw, solid_capstyle="round"))

# Top contact
ax.add_line(Line2D([ -xm,  xm], [ L,  L], color="black", linewidth=outline_lw, solid_capstyle="round"))
ax.add_line(Line2D([ -xm,  xm], [ L,  L], color=color_n, linewidth=contact_lw, solid_capstyle="round"))

# Side contact
ax.add_line(Line2D([xside, xside], [-ya, ya], color="black", linewidth=outline_lw, solid_capstyle="round"))
ax.add_line(Line2D([xside, xside], [-ya, ya], color=color_e, linewidth=contact_lw, solid_capstyle="round"))

# Framing
pad = 0.24
ax.set_xlim(-xm - pad, xside + pad)
ax.set_ylim(-L - pad, L + pad)
ax.set_aspect("equal", adjustable="box")
ax.axis("off")

# Save
out_pdf = "/Users/jfarrell/Desktop/ElectronKinetics.jl/projects/simple_geometries/meshes/intersection/intersection_channel.pdf"
out_png = "/Users/jfarrell/Desktop/ElectronKinetics.jl/projects/simple_geometries/meshes/intersection/intersection_channel.png"
pad_inches = 0.799
fig.savefig(out_pdf, bbox_inches="tight", pad_inches=pad_inches)
fig.savefig(out_png, bbox_inches="tight", pad_inches=pad_inches)
