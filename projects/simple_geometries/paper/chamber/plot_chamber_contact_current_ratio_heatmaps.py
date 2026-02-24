from pathlib import Path
import re

import h5py
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import seaborn as sns

# Research-style script: edit constants here as needed.
PROJECT_ROOT = Path(__file__).resolve().parents[2]
RESULTS_ROOT = PROJECT_ROOT / "results"
SWEEP_GLOB = "simple_geometries_chamber_sweep_*"
MESH_PATH = PROJECT_ROOT / "meshes" / "chamber" / "chamber.inp"

# Plot window in parameter space.
GAMMA_MR_MIN = 1e-2
GAMMA_MR_MAX = 1e2
GAMMA_MC_MIN = 1e-2
GAMMA_MC_MAX = 1e2

# Small amount of NaN-aware smoothing can help readability on sparse/failed points.
SMOOTHING_PASSES = 0
# Sample this many grid spacings into the bulk from each contact.
BULK_OFFSET_CELLS = 2.0

# Save under this new paper/chamber folder.
OUTFILE = Path(__file__).resolve().parent / "chamber_contact_current_ratio_heatmaps.png"

# Numerator, denominator, title, colormap
RATIO_DEFS = [
    ("drain_left", "contact_source", r"$|I_{DL}| / |I_S|$", "rocket"),
    ("drain_center", "contact_source", r"$|I_{DC}| / |I_S|$", "mako"),
    ("drain_right", "contact_source", r"$|I_{DR}| / |I_S|$", "crest"),
    ("drain_left", "drain_center", r"$|I_{DL}| / |I_{DC}|$", "flare"),
    ("drain_right", "drain_center", r"$|I_{DR}| / |I_{DC}|$", "viridis"),
    ("drain_left", "drain_right", r"$|I_{DL}| / |I_{DR}|$", "cubehelix"),
]


def integrate_over_valid_segments(coord: np.ndarray, values: np.ndarray, valid: np.ndarray) -> float:
    total = 0.0
    n = len(coord)
    i = 0
    while i < n:
        if not valid[i]:
            i += 1
            continue
        j = i
        while j < n and valid[j]:
            j += 1
        if j - i >= 2:
            total += np.trapezoid(values[i:j], coord[i:j])
        i = j
    return float(total)


def bilinear_interp(field: np.ndarray, x: np.ndarray, y: np.ndarray, xq: np.ndarray, yq: np.ndarray) -> np.ndarray:
    # Bilinear interpolation on a regular (x, y) grid where field.shape == (len(x), len(y)).
    xq = np.asarray(xq, dtype=float)
    yq = np.asarray(yq, dtype=float)
    out = np.full_like(xq, np.nan, dtype=float)

    inside = (xq >= x[0]) & (xq <= x[-1]) & (yq >= y[0]) & (yq <= y[-1])
    if not np.any(inside):
        return out

    xv = xq[inside]
    yv = yq[inside]

    ix = np.searchsorted(x, xv, side="right") - 1
    iy = np.searchsorted(y, yv, side="right") - 1
    ix = np.clip(ix, 0, len(x) - 2)
    iy = np.clip(iy, 0, len(y) - 2)

    x0 = x[ix]
    x1 = x[ix + 1]
    y0 = y[iy]
    y1 = y[iy + 1]
    tx = np.where(x1 > x0, (xv - x0) / (x1 - x0), 0.0)
    ty = np.where(y1 > y0, (yv - y0) / (y1 - y0), 0.0)

    f00 = field[ix, iy]
    f10 = field[ix + 1, iy]
    f01 = field[ix, iy + 1]
    f11 = field[ix + 1, iy + 1]

    out_inside = (1.0 - tx) * (1.0 - ty) * f00 + tx * (1.0 - ty) * f10 + (1.0 - tx) * ty * f01 + tx * ty * f11
    out[inside] = out_inside
    return out


def choose_horizontal_sample_y(
    y0: float, x_min: float, x_max: float, x: np.ndarray, y: np.ndarray, mask: np.ndarray
) -> float:
    dx = float(np.mean(np.diff(x)))
    dy = float(np.mean(np.diff(y)))
    n_probe = max(200, int(np.ceil((x_max - x_min) / max(0.25 * dx, 1e-12))))
    xq = np.linspace(x_min, x_max, n_probe)

    offsets = [0.0, 0.5 * dy, -0.5 * dy, 1.0 * dy, -1.0 * dy, 1.5 * dy, -1.5 * dy]
    best_y = y0
    best_count = -1
    for off in offsets:
        y_try = y0 + off
        if not (y[0] <= y_try <= y[-1]):
            continue
        mq = bilinear_interp(mask.astype(float), x, y, xq, np.full_like(xq, y_try))
        count = int(np.count_nonzero(mq > 0.5))
        if count > best_count:
            best_count = count
            best_y = y_try
    return best_y


def choose_vertical_sample_x(
    x0: float, y_min: float, y_max: float, x: np.ndarray, y: np.ndarray, mask: np.ndarray
) -> float:
    dx = float(np.mean(np.diff(x)))
    dy = float(np.mean(np.diff(y)))
    n_probe = max(200, int(np.ceil((y_max - y_min) / max(0.25 * dy, 1e-12))))
    yq = np.linspace(y_min, y_max, n_probe)

    offsets = [0.0, 0.5 * dx, -0.5 * dx, 1.0 * dx, -1.0 * dx, 1.5 * dx, -1.5 * dx]
    best_x = x0
    best_count = -1
    for off in offsets:
        x_try = x0 + off
        if not (x[0] <= x_try <= x[-1]):
            continue
        mq = bilinear_interp(mask.astype(float), x, y, np.full_like(yq, x_try), yq)
        count = int(np.count_nonzero(mq > 0.5))
        if count > best_count:
            best_count = count
            best_x = x_try
    return best_x


def smooth_nan_grid(grid: np.ndarray, passes: int = 1) -> np.ndarray:
    # Simple NaN-aware 3x3 box smoothing.
    out = grid.copy()
    k = np.ones((3, 3), dtype=float)
    for _ in range(max(0, passes)):
        padded = np.pad(out, 1, mode="edge")
        next_out = out.copy()
        for i in range(out.shape[0]):
            for j in range(out.shape[1]):
                w = padded[i:i + 3, j:j + 3]
                m = np.isfinite(w)
                if np.any(m):
                    next_out[i, j] = np.sum(w[m] * k[m]) / np.sum(k[m])
                else:
                    next_out[i, j] = np.nan
        out = next_out
    return out


def parse_mesh_nodes_and_nsets(mesh_path: Path):
    node_coords = {}
    nsets = {}

    lines = mesh_path.read_text(encoding="utf-8").splitlines()
    i = 0
    while i < len(lines):
        raw = lines[i].strip()
        upper = raw.upper()

        if upper.startswith("*NODE"):
            i += 1
            while i < len(lines):
                line = lines[i].strip()
                if not line:
                    i += 1
                    continue
                if line.startswith("*"):
                    break
                parts = [p.strip() for p in line.split(",")]
                if len(parts) >= 3:
                    nid = int(parts[0])
                    x = float(parts[1])
                    y = float(parts[2])
                    node_coords[nid] = (x, y)
                i += 1
            continue

        if upper.startswith("*NSET") and "NSET=" in upper:
            m = re.search(r"NSET\s*=\s*([A-Z0-9_]+)", upper)
            if m is None:
                i += 1
                continue
            name = m.group(1).lower()
            generate = "GENERATE" in upper
            ids = []
            i += 1
            while i < len(lines):
                line = lines[i].strip()
                if not line:
                    i += 1
                    continue
                if line.startswith("*"):
                    break
                nums = [int(x) for x in re.findall(r"-?\d+", line)]
                if generate and len(nums) >= 3:
                    start, stop, step = nums[:3]
                    ids.extend(list(range(start, stop + 1, step)))
                else:
                    ids.extend(nums)
                i += 1
            nsets[name] = sorted(set(ids))
            continue

        i += 1

    return node_coords, nsets


def infer_contact_segment(node_coords, node_ids):
    pts = np.array([node_coords[nid] for nid in node_ids], dtype=float)
    x = pts[:, 0]
    y = pts[:, 1]

    x_span = float(np.max(x) - np.min(x))
    y_span = float(np.max(y) - np.min(y))
    tol = 1e-9

    if y_span <= tol:
        return {
            "orientation": "horizontal",
            "coord": float(np.mean(y)),
            "s_min": float(np.min(x)),
            "s_max": float(np.max(x)),
        }
    if x_span <= tol:
        return {
            "orientation": "vertical",
            "coord": float(np.mean(x)),
            "s_min": float(np.min(y)),
            "s_max": float(np.max(y)),
        }

    raise ValueError(
        f"Contact nodes are not axis-aligned: x_span={x_span:.4g}, y_span={y_span:.4g}."
    )


def integrate_contact_current(segment, x, y, a1, b1, mask):
    if segment["orientation"] == "horizontal":
        dx = float(np.mean(np.diff(x)))
        dy = float(np.mean(np.diff(y)))
        n = max(300, int(np.ceil((segment["s_max"] - segment["s_min"]) / max(0.25 * dx, 1e-12))))
        xq = np.linspace(segment["s_min"], segment["s_max"], n)

        # Move into the bulk: contacts near lower boundary sample upward, near upper boundary sample downward.
        fluid_rows = np.any(mask, axis=0)
        y_fluid = y[fluid_rows]
        y_mid = 0.5 * (float(np.min(y_fluid)) + float(np.max(y_fluid)))
        inward_sign = 1.0 if segment["coord"] <= y_mid else -1.0
        base = inward_sign * BULK_OFFSET_CELLS * dy
        offsets = [base, base + 0.5 * dy, base - 0.5 * dy, base + 1.0 * dy, base - 1.0 * dy]
        best_vals = None
        best_valid = None
        best_count = -1
        best_abs_off = np.inf
        for off in offsets:
            y_try = segment["coord"] + off
            if not (y[0] <= y_try <= y[-1]):
                continue
            vals = bilinear_interp(b1, x, y, xq, np.full_like(xq, y_try))
            valid = np.isfinite(vals)
            count = int(np.count_nonzero(valid))
            if count > best_count or (count == best_count and abs(off) < best_abs_off):
                best_count = count
                best_abs_off = abs(off)
                best_vals = vals
                best_valid = valid

        if best_vals is None or best_count < 2:
            return np.nan
        return np.abs(integrate_over_valid_segments(xq, best_vals, best_valid))

    dx = float(np.mean(np.diff(x)))
    dy = float(np.mean(np.diff(y)))
    n = max(300, int(np.ceil((segment["s_max"] - segment["s_min"]) / max(0.25 * dy, 1e-12))))
    yq = np.linspace(segment["s_min"], segment["s_max"], n)

    # Move into the bulk: contacts near left boundary sample rightward, near right boundary sample leftward.
    fluid_cols = np.any(mask, axis=1)
    x_fluid = x[fluid_cols]
    x_mid = 0.5 * (float(np.min(x_fluid)) + float(np.max(x_fluid)))
    inward_sign = 1.0 if segment["coord"] <= x_mid else -1.0
    base = inward_sign * BULK_OFFSET_CELLS * dx
    offsets = [base, base + 0.5 * dx, base - 0.5 * dx, base + 1.0 * dx, base - 1.0 * dx]
    best_vals = None
    best_valid = None
    best_count = -1
    best_abs_off = np.inf
    for off in offsets:
        x_try = segment["coord"] + off
        if not (x[0] <= x_try <= x[-1]):
            continue
        vals = bilinear_interp(a1, x, y, np.full_like(yq, x_try), yq)
        valid = np.isfinite(vals)
        count = int(np.count_nonzero(valid))
        if count > best_count or (count == best_count and abs(off) < best_abs_off):
            best_count = count
            best_abs_off = abs(off)
            best_vals = vals
            best_valid = valid

    if best_vals is None or best_count < 2:
        return np.nan
    return np.abs(integrate_over_valid_segments(yq, best_vals, best_valid))


def safe_ratio(num, den):
    if not np.isfinite(num) or not np.isfinite(den) or np.isclose(den, 0.0):
        return np.nan
    return num / den


# Find most recent chamber sweep.
sweep_dirs = sorted([p for p in RESULTS_ROOT.glob(SWEEP_GLOB) if p.is_dir()], key=lambda p: p.stat().st_mtime)
if not sweep_dirs:
    raise RuntimeError(f"No sweep directory matching {SWEEP_GLOB!r} found in {RESULTS_ROOT}")

sweep_dir = sweep_dirs[-1]
data_dir = sweep_dir / "data"
files = sorted(data_dir.glob("*.h5"))
if not files:
    raise RuntimeError(f"No .h5 files found in {data_dir}")

print(f"Using sweep directory: {sweep_dir}")
print(f"Found {len(files)} files")

# Pull contact segments directly from mesh NSET definitions.
node_coords, nsets = parse_mesh_nodes_and_nsets(MESH_PATH)
required_contacts = ["contact_source", "drain_left", "drain_center", "drain_right"]
for name in required_contacts:
    if name not in nsets:
        raise RuntimeError(f"Missing NSET {name!r} in {MESH_PATH}")

contact_segments = {name: infer_contact_segment(node_coords, nsets[name]) for name in required_contacts}
for name, seg in contact_segments.items():
    print(
        f"{name}: {seg['orientation']} at {seg['coord']:.4g}, span=({seg['s_min']:.4g}, {seg['s_max']:.4g})"
    )

pattern = re.compile(r"gamma_mc=([0-9eE+\-.]+)_gamma_mr=([0-9eE+\-.]+)")
rows = []

for fpath in files:
    m = pattern.search(fpath.name)
    if m is None:
        continue

    gamma_mc = float(m.group(1))
    gamma_mr = float(m.group(2))

    with h5py.File(fpath, "r") as h5:
        # Match junction paper scripts: transpose so index 0 aligns with x and index 1 aligns with y.
        a1 = np.asarray(h5["a1"]).T
        b1 = np.asarray(h5["b1"]).T
        x = np.asarray(h5["x"])
        y = np.asarray(h5["y"])
        mask = np.asarray(h5["mask"]).astype(bool).T

    if a1.shape != b1.shape or a1.shape != mask.shape:
        raise RuntimeError(
            f"Inconsistent field shapes for {fpath.name}: "
            f"a1={a1.shape}, b1={b1.shape}, mask={mask.shape}."
        )
    if a1.ndim != 2:
        raise RuntimeError(f"Expected 2D fields for {fpath.name}, got a1.ndim={a1.ndim}.")

    currents = {
        name: integrate_contact_current(seg, x, y, a1, b1, mask)
        for name, seg in contact_segments.items()
    }
    rows.append((gamma_mr, gamma_mc, currents))

if not rows:
    raise RuntimeError("No parseable files found (filename regex did not match).")

rows = [r for r in rows if GAMMA_MR_MIN <= r[0] <= GAMMA_MR_MAX and GAMMA_MC_MIN <= r[1] <= GAMMA_MC_MAX]
if not rows:
    raise RuntimeError("No rows remain after gamma_mr/gamma_mc truncation.")

gamma_mr_vals = np.array(sorted({r[0] for r in rows}))
gamma_mc_vals = np.array(sorted({r[1] for r in rows}))

ratio_grids = {}
for num_name, den_name, _, _ in RATIO_DEFS:
    ratio_grids[(num_name, den_name)] = np.full((len(gamma_mr_vals), len(gamma_mc_vals)), np.nan)

for gamma_mr, gamma_mc, currents in rows:
    i = np.searchsorted(gamma_mr_vals, gamma_mr)
    j = np.searchsorted(gamma_mc_vals, gamma_mc)

    for num_name, den_name, _, _ in RATIO_DEFS:
        ratio = safe_ratio(currents[num_name], currents[den_name])
        ratio_grids[(num_name, den_name)][i, j] = ratio

fig, axes = plt.subplots(2, 3, figsize=(12.0, 7.2), constrained_layout=True, sharex=True, sharey=True)
axes = axes.flatten()

for ax, (num_name, den_name, title, palette) in zip(axes, RATIO_DEFS):
    grid = smooth_nan_grid(ratio_grids[(num_name, den_name)], passes=SMOOTHING_PASSES)
    finite = grid[np.isfinite(grid)]
    if finite.size == 0:
        raise RuntimeError(f"No finite values to plot for ratio {num_name}/{den_name}.")

    vmin = float(np.quantile(finite, 0.01))
    vmax = float(np.quantile(finite, 0.99))
    if not np.isfinite(vmin) or not np.isfinite(vmax) or vmax <= vmin:
        vmin, vmax = float(np.nanmin(finite)), float(np.nanmax(finite))
        if not np.isfinite(vmin) or not np.isfinite(vmax) or vmax <= vmin:
            vmin, vmax = 0.0, 1.0

    cmap = sns.color_palette(palette, as_cmap=True).copy()
    cmap.set_bad("0.85")

    z = np.ma.masked_invalid(grid.T)
    pcm = ax.pcolormesh(
        gamma_mr_vals,
        gamma_mc_vals,
        z,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        shading="auto",
    )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(GAMMA_MR_MIN, GAMMA_MR_MAX)
    ax.set_ylim(GAMMA_MC_MIN, GAMMA_MC_MAX)
    ax.set_xlabel(r"$\gamma_{mr}$")
    ax.set_title(title)

    cb = fig.colorbar(pcm, ax=ax)
    cb.formatter = mticker.FormatStrFormatter("%.2g")
    cb.update_ticks()
    cb.set_label(title)

for i in [0, 3]:
    axes[i].set_ylabel(r"$\gamma_{mc}$")

fig.suptitle(r"Chamber contact-current ratio heatmaps")
fig.savefig(OUTFILE, dpi=220)
print(f"Saved figure: {OUTFILE}")
