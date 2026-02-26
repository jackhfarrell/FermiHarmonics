"""
Plot |I_E| / |I_N| versus gamma_mc for selected gamma_mr curves
from an intersection full-sweep directory.

Defaults to the most recent sweep under:
  projects/simple_geometries/full_sweeps/simple_geometries_intersection_sweep_*
"""

from __future__ import annotations

import argparse
from pathlib import Path
import re

import h5py
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.cm import ScalarMappable
import numpy as np
import seaborn as sns

plt.rcParams.update(
    {
        "text.usetex": True,
        "font.family": "serif",
        "font.serif": ["Computer Modern Roman", "Times New Roman", "DejaVu Serif"],
        "font.size": 16,
        "axes.labelsize": 20,
        "xtick.labelsize": 15,
        "ytick.labelsize": 15,
    }
)

# Match intersection contact colors used in CPHT illustrations.
EAST_COLOR = "#FA6977"
NORTH_COLOR = "#87D5F8"


FILE_RE = re.compile(
    r"bias=(?P<bias>[-+0-9.eE]+)_gamma_mc=(?P<gamma_mc>[-+0-9.eE]+)_gamma_mr=(?P<gamma_mr>[-+0-9.eE]+)_p_scatter=(?P<p>[-+0-9.eE]+)\.h5$"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--sweep-dir",
        type=Path,
        default=None,
        help="Path to a full-sweep directory containing data/*.h5 and intersection.inp",
    )
    parser.add_argument(
        "--num-curves",
        type=int,
        default=8,
        help="Number of representative gamma_mr curves to plot when --gamma-mr is not provided",
    )
    parser.add_argument(
        "--gamma-mr",
        type=float,
        nargs="+",
        default=None,
        help="Explicit gamma_mr values to plot (snapped to nearest available)",
    )
    parser.add_argument(
        "--sample-depth",
        type=int,
        default=2,
        help="Number of grid points sampled inward from contact boundaries",
    )
    parser.add_argument(
        "--gamma-mr-max",
        type=float,
        default=10.0,
        help="Maximum gamma_mr included in plotted curves",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output PNG path (default: <sweep-dir>/ie_over_in_vs_gamma_mc.png)",
    )
    parser.add_argument(
        "--csv",
        type=Path,
        default=None,
        help="Output CSV path (default: <sweep-dir>/ie_over_in_vs_gamma_mc.csv)",
    )
    return parser.parse_args()


def find_latest_intersection_sweep(repo_root: Path) -> Path:
    root = repo_root / "projects" / "simple_geometries" / "full_sweeps"
    candidates = sorted(root.glob("simple_geometries_intersection_sweep_*"))
    if not candidates:
        raise FileNotFoundError(f"No intersection sweep directories found under {root}")
    return max(candidates, key=lambda p: p.stat().st_mtime)


def parse_nodes_and_nsets(inp_path: Path) -> tuple[dict[int, tuple[float, float]], dict[str, set[int]]]:
    nodes: dict[int, tuple[float, float]] = {}
    nsets: dict[str, set[int]] = {}
    section = ""
    active_nset = ""

    for raw in inp_path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith("*"):
            upper = line.upper()
            if upper.startswith("*NODE"):
                section = "node"
                active_nset = ""
            elif upper.startswith("*NSET"):
                match = re.search(r"NSET\s*=\s*([^,\s]+)", line, re.IGNORECASE)
                active_nset = match.group(1) if match else ""
                if active_nset:
                    nsets.setdefault(active_nset, set())
                    section = "nset"
                else:
                    section = ""
            else:
                section = ""
                active_nset = ""
            continue

        if section == "node":
            parts = [p.strip() for p in line.split(",") if p.strip()]
            if len(parts) >= 3:
                nodes[int(parts[0])] = (float(parts[1]), float(parts[2]))
        elif section == "nset" and active_nset:
            for tok in line.split(","):
                tok = tok.strip()
                if tok:
                    nsets[active_nset].add(int(tok))

    return nodes, nsets


def integrate_over_valid_segments(coord: np.ndarray, values: np.ndarray, valid: np.ndarray) -> float:
    idx = np.flatnonzero(valid)
    if idx.size < 2:
        return np.nan
    return float(np.trapezoid(values[idx], coord[idx]))


def format_gamma(v: float) -> str:
    return f"{v:.3g}"


def main() -> None:
    args = parse_args()
    repo_root = Path(__file__).resolve().parents[3]
    sweep_dir = args.sweep_dir if args.sweep_dir is not None else find_latest_intersection_sweep(repo_root)
    sweep_dir = sweep_dir.resolve()
    data_dir = sweep_dir / "data"
    inp_path = sweep_dir / "intersection.inp"

    if not data_dir.is_dir():
        raise FileNotFoundError(f"Missing data directory: {data_dir}")
    if not inp_path.is_file():
        raise FileNotFoundError(f"Missing archived mesh file: {inp_path}")

    files = sorted(data_dir.glob("*.h5"))
    if not files:
        raise FileNotFoundError(f"No .h5 files found in {data_dir}")

    # Parse sweep parameters from filenames.
    entries: list[tuple[Path, float, float]] = []
    for path in files:
        m = FILE_RE.match(path.name)
        if not m:
            continue
        gamma_mc = float(m.group("gamma_mc"))
        gamma_mr = float(m.group("gamma_mr"))
        entries.append((path, gamma_mr, gamma_mc))
    if not entries:
        raise RuntimeError(f"No files matching expected naming pattern in {data_dir}")

    # Read mesh contact geometry for east (contact_side) and north (contact_top).
    nodes, nsets = parse_nodes_and_nsets(inp_path)
    for required in ("contact_side", "contact_top"):
        if required not in nsets:
            raise RuntimeError(f"Missing NSET '{required}' in {inp_path}")

    side_pts = np.array([nodes[i] for i in sorted(nsets["contact_side"])], dtype=float)
    top_pts = np.array([nodes[i] for i in sorted(nsets["contact_top"])], dtype=float)

    east_x = float(np.mean(side_pts[:, 0]))
    east_y_min = float(np.min(side_pts[:, 1]))
    east_y_max = float(np.max(side_pts[:, 1]))

    north_y = float(np.mean(top_pts[:, 1]))
    north_x_min = float(np.min(top_pts[:, 0]))
    north_x_max = float(np.max(top_pts[:, 0]))

    # Grid from one representative file.
    with h5py.File(entries[0][0], "r") as h5:
        x = np.asarray(h5["x"])
        y = np.asarray(h5["y"])

    dx = float(np.mean(np.diff(x)))
    dy = float(np.mean(np.diff(y)))

    ix_east_boundary = int(np.argmin(np.abs(x - east_x)))
    iy_north_boundary = int(np.argmin(np.abs(y - north_y)))
    ix_east = max(0, ix_east_boundary - max(args.sample_depth, 0))
    iy_north = max(0, iy_north_boundary - max(args.sample_depth, 0))

    east_segment = (y >= east_y_min - 0.5 * dy) & (y <= east_y_max + 0.5 * dy)
    north_segment = (x >= north_x_min - 0.5 * dx) & (x <= north_x_max + 0.5 * dx)

    rows: list[tuple[float, float, float, float, float]] = []
    for path, gamma_mr, gamma_mc in entries:
        with h5py.File(path, "r") as h5:
            # Saved arrays are indexed as (y, x), i.e., first axis follows y-grid.
            a1_east = np.asarray(h5["a1"][:, ix_east])
            mask_east = np.asarray(h5["mask"][:, ix_east]).astype(bool)
            b1_north = np.asarray(h5["b1"][iy_north, :])
            mask_north = np.asarray(h5["mask"][iy_north, :]).astype(bool)

        valid_east = east_segment & mask_east & np.isfinite(a1_east)
        valid_north = north_segment & mask_north & np.isfinite(b1_north)

        i_east = integrate_over_valid_segments(y, a1_east, valid_east)
        i_north = integrate_over_valid_segments(x, b1_north, valid_north)
        ratio = np.nan
        if np.isfinite(i_east) and np.isfinite(i_north) and i_north != 0:
            ratio = i_east / i_north

        rows.append((gamma_mr, gamma_mc, i_east, i_north, ratio))

    gamma_mr_vals_all = np.array(sorted({r[0] for r in rows}))
    gamma_mc_vals = np.array(sorted({r[1] for r in rows}))
    gamma_mr_vals = gamma_mr_vals_all[gamma_mr_vals_all <= args.gamma_mr_max]
    if gamma_mr_vals.size == 0:
        raise RuntimeError(f"No gamma_mr values <= {args.gamma_mr_max} available in sweep")

    if args.gamma_mr:
        target = np.array(args.gamma_mr, dtype=float)
        picks = []
        for t in target:
            picks.append(gamma_mr_vals[int(np.argmin(np.abs(np.log10(gamma_mr_vals) - np.log10(t))))])
        selected_gamma_mr = np.array(sorted(set(picks)))
    else:
        n = max(1, min(args.num_curves, len(gamma_mr_vals)))
        idx = np.linspace(0, len(gamma_mr_vals) - 1, n).round().astype(int)
        selected_gamma_mr = gamma_mr_vals[np.unique(idx)]

    out_png = args.output if args.output else sweep_dir / "ie_over_in_vs_gamma_mc.png"
    out_csv = args.csv if args.csv else sweep_dir / "ie_over_in_vs_gamma_mc.csv"

    # Write CSV for all points.
    with out_csv.open("w", encoding="utf-8") as f:
        f.write("gamma_mr,gamma_mc,i_east,i_north,ie_over_in\n")
        for gamma_mr, gamma_mc, i_east, i_north, ratio in sorted(rows, key=lambda r: (r[0], r[1])):
            f.write(f"{gamma_mr:.16g},{gamma_mc:.16g},{i_east:.16g},{i_north:.16g},{ratio:.16g}\n")

    # Plot selected gamma_mr cuts.
    fig, ax = plt.subplots(figsize=(5.7, 3.6), constrained_layout=True)
    fig.patch.set_facecolor("black")
    ax.set_facecolor("black")

    flare_base = sns.color_palette("flare", as_cmap=True)
    # Trim opposite end of flare (per user preference) to avoid the undesired extreme.
    cmap = LinearSegmentedColormap.from_list(
        "flare_trimmed",
        flare_base(np.linspace(0.0, 0.78, 256)),
    )
    norm = colors.LogNorm(vmin=float(np.min(selected_gamma_mr)), vmax=float(np.max(selected_gamma_mr)))

    for gm in selected_gamma_mr[::-1]:
        subset = [(gmc, ratio) for gmr, gmc, _, _, ratio in rows if np.isclose(gmr, gm) and np.isfinite(ratio)]
        subset.sort(key=lambda t: t[0])
        if not subset:
            continue
        gmc_vals = np.array([p[0] for p in subset])
        ratio_vals = np.array([p[1] for p in subset])
        ax.plot(
            gmc_vals,
            ratio_vals,
            "-",
            lw=2.3,
            color=cmap(norm(float(gm))),
        )

    ax.set_xscale("log")
    ax.set_xlabel(r"$\gamma_{mc}$")
    ax.set_ylabel("")
    ax.grid(alpha=0.22, which="both", color="white")
    ax.tick_params(colors="white")
    ax.xaxis.label.set_color("white")
    for spine in ax.spines.values():
        spine.set_color("white")

    # Colored y-label matching intersection contact colors.
    xlab = -0.16
    ax.text(
        xlab,
        0.58,
        r"$I_E$",
        transform=ax.transAxes,
        rotation=90,
        va="center",
        ha="center",
        color=EAST_COLOR,
        fontsize=20,
    )
    ax.text(
        xlab,
        0.50,
        r"$/$",
        transform=ax.transAxes,
        rotation=90,
        va="center",
        ha="center",
        color="white",
        fontsize=18,
    )
    ax.text(
        xlab,
        0.42,
        r"$I_N$",
        transform=ax.transAxes,
        rotation=90,
        va="center",
        ha="center",
        color=NORTH_COLOR,
        fontsize=20,
    )

    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, pad=0.02)
    cbar.set_label(r"$\gamma_{mr}$", color="white")
    cbar.ax.yaxis.set_tick_params(color="white")
    plt.setp(cbar.ax.get_yticklabels(), color="white")
    cbar.outline.set_edgecolor("white")

    fig.savefig(out_png, dpi=700)

    print(f"sweep_dir: {sweep_dir}")
    print(f"data_files: {len(entries)}")
    print(
        f"sampling: x_east={x[ix_east]:.6g} (boundary {x[ix_east_boundary]:.6g}), "
        f"y_north={y[iy_north]:.6g} (boundary {y[iy_north_boundary]:.6g}), depth={args.sample_depth}"
    )
    print(f"gamma_mr_count_total: {len(gamma_mr_vals_all)}, gamma_mr_count_plotted: {len(gamma_mr_vals)}, gamma_mc_count: {len(gamma_mc_vals)}")
    print(f"gamma_mr_max_plot: {args.gamma_mr_max}")
    print("selected_gamma_mr:", ", ".join(format_gamma(float(v)) for v in selected_gamma_mr))
    print(f"saved_csv: {out_csv}")
    print(f"saved_plot: {out_png}")


if __name__ == "__main__":
    main()
