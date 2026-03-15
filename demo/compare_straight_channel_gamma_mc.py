from pathlib import Path
import csv
import math

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover
    import tomli as tomllib

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


COLORS = ["navy", "firebrick", "darkgreen", "goldenrod", "purple"]
MARKERS = ["o", "s", "^", "D", "v"]


def load_run(directory: Path):
    metadata_path = directory / "probe_sweep_metadata.toml"
    csv_path = directory / "probe_sweep.csv"

    with metadata_path.open("rb") as f:
        metadata = tomllib.load(f)
    with csv_path.open() as f:
        rows = list(csv.DictReader(f))

    bias = np.array([float(row["bias"]) for row in rows])
    current = np.array([float(row["integrated_jx"]) for row in rows])
    return metadata, bias, current


def validate_runs(run_metadata):
    shared_float_keys = ("channel_length", "mu0", "mass", "vF", "gamma_mr")
    shared_string_keys = ("convention_name", "transport")
    reference_metadata = run_metadata[0]
    for metadata in run_metadata[1:]:
        for key in shared_string_keys:
            if reference_metadata[key] != metadata[key]:
                raise RuntimeError(f"Metadata mismatch for '{key}'.")
        for key in shared_float_keys:
            if not math.isclose(
                float(reference_metadata[key]),
                float(metadata[key]),
                rel_tol=1e-9,
                abs_tol=1e-9,
            ):
                raise RuntimeError(f"Metadata mismatch for '{key}'.")


def fit_linear_through_origin_first_two(drive, current):
    return float(np.dot(drive[:2], current[:2]) / np.dot(drive[:2], drive[:2]))


def cubic_model(x, a1, a3):
    return a1 * x + a3 * x**3


def fit_residual(bias, current):
    slope = fit_linear_through_origin_first_two(bias, current)
    residual = current - slope * bias
    params, _ = curve_fit(cubic_model, bias, residual)
    fit_x = np.linspace(bias.min(), bias.max(), 300)
    fit_y = cubic_model(fit_x, *params)
    return {
        "slope": slope,
        "residual": residual,
        "a1": params[0],
        "a3": params[1],
        "fit_x": fit_x,
        "fit_y": fit_y,
    }


def main():
    project_root = Path(__file__).resolve().parents[1]
    run_dirs = [
        project_root / "demo" / "data_straight_channel_linearity",
        project_root / "demo" / "data_straight_channel_linearity_gamma_mc_1p0",
        project_root / "demo" / "data_straight_channel_linearity_gamma_mc_100p0",
    ]
    output_dir = project_root / "demo" / "data_straight_channel_linearity_gamma_mc_comparison"
    output_dir.mkdir(parents=True, exist_ok=True)

    loaded = [load_run(directory) for directory in run_dirs]
    run_metadata = [item[0] for item in loaded]
    validate_runs(run_metadata)

    reference_bias = loaded[0][1]
    if not all(np.allclose(reference_bias, bias) for _, bias, _ in loaded[1:]):
        raise RuntimeError("Bias grids do not match across the runs.")

    fit_results = []
    for metadata, bias, current in loaded:
        fit_results.append((float(metadata.get("sweep_gamma_mc", 0.0)), fit_residual(bias, current)))

    overlay_fig, overlay_ax = plt.subplots(figsize=(7.8, 5.2))
    overlay_ax.axhline(0.0, color="black", linewidth=1.2, linestyle="--")
    for (gamma_mc, fit), color, marker in zip(fit_results, COLORS, MARKERS):
        overlay_ax.plot(
            reference_bias,
            fit["residual"],
            marker,
            markersize=6,
            color=color,
            linestyle="none",
            label=rf"data: $\gamma_{{mc}}={gamma_mc:g}$",
        )
        overlay_ax.plot(
            fit["fit_x"],
            fit["fit_y"],
            "-",
            linewidth=2,
            color=color,
            alpha=0.9,
            label=rf"fit: $\gamma_{{mc}}={gamma_mc:g}$",
        )
    overlay_ax.set_xlabel(r"applied bias ($\Delta\mu$)")
    overlay_ax.set_ylabel(r"$I - G_{0,2} x$")
    overlay_ax.set_title("First-two-point residual cubic-fit comparison")
    overlay_ax.grid(True, alpha=0.25)
    overlay_ax.legend(loc="best", ncol=2)
    overlay_fig.tight_layout()
    overlay_path = output_dir / "straight_channel_first_two_residual_cubic_comparison.png"
    overlay_fig.savefig(overlay_path, dpi=200)
    plt.close(overlay_fig)

    reference_gamma, reference_fit = fit_results[0]
    diff_fig, diff_ax = plt.subplots(figsize=(7.8, 5.2))
    diff_ax.axhline(0.0, color="black", linewidth=1.2, linestyle="--")
    for (gamma_mc, fit), color, marker in zip(fit_results[1:], COLORS[1:], MARKERS[1:]):
        residual_difference = fit["residual"] - reference_fit["residual"]
        diff_ax.plot(
            reference_bias,
            residual_difference,
            marker + "-",
            linewidth=2,
            markersize=6,
            color=color,
            label=rf"$\gamma_{{mc}}={gamma_mc:g}$ minus $\gamma_{{mc}}={reference_gamma:g}$",
        )
    diff_ax.set_xlabel(r"applied bias ($\Delta\mu$)")
    diff_ax.set_ylabel(r"$\Delta[(I - G_{0,2}x)]$")
    diff_ax.set_title("Residual difference relative to the reference case")
    diff_ax.grid(True, alpha=0.25)
    diff_ax.legend(loc="best")
    diff_fig.tight_layout()
    diff_path = output_dir / "straight_channel_first_two_residual_difference_vs_reference.png"
    diff_fig.savefig(diff_path, dpi=200)
    plt.close(diff_fig)

    summary_path = output_dir / "straight_channel_first_two_residual_gamma_mc_comparison.txt"
    with summary_path.open("w") as f:
        for gamma_mc, fit in fit_results:
            label = f"gamma_mc_{gamma_mc:g}".replace(".", "p")
            f.write(f"{label}_G0_2={fit['slope']}\n")
            f.write(f"{label}_a1={fit['a1']}\n")
            f.write(f"{label}_a3={fit['a3']}\n")

    print(f"wrote {overlay_path}")
    print(f"wrote {diff_path}")
    print(f"wrote {summary_path}")


if __name__ == "__main__":
    main()
