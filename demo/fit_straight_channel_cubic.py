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


EXPECTED_REFERENCE_METADATA = {
    "convention_name": "blg_reference_dimensionless",
    "convention_version": 1,
    "channel_length": 1.0,
    "mu0": 1.0,
    "mass": 2.0,
    "vF": 1.0,
    "gamma_mr": 0.0,
    "gamma_mc": 0.0,
    "transport": "parabolic_nonlinear",
    "analysis_case": "reference",
}


def load_reference_metadata(output_dir):
    metadata_path = output_dir / "probe_sweep_metadata.toml"
    if not metadata_path.exists():
        raise RuntimeError(
            f"Missing {metadata_path.name}; rerun the reference straight-channel sweep before plotting."
        )

    with metadata_path.open("rb") as f:
        metadata = tomllib.load(f)

    for key, expected in EXPECTED_REFERENCE_METADATA.items():
        if key not in metadata:
            raise RuntimeError(f"Metadata file is missing required key '{key}'.")
        actual = metadata[key]
        if isinstance(expected, float):
            if not math.isclose(float(actual), expected, rel_tol=1e-9, abs_tol=1e-9):
                raise RuntimeError(
                    f"Metadata mismatch for '{key}': expected {expected}, got {actual}. "
                    "Rerun the reference straight-channel sweep before plotting."
                )
        elif actual != expected:
            raise RuntimeError(
                f"Metadata mismatch for '{key}': expected {expected!r}, got {actual!r}. "
                "Rerun the reference straight-channel sweep before plotting."
            )
    return metadata_path, metadata


def cubic_current(delta_mu, c1, c3):
    return c1 * delta_mu + c3 * delta_mu**3


def differential_conductance(delta_mu, c1, c3):
    return c1 + 3.0 * c3 * delta_mu**2


def differential_resistance(delta_mu, c1, c3):
    return 1.0 / differential_conductance(delta_mu, c1, c3)


def fit_linear_through_origin(drive, integrated_current):
    slope = float(np.dot(drive, integrated_current) / np.dot(drive, drive))
    fit = slope * drive
    residual = integrated_current - fit
    return slope, fit, residual


def fit_linear_through_origin_masked(drive, integrated_current, mask):
    masked_drive = drive[mask]
    masked_current = integrated_current[mask]
    return fit_linear_through_origin(masked_drive, masked_current)[0]


def fit_linear_through_origin_first_n(drive, integrated_current, count):
    return fit_linear_through_origin(drive[:count], integrated_current[:count])[0]


def fit_and_plot(drive, integrated_current, output_dir, stem, xlabel, title_drive, metadata):
    params, covariance = curve_fit(cubic_current, drive, integrated_current)
    c1, c3 = params

    fit_x = np.linspace(drive.min(), drive.max(), 300)
    fit_y = cubic_current(fit_x, c1, c3)
    data_fit = cubic_current(drive, c1, c3)

    residuals = integrated_current - data_fit
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((integrated_current - np.mean(integrated_current)) ** 2)
    r_squared = 1.0 - ss_res / ss_tot

    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    ax.plot(drive, integrated_current, "o", label="integrated current data", markersize=6)
    ax.plot(
        fit_x,
        fit_y,
        "-",
        linewidth=2,
        color="crimson",
        label=rf"fit: $I=c_1x+c_3x^3$" "\n" rf"$c_1={c1:.6f},\ c_3={c3:.6f}$",
    )
    ax.set_xlabel(xlabel)
    ax.set_ylabel("integrated current $I(x_*)$")
    ax.set_title(f"Reference straight-channel integrated current vs {title_drive}")
    ax.legend(loc="upper left")
    ax.grid(True, alpha=0.25)

    figure_path = output_dir / f"{stem}_cubic_curve_fit.png"
    fig.tight_layout()
    fig.savefig(figure_path, dpi=200)
    plt.close(fig)

    diff_x = np.linspace(drive.min(), drive.max(), 300)
    diff_cond = differential_conductance(diff_x, c1, c3)
    diff_res = differential_resistance(diff_x, c1, c3)

    cond_fig, cond_ax = plt.subplots(figsize=(7.2, 5.0))
    cond_ax.plot(diff_x, diff_cond, linewidth=2, color="darkgreen")
    cond_ax.set_xlabel(xlabel)
    cond_ax.set_ylabel(r"$dI/dx$")
    cond_ax.set_title(f"Straight-channel differential conductance vs {title_drive}")
    cond_ax.grid(True, alpha=0.25)
    cond_fig.tight_layout()
    cond_path = output_dir / f"{stem}_differential_conductance.png"
    cond_fig.savefig(cond_path, dpi=200)
    plt.close(cond_fig)

    res_fig, res_ax = plt.subplots(figsize=(7.2, 5.0))
    res_ax.plot(diff_x, diff_res, linewidth=2, color="darkorange")
    res_ax.set_xlabel(xlabel)
    res_ax.set_ylabel(r"$dx/dI$")
    res_ax.set_title(f"Straight-channel differential resistance vs {title_drive}")
    res_ax.grid(True, alpha=0.25)
    res_fig.tight_layout()
    res_path = output_dir / f"{stem}_differential_resistance.png"
    res_fig.savefig(res_path, dpi=200)
    plt.close(res_fig)

    summary_path = output_dir / f"{stem}_cubic_curve_fit.txt"
    with summary_path.open("w") as f:
        f.write("Model: I(x) = c1 * x + c3 * x^3\n")
        f.write(f"convention_name={metadata['convention_name']}\n")
        f.write(f"channel_length={metadata['channel_length']}\n")
        f.write(f"mu0={metadata['mu0']}\n")
        f.write(f"mass={metadata['mass']}\n")
        f.write(f"gamma_mr={metadata['gamma_mr']}\n")
        f.write(f"gamma_mc={metadata['gamma_mc']}\n")
        f.write(f"c1={c1}\n")
        f.write(f"c3={c3}\n")
        f.write(f"r_squared={r_squared}\n")
        f.write(f"residual_l2={np.sqrt(np.mean(residuals**2))}\n")

    return {
        "figure_path": figure_path,
        "cond_path": cond_path,
        "res_path": res_path,
        "summary_path": summary_path,
        "c1": c1,
        "c3": c3,
        "r_squared": r_squared,
    }


def make_linear_subtracted_plot(drive, integrated_current, output_dir, stem, xlabel, title_drive):
    slope, fit, residual = fit_linear_through_origin(drive, integrated_current)

    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    ax.axhline(0.0, color="black", linewidth=1.2, linestyle="--")
    ax.plot(drive, residual, "o-", linewidth=2, markersize=6, color="teal")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(r"$I - G_{\rm lin} x$")
    ax.set_title(f"Linear-fit-subtracted current vs {title_drive}")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()

    figure_path = output_dir / f"{stem}_linear_subtracted.png"
    fig.savefig(figure_path, dpi=200)
    plt.close(fig)

    summary_path = output_dir / f"{stem}_linear_subtracted.txt"
    with summary_path.open("w") as f:
        f.write("Model: I_linear(x) = G_lin * x\n")
        f.write(f"G_lin={slope}\n")
        f.write(f"max_abs_residual={np.max(np.abs(residual))}\n")

    return {
        "figure_path": figure_path,
        "summary_path": summary_path,
        "slope": slope,
    }


def make_low_drive_subtracted_plot(drive, integrated_current, low_drive_mask, output_dir, stem, xlabel, title_drive):
    slope = fit_linear_through_origin_masked(drive, integrated_current, low_drive_mask)
    residual = integrated_current - slope * drive

    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    ax.axhline(0.0, color="black", linewidth=1.2, linestyle="--")
    ax.plot(drive, residual, "o-", linewidth=2, markersize=6, color="firebrick")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(r"$I - G_0 x$")
    ax.set_title(f"Low-drive-subtracted current vs {title_drive}")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()

    figure_path = output_dir / f"{stem}_low_drive_subtracted.png"
    fig.savefig(figure_path, dpi=200)
    plt.close(fig)

    summary_path = output_dir / f"{stem}_low_drive_subtracted.txt"
    with summary_path.open("w") as f:
        f.write("Model: I_low_drive(x) = G0 * x\n")
        f.write(f"G0={slope}\n")
        f.write(f"low_drive_points={int(np.count_nonzero(low_drive_mask))}\n")
        f.write(f"max_abs_residual={np.max(np.abs(residual))}\n")

    return {
        "figure_path": figure_path,
        "summary_path": summary_path,
        "slope": slope,
    }


def make_cubic_scaled_plot(drive, integrated_current, slope, output_dir, stem, xlabel, title_drive):
    cubic_scaled = (integrated_current - slope * drive) / drive**3

    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    ax.plot(drive, cubic_scaled, "o-", linewidth=2, markersize=6, color="purple")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(r"$(I - G_0 x) / x^3$")
    ax.set_title(f"Cubic-scaled nonlinear signal vs {title_drive}")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()

    figure_path = output_dir / f"{stem}_cubic_scaled.png"
    fig.savefig(figure_path, dpi=200)
    plt.close(fig)

    summary_path = output_dir / f"{stem}_cubic_scaled.txt"
    with summary_path.open("w") as f:
        f.write("Model: nonlinear_signal(x) = (I - G0 * x) / x^3\n")
        f.write(f"G0={slope}\n")
        f.write(f"min_value={np.min(cubic_scaled)}\n")
        f.write(f"max_value={np.max(cubic_scaled)}\n")

    return {
        "figure_path": figure_path,
        "summary_path": summary_path,
    }


def make_first_two_subtracted_plot(drive, integrated_current, output_dir, stem, xlabel, title_drive):
    slope = fit_linear_through_origin_first_n(drive, integrated_current, 2)
    residual = integrated_current - slope * drive

    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    ax.axhline(0.0, color="black", linewidth=1.2, linestyle="--")
    ax.plot(drive, residual, "o-", linewidth=2, markersize=6, color="darkmagenta")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(r"$I - G_{0,2} x$")
    ax.set_title(f"First-two-point-subtracted current vs {title_drive}")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()

    figure_path = output_dir / f"{stem}_first_two_subtracted.png"
    fig.savefig(figure_path, dpi=200)
    plt.close(fig)

    summary_path = output_dir / f"{stem}_first_two_subtracted.txt"
    with summary_path.open("w") as f:
        f.write("Model: I_first_two(x) = G0_2 * x\n")
        f.write(f"G0_2={slope}\n")
        f.write("fit_points=2\n")
        f.write(f"max_abs_residual={np.max(np.abs(residual))}\n")

    return {
        "figure_path": figure_path,
        "summary_path": summary_path,
        "slope": slope,
    }


def make_residual_cubic_fit_plot(drive, integrated_current, linear_slope, output_dir, stem, xlabel, title_drive):
    residual = integrated_current - linear_slope * drive

    def residual_cubic(x, a1, a3):
        return a1 * x + a3 * x**3

    params, _ = curve_fit(residual_cubic, drive, residual)
    a1, a3 = params
    fit_x = np.linspace(drive.min(), drive.max(), 300)
    fit_y = residual_cubic(fit_x, a1, a3)
    data_fit = residual_cubic(drive, a1, a3)

    ss_res = np.sum((residual - data_fit) ** 2)
    ss_tot = np.sum((residual - np.mean(residual)) ** 2)
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")

    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    ax.axhline(0.0, color="black", linewidth=1.2, linestyle="--")
    ax.plot(drive, residual, "o", markersize=6, color="darkmagenta", label="residual data")
    ax.plot(
        fit_x,
        fit_y,
        "-",
        linewidth=2,
        color="black",
        label=rf"fit: $a_1 x + a_3 x^3$" "\n" rf"$a_1={a1:.6g},\ a_3={a3:.6g}$",
    )
    ax.set_xlabel(xlabel)
    ax.set_ylabel(r"$I - G_{0,2} x$")
    ax.set_title(f"First-two-point residual cubic fit vs {title_drive}")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best")
    fig.tight_layout()

    figure_path = output_dir / f"{stem}_first_two_residual_cubic_fit.png"
    fig.savefig(figure_path, dpi=200)
    plt.close(fig)

    summary_path = output_dir / f"{stem}_first_two_residual_cubic_fit.txt"
    with summary_path.open("w") as f:
        f.write("Model: residual(x) = a1 * x + a3 * x^3\n")
        f.write(f"G0_2={linear_slope}\n")
        f.write(f"a1={a1}\n")
        f.write(f"a3={a3}\n")
        f.write(f"r_squared={r_squared}\n")

    return {
        "figure_path": figure_path,
        "summary_path": summary_path,
        "a1": a1,
        "a3": a3,
        "r_squared": r_squared,
    }


def make_comparison_plot(applied_bias, measured_drop, output_dir):
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    ax.plot(applied_bias, measured_drop, "o-", linewidth=2, markersize=6, color="navy")
    ax.plot(applied_bias, applied_bias, "--", linewidth=1.5, color="gray", label="measured drop = applied bias")
    ax.set_xlabel("applied bias ($\\Delta\\mu$)")
    ax.set_ylabel("measured average drop ($\\Delta a_0$)")
    ax.set_title("Measured average drop vs applied bias")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="upper left")
    fig.tight_layout()

    comparison_path = output_dir / "straight_channel_measured_drop_vs_bias.png"
    fig.savefig(comparison_path, dpi=200)
    plt.close(fig)
    return comparison_path


def make_conductance_summary_plot(applied_bias, measured_drop, integrated_current, output_dir):
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    ax.plot(
        applied_bias,
        integrated_current / applied_bias,
        "o-",
        linewidth=2,
        markersize=6,
        color="firebrick",
        label=r"$I / \Delta\mu$",
    )
    ax.plot(
        measured_drop,
        integrated_current / measured_drop,
        "s-",
        linewidth=2,
        markersize=6,
        color="navy",
        label=r"$I / \Delta a_0$",
    )
    ax.set_xlabel("drive")
    ax.set_ylabel("effective conductance")
    ax.set_title("Effective conductance summary")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best")
    fig.tight_layout()

    figure_path = output_dir / "straight_channel_effective_conductance_summary.png"
    fig.savefig(figure_path, dpi=200)
    plt.close(fig)
    return figure_path


def main():
    project_root = Path(__file__).resolve().parents[1]
    data_path = project_root / "demo" / "data_straight_channel_linearity" / "probe_sweep.csv"
    output_dir = project_root / "demo" / "data_straight_channel_linearity"
    metadata_path, metadata = load_reference_metadata(output_dir)

    with data_path.open() as f:
        rows = list(csv.DictReader(f))

    applied_bias = np.array([float(row["bias"]) for row in rows])
    integrated_current = np.array([float(row["integrated_jx"]) for row in rows])
    measured_drop = np.array([float(row["measured_avg_drop"]) for row in rows])
    low_drive_mask = applied_bias <= 0.2154434690031884

    bias_fit = fit_and_plot(
        applied_bias,
        integrated_current,
        output_dir,
        "straight_channel_integrated_current_vs_bias",
        "applied bias ($\\Delta\\mu$)",
        "applied bias",
        metadata,
    )
    measured_fit = fit_and_plot(
        measured_drop,
        integrated_current,
        output_dir,
        "straight_channel_integrated_current_vs_measured_drop",
        "measured average drop ($\\Delta a_0$)",
        "measured average drop",
        metadata,
    )
    bias_linear_subtracted = make_linear_subtracted_plot(
        applied_bias,
        integrated_current,
        output_dir,
        "straight_channel_integrated_current_vs_bias",
        "applied bias ($\\Delta\\mu$)",
        "applied bias",
    )
    measured_linear_subtracted = make_linear_subtracted_plot(
        measured_drop,
        integrated_current,
        output_dir,
        "straight_channel_integrated_current_vs_measured_drop",
        "measured average drop ($\\Delta a_0$)",
        "measured average drop",
    )
    bias_low_drive_subtracted = make_low_drive_subtracted_plot(
        applied_bias,
        integrated_current,
        low_drive_mask,
        output_dir,
        "straight_channel_integrated_current_vs_bias",
        "applied bias ($\\Delta\\mu$)",
        "applied bias",
    )
    measured_low_drive_subtracted = make_low_drive_subtracted_plot(
        measured_drop,
        integrated_current,
        low_drive_mask,
        output_dir,
        "straight_channel_integrated_current_vs_measured_drop",
        "measured average drop ($\\Delta a_0$)",
        "measured average drop",
    )
    bias_cubic_scaled = make_cubic_scaled_plot(
        applied_bias,
        integrated_current,
        bias_low_drive_subtracted["slope"],
        output_dir,
        "straight_channel_integrated_current_vs_bias",
        "applied bias ($\\Delta\\mu$)",
        "applied bias",
    )
    measured_cubic_scaled = make_cubic_scaled_plot(
        measured_drop,
        integrated_current,
        measured_low_drive_subtracted["slope"],
        output_dir,
        "straight_channel_integrated_current_vs_measured_drop",
        "measured average drop ($\\Delta a_0$)",
        "measured average drop",
    )
    bias_first_two_subtracted = make_first_two_subtracted_plot(
        applied_bias,
        integrated_current,
        output_dir,
        "straight_channel_integrated_current_vs_bias",
        "applied bias ($\\Delta\\mu$)",
        "applied bias",
    )
    measured_first_two_subtracted = make_first_two_subtracted_plot(
        measured_drop,
        integrated_current,
        output_dir,
        "straight_channel_integrated_current_vs_measured_drop",
        "measured average drop ($\\Delta a_0$)",
        "measured average drop",
    )
    bias_residual_cubic_fit = make_residual_cubic_fit_plot(
        applied_bias,
        integrated_current,
        bias_first_two_subtracted["slope"],
        output_dir,
        "straight_channel_integrated_current_vs_bias",
        "applied bias ($\\Delta\\mu$)",
        "applied bias",
    )
    measured_residual_cubic_fit = make_residual_cubic_fit_plot(
        measured_drop,
        integrated_current,
        measured_first_two_subtracted["slope"],
        output_dir,
        "straight_channel_integrated_current_vs_measured_drop",
        "measured average drop ($\\Delta a_0$)",
        "measured average drop",
    )
    comparison_path = make_comparison_plot(applied_bias, measured_drop, output_dir)
    conductance_summary_path = make_conductance_summary_plot(
        applied_bias,
        measured_drop,
        integrated_current,
        output_dir,
    )

    print(f"wrote {bias_fit['figure_path']}")
    print(f"wrote {bias_fit['cond_path']}")
    print(f"wrote {bias_fit['res_path']}")
    print(f"wrote {bias_fit['summary_path']}")
    print(f"wrote {measured_fit['figure_path']}")
    print(f"wrote {measured_fit['cond_path']}")
    print(f"wrote {measured_fit['res_path']}")
    print(f"wrote {measured_fit['summary_path']}")
    print(f"validated {metadata_path}")
    print(f"wrote {bias_linear_subtracted['figure_path']}")
    print(f"wrote {bias_linear_subtracted['summary_path']}")
    print(f"wrote {measured_linear_subtracted['figure_path']}")
    print(f"wrote {measured_linear_subtracted['summary_path']}")
    print(f"wrote {bias_low_drive_subtracted['figure_path']}")
    print(f"wrote {bias_low_drive_subtracted['summary_path']}")
    print(f"wrote {measured_low_drive_subtracted['figure_path']}")
    print(f"wrote {measured_low_drive_subtracted['summary_path']}")
    print(f"wrote {bias_cubic_scaled['figure_path']}")
    print(f"wrote {bias_cubic_scaled['summary_path']}")
    print(f"wrote {measured_cubic_scaled['figure_path']}")
    print(f"wrote {measured_cubic_scaled['summary_path']}")
    print(f"wrote {bias_first_two_subtracted['figure_path']}")
    print(f"wrote {bias_first_two_subtracted['summary_path']}")
    print(f"wrote {measured_first_two_subtracted['figure_path']}")
    print(f"wrote {measured_first_two_subtracted['summary_path']}")
    print(f"wrote {bias_residual_cubic_fit['figure_path']}")
    print(f"wrote {bias_residual_cubic_fit['summary_path']}")
    print(f"wrote {measured_residual_cubic_fit['figure_path']}")
    print(f"wrote {measured_residual_cubic_fit['summary_path']}")
    print(f"wrote {comparison_path}")
    print(f"wrote {conductance_summary_path}")
    print(f"bias fit: c1={bias_fit['c1']}, c3={bias_fit['c3']}, R^2={bias_fit['r_squared']}")
    print(
        "measured-drop fit: "
        f"c1={measured_fit['c1']}, c3={measured_fit['c3']}, R^2={measured_fit['r_squared']}"
    )
    print(
        "linear-through-origin slopes: "
        f"bias={bias_linear_subtracted['slope']}, "
        f"measured_drop={measured_linear_subtracted['slope']}"
    )
    print(
        "low-drive slopes: "
        f"bias={bias_low_drive_subtracted['slope']}, "
        f"measured_drop={measured_low_drive_subtracted['slope']}"
    )
    print(
        "first-two-point slopes: "
        f"bias={bias_first_two_subtracted['slope']}, "
        f"measured_drop={measured_first_two_subtracted['slope']}"
    )
    print(
        "first-two residual cubic fits: "
        f"bias(a1={bias_residual_cubic_fit['a1']}, a3={bias_residual_cubic_fit['a3']}, R^2={bias_residual_cubic_fit['r_squared']}), "
        f"measured_drop(a1={measured_residual_cubic_fit['a1']}, a3={measured_residual_cubic_fit['a3']}, R^2={measured_residual_cubic_fit['r_squared']})"
    )


if __name__ == "__main__":
    main()
