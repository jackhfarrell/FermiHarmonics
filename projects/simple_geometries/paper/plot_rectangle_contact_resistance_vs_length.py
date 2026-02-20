#!/usr/bin/env python3
"""
Plot contact resistance metrics vs device length from
`run_rectangle_contact_resistance.jl` output.

Usage:
  python projects/simple_geometries/paper/plot_rectangle_contact_resistance_vs_length.py \
    --summary projects/simple_geometries/data/contact_resistance_rectangles/run_YYYY-mm-dd_HHMMSS/contact_resistance_summary.csv
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary", type=Path, required=True, help="Path to contact_resistance_summary.csv")
    parser.add_argument("--out", type=Path, default=None, help="Optional output PNG path")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    summary = args.summary
    if not summary.exists():
        raise FileNotFoundError(f"Summary CSV not found: {summary}")

    df = pd.read_csv(summary).sort_values("length")
    out = args.out or summary.with_name("contact_resistance_vs_length.png")

    fig, ax = plt.subplots(figsize=(7.0, 4.5))
    ax.plot(df["length"], df["Rc_total"], "o-", label="Rc_total")
    ax.plot(df["length"], df["Rc_bottom"], "s--", label="Rc_bottom")
    ax.plot(df["length"], df["Rc_top"], "d--", label="Rc_top")
    ax.plot(df["length"], df["R_bulk"], "^-", label="R_bulk")
    ax.plot(df["length"], df["R_total"], "x-", label="R_total")

    for _, row in df.iterrows():
        ax.annotate(row["device"], (row["length"], row["Rc_total"]), textcoords="offset points", xytext=(0, 8), ha="center")

    ax.set_xlabel("Device length")
    ax.set_ylabel("Resistance proxy (arb. units)")
    ax.set_title("Rectangle Contact Resistance vs Device Length")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    print(f"Saved plot: {out}")


if __name__ == "__main__":
    main()
