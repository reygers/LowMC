#!/usr/bin/env python3
"""Plot weighted-sampler study summary CSV output."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path


REQUIRED_COLUMNS = {
    "momentum",
    "meanWeight",
    "varianceWeight",
    "effectiveSampleFraction",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create a PNG summary plot from a LowMC weighted-sampler study CSV."
    )
    parser.add_argument(
        "input_csv",
        type=Path,
        help="CSV file written by lowmc_weighted_sampler_study",
    )
    parser.add_argument("output_png", type=Path, help="PNG file to write")
    parser.add_argument(
        "--linear-x",
        action="store_true",
        help="use a linear momentum axis instead of the default log scale",
    )
    return parser.parse_args()


def read_rows(path: Path) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []

    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = set(reader.fieldnames or [])
        missing = REQUIRED_COLUMNS - fieldnames
        if missing:
            names = ", ".join(sorted(missing))
            raise ValueError(f"input CSV is missing required column(s): {names}")

        for row in reader:
            momentum = float(row["momentum"])
            mean_weight = float(row["meanWeight"])
            variance_weight = float(row["varianceWeight"])
            effective_fraction = float(row["effectiveSampleFraction"])
            if momentum <= 0.0:
                raise ValueError("momentum values must be positive")
            if not all(
                math.isfinite(value)
                for value in (momentum, mean_weight, variance_weight, effective_fraction)
            ):
                raise ValueError("input CSV contains non-finite values")
            rows.append(
                {
                    "momentum": momentum,
                    "meanWeight": mean_weight,
                    "varianceWeight": variance_weight,
                    "effectiveSampleFraction": effective_fraction,
                }
            )

    if not rows:
        raise ValueError("input CSV is empty")

    return sorted(rows, key=lambda row: row["momentum"])


def main() -> int:
    args = parse_args()

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print(
            "matplotlib is required for plotting. Install it or run only the C++ diagnostic.",
            file=sys.stderr,
        )
        return 2

    try:
        rows = read_rows(args.input_csv)
    except (OSError, ValueError) as error:
        print(f"Could not read sampler-study CSV: {error}", file=sys.stderr)
        return 1

    momenta = [row["momentum"] for row in rows]
    mean_weights = [row["meanWeight"] for row in rows]
    rms_weights = [math.sqrt(max(row["varianceWeight"], 0.0)) for row in rows]
    effective_fractions = [row["effectiveSampleFraction"] for row in rows]

    fig, axes = plt.subplots(3, 1, figsize=(7, 8), sharex=True, constrained_layout=True)

    axes[0].plot(momenta, mean_weights, marker="o")
    axes[0].set_ylabel("mean weight")
    axes[0].set_title("Weighted sampler study")

    axes[1].plot(momenta, rms_weights, marker="o", color="tab:orange")
    axes[1].set_ylabel("weight RMS")

    axes[2].plot(momenta, effective_fractions, marker="o", color="tab:green")
    axes[2].set_ylabel("ESS / samples")
    axes[2].set_xlabel("charged-pion momentum [GeV]")
    axes[2].set_ylim(bottom=0.0)

    if not args.linear_x:
        for ax in axes:
            ax.set_xscale("log")

    for ax in axes:
        ax.grid(True, which="both", alpha=0.25)

    fig.savefig(args.output_png, dpi=160)
    plt.close(fig)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
