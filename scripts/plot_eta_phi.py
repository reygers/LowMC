#!/usr/bin/env python3
"""Plot an eta-phi heatmap from the two-particle LowMC diagnostic CSV."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create a PNG heatmap from a LowMC eta,phi,S diagnostic CSV."
    )
    parser.add_argument("input_csv", type=Path, help="CSV file with eta,phi,S columns")
    parser.add_argument("output_png", type=Path, help="PNG file to write")
    parser.add_argument(
        "--log",
        action="store_true",
        help="plot log10(S) instead of S",
    )
    return parser.parse_args()


def read_grid(path: Path) -> tuple[list[float], list[float], list[list[float]]]:
    rows: list[tuple[float, float, float]] = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"eta", "phi", "S"}
        if set(reader.fieldnames or []) != required:
            raise ValueError("input CSV must have exactly the columns eta,phi,S")

        for row in reader:
            rows.append((float(row["eta"]), float(row["phi"]), float(row["S"])))

    if not rows:
        raise ValueError("input CSV is empty")

    etas = sorted({row[0] for row in rows})
    phis = sorted({row[1] for row in rows})
    eta_index = {eta: index for index, eta in enumerate(etas)}
    phi_index = {phi: index for index, phi in enumerate(phis)}
    values = [[math.nan for _ in phis] for _ in etas]

    for eta, phi, intensity in rows:
        values[eta_index[eta]][phi_index[phi]] = intensity

    return etas, phis, values


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
        etas, phis, values = read_grid(args.input_csv)
    except (OSError, ValueError) as error:
        print(f"Could not read diagnostic CSV: {error}", file=sys.stderr)
        return 1

    plotted_values = values
    label = "S"
    if args.log:
        plotted_values = [
            [math.log10(value) if value > 0 else math.nan for value in row]
            for row in values
        ]
        label = "log10(S)"

    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
    image = ax.imshow(
        plotted_values,
        origin="lower",
        aspect="auto",
        extent=[min(phis), max(phis), min(etas), max(etas)],
        interpolation="nearest",
    )
    ax.set_xlabel("phi")
    ax.set_ylabel("eta")
    ax.set_title("Two-particle Low-current diagnostic")
    fig.colorbar(image, ax=ax, label=label)
    fig.savefig(args.output_png, dpi=160)
    plt.close(fig)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
