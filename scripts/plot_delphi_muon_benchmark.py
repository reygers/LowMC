#!/usr/bin/env python3
"""Overlay LowMC muon-FSR angle samples with digitized DELPHI Fig. 1a."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot a DELPHI Fig. 1a shape comparison from LowMC theta samples."
    )
    parser.add_argument(
        "samples_csv",
        type=Path,
        help="CSV from lowmc_delphi_muon_fsr_diagnostic",
    )
    parser.add_argument("output_png", type=Path, help="PNG file to write")
    parser.add_argument(
        "--reference-csv",
        type=Path,
        default=Path("fixtures/reference/delphi_fig1a_muon_fsr.csv"),
        help="digitized DELPHI Fig. 1a CSV",
    )
    parser.add_argument("--bins", type=int, default=100, help="histogram bins")
    parser.add_argument(
        "--theta-max-mrad",
        type=float,
        default=10.0,
        help="maximum theta shown in mrad",
    )
    parser.add_argument(
        "--scale",
        choices=("peak", "area", "none"),
        default="area",
        help="normalization applied to the LowMC histogram",
    )
    return parser.parse_args()


def read_reference(path: Path) -> tuple[list[float], list[float]]:
    theta: list[float] = []
    entries: list[float] = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"theta_mrad", "entries_per_0p1_mrad"}
        if not required.issubset(set(reader.fieldnames or [])):
            raise ValueError(
                "reference CSV must include theta_mrad and entries_per_0p1_mrad"
            )
        for row in reader:
            theta.append(float(row["theta_mrad"]))
            entries.append(float(row["entries_per_0p1_mrad"]))
    if not theta:
        raise ValueError("reference CSV is empty")
    return theta, entries


def read_samples(path: Path, theta_max_mrad: float) -> list[float]:
    samples: list[float] = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if "theta_mrad" not in set(reader.fieldnames or []):
            raise ValueError("samples CSV must include theta_mrad")
        for row in reader:
            theta = float(row["theta_mrad"])
            if 0.0 <= theta <= theta_max_mrad:
                samples.append(theta)
    if not samples:
        raise ValueError("samples CSV contains no theta_mrad values in the plot range")
    return samples


def histogram_density(samples: list[float], bins: int, theta_max_mrad: float) -> tuple[list[float], list[float]]:
    if bins <= 0:
        raise ValueError("--bins must be positive")
    width = theta_max_mrad / bins
    counts = [0 for _ in range(bins)]
    for theta in samples:
        index = min(bins - 1, int(theta / width))
        counts[index] += 1
    centers = [(index + 0.5) * width for index in range(bins)]
    density = [count / (len(samples) * width) for count in counts]
    return centers, density


def trapezoid_area(x_values: list[float], y_values: list[float]) -> float:
    return sum(
        0.5 * (right_y + left_y) * (right_x - left_x)
        for left_x, right_x, left_y, right_y in zip(
            x_values, x_values[1:], y_values, y_values[1:]
        )
    )


def scale_histogram(
    hist_x: list[float],
    hist_y: list[float],
    ref_x: list[float],
    ref_y: list[float],
    mode: str,
) -> list[float]:
    if mode == "none":
        return hist_y
    if mode == "peak":
        hist_peak = max(hist_y)
        ref_peak = max(ref_y)
        if hist_peak <= 0.0:
            raise ValueError("LowMC histogram peak is zero")
        return [value * ref_peak / hist_peak for value in hist_y]
    if mode == "area":
        hist_area = trapezoid_area(hist_x, hist_y)
        ref_area = trapezoid_area(ref_x, ref_y)
        if hist_area <= 0.0:
            raise ValueError("LowMC histogram area is zero")
        return [value * ref_area / hist_area for value in hist_y]
    raise ValueError(f"unknown scale mode: {mode}")


def main() -> int:
    args = parse_args()
    if args.theta_max_mrad <= 0.0:
        print("--theta-max-mrad must be positive", file=sys.stderr)
        return 2

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is required for plotting.", file=sys.stderr)
        return 2

    try:
        ref_x, ref_y = read_reference(args.reference_csv)
        samples = read_samples(args.samples_csv, args.theta_max_mrad)
        hist_x, hist_y_raw = histogram_density(samples, args.bins, args.theta_max_mrad)
        hist_y = scale_histogram(hist_x, hist_y_raw, ref_x, ref_y, args.scale)
    except (OSError, ValueError) as error:
        print(f"Could not prepare DELPHI benchmark plot: {error}", file=sys.stderr)
        return 1

    fig, ax = plt.subplots(figsize=(7.0, 4.6), constrained_layout=True)
    ax.plot(ref_x, ref_y, color="#2563eb", linewidth=2.0, label="DELPHI Fig. 1a digitized")
    ax.step(
        hist_x,
        hist_y,
        where="mid",
        color="#16a34a",
        linewidth=1.6,
        label=f"LowMC samples ({len(samples)} in range)",
    )
    ax.axvline(
        math.sqrt(3.0) / 430.0 * 1000.0,
        color="#525252",
        linestyle=":",
        linewidth=1.2,
        label=r"$\sqrt{3}/\Gamma$, $\Gamma=430$",
    )
    ax.set_xlim(0.0, args.theta_max_mrad)
    ax.set_ylim(bottom=0.0)
    ax.set_xlabel(r"$\theta_\gamma$ relative to $\mu^+$ [mrad]")
    ax.set_ylabel("entries / 0.1 mrad, arbitrary normalization")
    ax.set_title(
        "Muon final-state radiation angular shape\n"
        "Reference: DELPHI, Eur. Phys. J. C 57, 499-514 (2008)"
    )
    ax.legend(frameon=False)
    fig.savefig(args.output_png, dpi=180)
    plt.close(fig)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
