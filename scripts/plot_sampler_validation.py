#!/usr/bin/env python3
"""Plot fixed-kT truth, proposal occupancy, and weighted-sample validation views."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create a fixed-kT sampler validation plot with 2D heatmaps and "
            "1D eta/phi projections."
        )
    )
    parser.add_argument("truth_csv", type=Path, help="CSV with eta,phi,S columns")
    parser.add_argument(
        "sample_csv",
        type=Path,
        help="CSV from lowmc_fixed_kt_weighted_samples with eta,phi,weight columns",
    )
    parser.add_argument("output_png", type=Path, help="PNG plot to write")
    parser.add_argument("--title", default="", help="optional plot title")
    parser.add_argument(
        "--truth-label",
        default="truth grid",
        help="label for the deterministic reference shown in blue",
    )
    parser.add_argument(
        "--reference-label",
        default="",
        help="optional label for an external integral reference",
    )
    parser.add_argument(
        "--reference-integral",
        type=float,
        default=None,
        help="optional external reference integral to print in the title",
    )
    parser.add_argument(
        "--reference-uncertainty",
        type=float,
        default=None,
        help="optional one-sigma uncertainty for --reference-integral",
    )
    parser.add_argument(
        "--warning",
        default="",
        help="optional warning text printed below the title",
    )
    parser.add_argument(
        "--zoom-half-width",
        type=float,
        default=0.35,
        help="half-width of the automatic projection zoom windows; use 0 to disable zoom panels",
    )
    parser.add_argument(
        "--alpha",
        type=float,
        default=1.0 / 137.035999084,
        help="fine-structure constant used to convert S to Low density",
    )
    return parser.parse_args()


def read_truth(path: Path) -> tuple[list[float], list[float], list[list[float]]]:
    rows: list[tuple[float, float, float]] = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"eta", "phi", "S"}
        if set(reader.fieldnames or []) != required:
            raise ValueError("truth CSV must have exactly the columns eta,phi,S")
        for row in reader:
            rows.append((float(row["eta"]), float(row["phi"]), float(row["S"])))

    if not rows:
        raise ValueError("truth CSV is empty")

    etas = sorted({row[0] for row in rows})
    phis = sorted({row[1] for row in rows})
    eta_index = {eta: index for index, eta in enumerate(etas)}
    phi_index = {phi: index for index, phi in enumerate(phis)}
    values = [[math.nan for _ in phis] for _ in etas]
    for eta, phi, intensity in rows:
        values[eta_index[eta]][phi_index[phi]] = intensity

    return etas, phis, values


def read_samples(path: Path) -> list[tuple[float, float, float]]:
    rows: list[tuple[float, float, float]] = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = set(reader.fieldnames or [])
        required = {"eta", "phi", "weight"}
        if not required.issubset(fieldnames):
            raise ValueError("sample CSV must contain eta,phi,weight columns")
        for row in reader:
            rows.append((float(row["eta"]), float(row["phi"]), float(row["weight"])))
    if not rows:
        raise ValueError("sample CSV is empty")
    return rows


def bin_edges(centers: list[float]) -> list[float]:
    if len(centers) < 2:
        raise ValueError("need at least two grid centers")
    edges = [centers[0] - 0.5 * (centers[1] - centers[0])]
    for index in range(len(centers) - 1):
        edges.append(0.5 * (centers[index] + centers[index + 1]))
    edges.append(centers[-1] + 0.5 * (centers[-1] - centers[-2]))
    return edges


def normalize_shape(values):
    import numpy as np

    total = float(np.nansum(values))
    if total <= 0.0:
        return values
    return values / total


def positive_log(values):
    import numpy as np

    positive = np.asarray(values, dtype=float)
    positive[positive <= 0.0] = np.nan
    return np.log10(positive)


def set_zoom_limits(axis, centers, peak_index: int, half_width: float) -> None:
    if half_width <= 0.0:
        return
    lower = max(float(centers[0]), float(centers[peak_index]) - half_width)
    upper = min(float(centers[-1]), float(centers[peak_index]) + half_width)
    if upper > lower:
        axis.set_xlim(lower, upper)


def main() -> int:
    args = parse_args()
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        print("matplotlib and numpy are required for sampler validation plots", file=sys.stderr)
        return 2

    try:
        etas, phis, truth_values = read_truth(args.truth_csv)
        samples = read_samples(args.sample_csv)
    except (OSError, ValueError) as error:
        print(f"Could not read validation inputs: {error}", file=sys.stderr)
        return 1

    eta_edges = np.asarray(bin_edges(etas))
    phi_edges = np.asarray(bin_edges(phis))
    truth = np.asarray(truth_values, dtype=float)
    sample_eta = np.asarray([row[0] for row in samples])
    sample_phi = np.asarray([row[1] for row in samples])
    sample_weight = np.asarray([row[2] for row in samples])

    occupancy_phi_eta, _, _ = np.histogram2d(
        sample_phi,
        sample_eta,
        bins=[phi_edges, eta_edges],
    )
    weighted_phi_eta, _, _ = np.histogram2d(
        sample_phi,
        sample_eta,
        bins=[phi_edges, eta_edges],
        weights=sample_weight,
    )
    occupancy = occupancy_phi_eta.T
    weighted_sum = weighted_phi_eta.T
    weighted_estimate = weighted_sum / float(len(samples))

    eta_width = float(np.mean(np.diff(eta_edges)))
    phi_width = float(np.mean(np.diff(phi_edges)))
    bin_area = eta_width * phi_width
    alpha_factor = args.alpha / ((2.0 * math.pi) * (2.0 * math.pi))
    truth_density = alpha_factor * truth
    truth_integral = float(np.nansum(truth_density) * bin_area)
    weighted_integral = float(np.nansum(weighted_estimate))

    truth_shape = normalize_shape(truth_density)
    occupancy_shape = normalize_shape(occupancy)
    weighted_shape = normalize_shape(weighted_estimate)

    truth_eta = normalize_shape(np.nansum(truth_density, axis=1))
    occupancy_eta = normalize_shape(np.nansum(occupancy, axis=1))
    weighted_eta = normalize_shape(np.nansum(weighted_estimate, axis=1))
    truth_phi = normalize_shape(np.nansum(truth_density, axis=0))
    occupancy_phi = normalize_shape(np.nansum(occupancy, axis=0))
    weighted_phi = normalize_shape(np.nansum(weighted_estimate, axis=0))
    eta_peak_index = int(np.nanargmax(truth_eta))
    phi_peak_index = int(np.nanargmax(truth_phi))

    if args.zoom_half_width > 0.0:
        fig = plt.figure(figsize=(14, 10), constrained_layout=True)
        axes = fig.subplot_mosaic(
            [
                ["truth", "occupancy", "weighted"],
                ["eta", "eta", "phi"],
                ["eta_zoom", "eta_zoom", "phi_zoom"],
            ],
            height_ratios=[2.0, 1.0, 1.0],
        )
    else:
        fig = plt.figure(figsize=(14, 8), constrained_layout=True)
        axes = fig.subplot_mosaic(
            [
                ["truth", "occupancy", "weighted"],
                ["eta", "eta", "phi"],
            ],
            height_ratios=[2.0, 1.0],
        )

    extent = [float(phi_edges[0]), float(phi_edges[-1]), float(eta_edges[0]), float(eta_edges[-1])]
    heatmaps = [
        ("truth", f"{args.truth_label} shape log10", truth_shape),
        ("occupancy", "Proposal occupancy log10", occupancy_shape),
        ("weighted", "Weighted sample shape log10", weighted_shape),
    ]
    for key, title, values in heatmaps:
        image = axes[key].imshow(
            positive_log(values),
            origin="lower",
            aspect="auto",
            extent=extent,
            interpolation="nearest",
        )
        axes[key].set_title(title)
        axes[key].set_xlabel("phi")
        axes[key].set_ylabel("eta")
        fig.colorbar(image, ax=axes[key], fraction=0.046, pad=0.04)

    axes["eta"].plot(etas, truth_eta, label=args.truth_label, linewidth=2)
    axes["eta"].plot(etas, occupancy_eta, label="occupancy", linewidth=1.5)
    axes["eta"].plot(etas, weighted_eta, label="weighted", linewidth=1.5)
    axes["eta"].set_xlabel("eta")
    axes["eta"].set_ylabel("normalized projection")
    axes["eta"].set_title("eta projection")
    axes["eta"].legend()

    axes["phi"].plot(phis, truth_phi, label=args.truth_label, linewidth=2)
    axes["phi"].plot(phis, occupancy_phi, label="occupancy", linewidth=1.5)
    axes["phi"].plot(phis, weighted_phi, label="weighted", linewidth=1.5)
    axes["phi"].set_xlabel("phi")
    axes["phi"].set_ylabel("normalized projection")
    axes["phi"].set_title("phi projection")
    axes["phi"].legend()

    if args.zoom_half_width > 0.0:
        axes["eta_zoom"].plot(etas, truth_eta, label=args.truth_label, linewidth=2)
        axes["eta_zoom"].plot(etas, occupancy_eta, label="occupancy", linewidth=1.5)
        axes["eta_zoom"].plot(etas, weighted_eta, label="weighted", linewidth=1.5)
        axes["eta_zoom"].set_xlabel("eta")
        axes["eta_zoom"].set_ylabel("normalized projection")
        axes["eta_zoom"].set_title("eta projection zoom")
        set_zoom_limits(axes["eta_zoom"], etas, eta_peak_index, args.zoom_half_width)
        axes["eta_zoom"].legend()

        axes["phi_zoom"].plot(phis, truth_phi, label=args.truth_label, linewidth=2)
        axes["phi_zoom"].plot(phis, occupancy_phi, label="occupancy", linewidth=1.5)
        axes["phi_zoom"].plot(phis, weighted_phi, label="weighted", linewidth=1.5)
        axes["phi_zoom"].set_xlabel("phi")
        axes["phi_zoom"].set_ylabel("normalized projection")
        axes["phi_zoom"].set_title("phi projection zoom")
        set_zoom_limits(axes["phi_zoom"], phis, phi_peak_index, args.zoom_half_width)
        axes["phi_zoom"].legend()

    title = args.title or "Fixed-kT sampler validation"
    title_lines = [
        title,
        f"{args.truth_label} integral ~= {truth_integral:.6g}, "
        f"weighted integral ~= {weighted_integral:.6g}, samples = {len(samples)}",
    ]
    if args.reference_integral is not None:
        reference_text = args.reference_label or "external reference"
        reference_text += f" ~= {args.reference_integral:.6g}"
        if args.reference_uncertainty is not None:
            reference_text += f" +/- {args.reference_uncertainty:.2g}"
        title_lines.append(reference_text)
    if args.warning:
        title_lines.append(args.warning)
    fig.suptitle(
        "\n".join(title_lines),
        fontsize=12,
    )
    fig.savefig(args.output_png, dpi=160)
    plt.close(fig)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
