#!/usr/bin/env python3
"""Plot normalized PHOTOS and weighted-LowMC Z-to-mumu photon shapes."""

from __future__ import annotations

import argparse
import csv
import math
import os
import sys
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Observable:
    name: str
    label: str
    low: float
    high: float
    bins: int


OBSERVABLES = (
    Observable("kT", r"$k_T$ [GeV]", 0.05, 0.2, 24),
    Observable("eta", r"$\eta$", -2.0, 2.0, 24),
    Observable("phi", r"$\phi$", -math.pi, math.pi, 24),
    Observable("theta_nearest_born", r"angle to nearest Born muon [rad]", 0.0, 0.5 * math.pi, 24),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot PHOTOS versus weighted-LowMC shapes for the fixed Z-to-mumu pilot."
    )
    parser.add_argument("photos_csv", type=Path)
    parser.add_argument("lowmc_csv", type=Path)
    parser.add_argument("output_png", type=Path)
    parser.add_argument(
        "--bins-csv",
        type=Path,
        default=Path("photos_zmumu_shape_bins.csv"),
        help="machine-readable binned comparison output",
    )
    parser.add_argument(
        "--summary-csv",
        type=Path,
        default=Path("photos_zmumu_shape_summary.csv"),
        help="per-observable approximate chi-square summary",
    )
    return parser.parse_args()


def read_photos(path: Path) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            if int(row["accepted"]) != 1:
                continue
            theta = min(float(row["theta_plus_born"]), float(row["theta_minus_born"]))
            rows.append(
                {
                    "kT": float(row["kT"]),
                    "eta": float(row["eta"]),
                    "phi": float(row["phi"]),
                    "theta_nearest_born": theta,
                    "weight": 1.0,
                }
            )
    if not rows:
        raise ValueError("PHOTOS CSV has no accepted photons")
    return rows


def read_lowmc(path: Path) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            rows.append(
                {
                    "kT": float(row["kT"]),
                    "eta": float(row["eta"]),
                    "phi": float(row["phi"]),
                    "theta_nearest_born": float(row["theta_nearest_born"]),
                    "weight": float(row["weight"]),
                }
            )
    if not rows:
        raise ValueError("LowMC CSV is empty")
    return rows


def normalized_histogram(
    rows: list[dict[str, float]], observable: Observable, *, poisson_floor: bool = False
) -> tuple[list[float], list[float], list[float]]:
    width = (observable.high - observable.low) / observable.bins
    sum_w = [0.0] * observable.bins
    sum_w2 = [0.0] * observable.bins
    total = 0.0
    for row in rows:
        value = row[observable.name]
        if not (observable.low <= value < observable.high):
            continue
        index = min(observable.bins - 1, int((value - observable.low) / width))
        weight = row["weight"]
        sum_w[index] += weight
        sum_w2[index] += weight * weight
        total += weight
    if total <= 0.0:
        raise ValueError(f"no positive weight for {observable.name}")
    centers = [observable.low + (index + 0.5) * width for index in range(observable.bins)]
    density = [value / (total * width) for value in sum_w]
    error = [
        math.sqrt(value + (0.25 if poisson_floor else 0.0)) / (total * width)
        for value in sum_w2
    ]
    return centers, density, error


def ratios(
    photos: list[float],
    photos_error: list[float],
    lowmc: list[float],
    lowmc_error: list[float],
) -> tuple[list[float], list[float]]:
    values: list[float] = []
    errors: list[float] = []
    for p_value, p_error, l_value, l_error in zip(
        photos, photos_error, lowmc, lowmc_error
    ):
        if p_value <= 0.0 or l_value <= 0.0:
            values.append(math.nan)
            errors.append(math.nan)
            continue
        ratio = p_value / l_value
        error = ratio * math.sqrt((p_error / p_value) ** 2 + (l_error / l_value) ** 2)
        values.append(ratio)
        errors.append(error)
    return values, errors


def main() -> int:
    args = parse_args()
    cache_directory = args.output_png.parent / ".matplotlib"
    cache_directory.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache_directory))
    os.environ.setdefault("XDG_CACHE_HOME", str(cache_directory))
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is required for plotting", file=sys.stderr)
        return 2

    try:
        photos_rows = read_photos(args.photos_csv)
        lowmc_rows = read_lowmc(args.lowmc_csv)
    except (OSError, KeyError, ValueError) as error:
        print(f"Could not read benchmark samples: {error}", file=sys.stderr)
        return 1

    figure, axes = plt.subplots(
        len(OBSERVABLES), 2, figsize=(11.0, 12.0), constrained_layout=True
    )
    bin_rows: list[dict[str, float | str]] = []
    summary_rows: list[dict[str, float | int | str]] = []
    for row_index, observable in enumerate(OBSERVABLES):
        centers, photos_density, photos_error = normalized_histogram(
            photos_rows, observable, poisson_floor=True
        )
        _, lowmc_density, lowmc_error = normalized_histogram(lowmc_rows, observable)
        ratio, ratio_error = ratios(
            photos_density, photos_error, lowmc_density, lowmc_error
        )
        chi_square = 0.0
        bins_used = 0
        for p_value, p_error, l_value, l_error in zip(
            photos_density, photos_error, lowmc_density, lowmc_error
        ):
            variance = p_error * p_error + l_error * l_error
            if variance <= 0.0:
                continue
            chi_square += (p_value - l_value) ** 2 / variance
            bins_used += 1
        ndf = max(1, bins_used - 1)
        summary_rows.append(
            {
                "observable": observable.name,
                "bins_used": bins_used,
                "chi_square": chi_square,
                "ndf": ndf,
                "chi_square_per_ndf": chi_square / ndf,
            }
        )
        width = (observable.high - observable.low) / observable.bins

        shape_axis = axes[row_index][0]
        ratio_axis = axes[row_index][1]
        shape_axis.errorbar(
            centers,
            photos_density,
            yerr=photos_error,
            fmt="o",
            markersize=3.0,
            color="#2563eb",
            label=f"PHOTOS ({len(photos_rows)} accepted)",
        )
        shape_axis.step(
            centers,
            lowmc_density,
            where="mid",
            color="#dc2626",
            linewidth=1.5,
            label=f"LowMC weighted ({len(lowmc_rows)} trials)",
        )
        shape_axis.set_ylabel("normalized density")
        shape_axis.set_xlabel(observable.label)
        shape_axis.set_xlim(observable.low, observable.high)
        shape_axis.set_ylim(bottom=0.0)
        if row_index == 0:
            shape_axis.legend(frameon=False)

        ratio_axis.errorbar(
            centers,
            ratio,
            yerr=ratio_error,
            fmt="o",
            markersize=3.0,
            color="#374151",
        )
        ratio_axis.axhline(1.0, color="#dc2626", linestyle="--", linewidth=1.0)
        ratio_axis.set_ylabel("PHOTOS / LowMC")
        ratio_axis.set_xlabel(observable.label)
        ratio_axis.set_xlim(observable.low, observable.high)

        for index, center in enumerate(centers):
            bin_rows.append(
                {
                    "observable": observable.name,
                    "bin_low": center - 0.5 * width,
                    "bin_high": center + 0.5 * width,
                    "photos_density": photos_density[index],
                    "photos_error": photos_error[index],
                    "lowmc_density": lowmc_density[index],
                    "lowmc_error": lowmc_error[index],
                    "photos_over_lowmc": ratio[index],
                    "ratio_error": ratio_error[index],
                }
            )

    figure.suptitle(
        r"Fixed $Z\to\mu^+\mu^-$ soft-photon shape pilot: PHOTOS 3.64 vs LowMC"
    )
    figure.savefig(args.output_png, dpi=180)
    plt.close(figure)

    args.bins_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.bins_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(bin_rows[0]))
        writer.writeheader()
        writer.writerows(bin_rows)
    args.summary_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.summary_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0]))
        writer.writeheader()
        writer.writerows(summary_rows)
    for row in summary_rows:
        print(
            f"observable={row['observable']} chi2_ndf={row['chi_square_per_ndf']:.3f} "
            f"chi2={row['chi_square']:.3f} ndf={row['ndf']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
