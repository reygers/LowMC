#!/usr/bin/env python3
"""Run fixed-grid versus weighted-sample stress comparisons."""

from __future__ import annotations

import argparse
import csv
import math
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

from compare_fixed_grid_weighted import (
    DEFAULT_ALPHA,
    integrate_grid,
    read_weighted_sample_weights,
    sample_mean_and_error,
)


@dataclass(frozen=True)
class Case:
    name: str
    momentum: float
    particle_eta: float = 0.0
    particle_phi: float = 0.5 * math.pi
    eta_min: float = -5.0
    eta_max: float = 5.0
    phi_min: float = 0.0
    phi_max: float = 2.0 * math.pi


DEFAULT_CASES = [
    Case(name="mid_p1", momentum=1.0),
    Case(name="mid_p10", momentum=10.0),
    Case(name="mid_p100", momentum=100.0),
    Case(name="mid_p1000", momentum=1000.0),
    Case(
        name="forward_eta4_p100",
        momentum=100.0,
        particle_eta=4.0,
        eta_min=2.0,
        eta_max=5.0,
    ),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate fixed-grid and weighted-sample diagnostics for boosted "
            "two-particle toy cases, then summarize their agreement."
        )
    )
    parser.add_argument(
        "summary_csv",
        nargs="?",
        type=Path,
        default=Path("build/grid_weighted_stress_summary.csv"),
        help="summary CSV to write",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("build/grid_weighted_stress"),
        help="directory for per-case grid and weighted CSV files",
    )
    parser.add_argument(
        "--grid-exe",
        type=Path,
        default=Path("build/lowmc_two_particle_diagnostic"),
        help="path to lowmc_two_particle_diagnostic",
    )
    parser.add_argument(
        "--weighted-exe",
        type=Path,
        default=Path("build/lowmc_fixed_kt_weighted_samples"),
        help="path to lowmc_fixed_kt_weighted_samples",
    )
    parser.add_argument(
        "--samples",
        type=int,
        default=100000,
        help="weighted samples per case",
    )
    parser.add_argument("--seed", type=int, default=12345, help="base random seed")
    parser.add_argument(
        "--alpha",
        type=float,
        default=DEFAULT_ALPHA,
        help="fine-structure constant used for grid conversion",
    )
    parser.add_argument(
        "--case",
        choices=[case.name for case in DEFAULT_CASES],
        action="append",
        help="run only the named default case; may be repeated",
    )
    parser.add_argument(
        "--eta-bins",
        type=int,
        default=101,
        help="eta grid points passed to lowmc_two_particle_diagnostic",
    )
    parser.add_argument(
        "--phi-bins",
        type=int,
        default=121,
        help="phi grid points passed to lowmc_two_particle_diagnostic",
    )
    return parser.parse_args()


def run_command(command: list[str]) -> None:
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        if result.stdout:
            print(result.stdout, file=sys.stdout, end="")
        if result.stderr:
            print(result.stderr, file=sys.stderr, end="")
        raise subprocess.CalledProcessError(result.returncode, command)


def case_arguments(case: Case) -> list[str]:
    return [
        "--momentum",
        f"{case.momentum:.17g}",
        "--particle-eta",
        f"{case.particle_eta:.17g}",
        "--particle-phi",
        f"{case.particle_phi:.17g}",
        "--eta-min",
        f"{case.eta_min:.17g}",
        "--eta-max",
        f"{case.eta_max:.17g}",
        "--phi-min",
        f"{case.phi_min:.17g}",
        "--phi-max",
        f"{case.phi_max:.17g}",
    ]


def summarize_case(
    case: Case,
    grid_csv: Path,
    weighted_csv: Path,
    samples: int,
    alpha: float,
) -> dict[str, float | str | int]:
    grid_estimate = integrate_grid(grid_csv, alpha)
    weights = read_weighted_sample_weights(weighted_csv)
    mean, error = sample_mean_and_error(weights)
    difference = mean - grid_estimate
    relative_difference = difference / grid_estimate if grid_estimate != 0.0 else math.nan
    pull = difference / error if error > 0.0 else math.nan

    sum_weights = sum(weights)
    sum_weight_squares = sum(weight * weight for weight in weights)
    effective_sample_size = (
        (sum_weights * sum_weights) / sum_weight_squares if sum_weight_squares > 0.0 else 0.0
    )
    effective_sample_fraction = effective_sample_size / len(weights)

    return {
        "case": case.name,
        "momentum": case.momentum,
        "particleEta": case.particle_eta,
        "particlePhi": case.particle_phi,
        "etaMin": case.eta_min,
        "etaMax": case.eta_max,
        "phiMin": case.phi_min,
        "phiMax": case.phi_max,
        "samples": samples,
        "gridIntegral": grid_estimate,
        "weightedMean": mean,
        "weightedStandardError": error,
        "difference": difference,
        "relativeDifference": relative_difference,
        "pull": pull,
        "effectiveSampleSize": effective_sample_size,
        "effectiveSampleFraction": effective_sample_fraction,
        "gridCsv": str(grid_csv),
        "weightedCsv": str(weighted_csv),
    }


def main() -> int:
    args = parse_args()
    if args.samples <= 0:
        print("--samples must be positive", file=sys.stderr)
        return 1
    if args.eta_bins < 2 or args.phi_bins < 2:
        print("--eta-bins and --phi-bins must be >= 2", file=sys.stderr)
        return 1

    selected_cases = DEFAULT_CASES
    if args.case:
        selected = set(args.case)
        selected_cases = [case for case in DEFAULT_CASES if case.name in selected]

    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.summary_csv.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    for index, case in enumerate(selected_cases):
        grid_csv = args.output_dir / f"{case.name}_grid.csv"
        weighted_csv = args.output_dir / f"{case.name}_weighted.csv"

        run_command(
            [
                str(args.grid_exe),
                str(grid_csv),
                "--eta-bins",
                str(args.eta_bins),
                "--phi-bins",
                str(args.phi_bins),
                *case_arguments(case),
            ]
        )
        run_command(
            [
                str(args.weighted_exe),
                str(weighted_csv),
                "--samples",
                str(args.samples),
                "--seed",
                str(args.seed + index),
                *case_arguments(case),
            ]
        )

        row = summarize_case(case, grid_csv, weighted_csv, args.samples, args.alpha)
        rows.append(row)
        print(
            f"{case.name}: grid={row['gridIntegral']:.6g} "
            f"weighted={row['weightedMean']:.6g} +/- {row['weightedStandardError']:.3g} "
            f"pull={row['pull']:.3g} ESS/N={row['effectiveSampleFraction']:.3g}"
        )

    fieldnames = [
        "case",
        "momentum",
        "particleEta",
        "particlePhi",
        "etaMin",
        "etaMax",
        "phiMin",
        "phiMax",
        "etaBins",
        "phiBins",
        "samples",
        "gridIntegral",
        "weightedMean",
        "weightedStandardError",
        "difference",
        "relativeDifference",
        "pull",
        "effectiveSampleSize",
        "effectiveSampleFraction",
        "gridCsv",
        "weightedCsv",
    ]
    with args.summary_csv.open("w", newline="") as handle:
        for row in rows:
            row["etaBins"] = args.eta_bins
            row["phiBins"] = args.phi_bins
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
