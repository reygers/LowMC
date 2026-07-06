#!/usr/bin/env python3
"""Scan adaptive fixed-kT eta-phi proposal parameters."""

from __future__ import annotations

import argparse
import csv
import math
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

from compare_fixed_grid_weighted import DEFAULT_ALPHA, integrate_grid


@dataclass(frozen=True)
class Variant:
    label: str
    proposal: str
    uniform_fraction: float
    width_scale: float
    min_half_width: float
    max_half_width: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run uniform and adaptive weighted-sampler studies for a small "
            "parameter scan and write one comparison CSV."
        )
    )
    parser.add_argument(
        "summary_csv",
        nargs="?",
        type=Path,
        default=Path("build/adaptive_tuning_summary.csv"),
        help="summary CSV to write",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("build/adaptive_tuning"),
        help="directory for intermediate CSV files",
    )
    parser.add_argument(
        "--study-exe",
        type=Path,
        default=Path("build/lowmc_weighted_sampler_study"),
        help="path to lowmc_weighted_sampler_study",
    )
    parser.add_argument(
        "--grid-exe",
        type=Path,
        default=Path("build/lowmc_two_particle_diagnostic"),
        help="path to lowmc_two_particle_diagnostic",
    )
    parser.add_argument("--samples", type=int, default=50000, help="samples per variant")
    parser.add_argument("--seed", type=int, default=12345, help="base random seed")
    parser.add_argument("--kt", type=float, default=0.02, help="fixed kT in GeV")
    parser.add_argument(
        "--momentum",
        type=float,
        action="append",
        help="charged-pion momentum to include; may be repeated",
    )
    parser.add_argument(
        "--uniform-fraction",
        type=float,
        default=0.2,
        help="adaptive uniform floor fraction",
    )
    parser.add_argument(
        "--width-scale",
        type=float,
        action="append",
        help="adaptive width scale to scan; may be repeated",
    )
    parser.add_argument(
        "--min-half-width",
        type=float,
        action="append",
        help="adaptive minimum half-width to scan; may be repeated",
    )
    parser.add_argument(
        "--max-half-width",
        type=float,
        default=0.8,
        help="adaptive maximum half-width",
    )
    parser.add_argument(
        "--eta-bins",
        type=int,
        default=401,
        help="eta grid points for fine-grid references",
    )
    parser.add_argument(
        "--phi-bins",
        type=int,
        default=481,
        help="phi grid points for fine-grid references",
    )
    parser.add_argument(
        "--skip-reference",
        action="store_true",
        help="do not compute fine-grid reference integrals",
    )
    parser.add_argument(
        "--alpha",
        type=float,
        default=DEFAULT_ALPHA,
        help="fine-structure constant used for grid conversion",
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


def default_momenta(args: argparse.Namespace) -> list[float]:
    return args.momentum if args.momentum else [10.0, 100.0, 1000.0]


def scan_variants(args: argparse.Namespace) -> list[Variant]:
    width_scales = args.width_scale if args.width_scale else [2.0, 5.0, 10.0, 20.0]
    min_half_widths = args.min_half_width if args.min_half_width else [0.01, 0.02, 0.05]

    variants = [
        Variant(
            label="uniform",
            proposal="uniform",
            uniform_fraction=args.uniform_fraction,
            width_scale=5.0,
            min_half_width=0.02,
            max_half_width=args.max_half_width,
        )
    ]
    for width_scale in width_scales:
        for min_half_width in min_half_widths:
            variants.append(
                Variant(
                    label=f"adaptive_w{width_scale:g}_min{min_half_width:g}",
                    proposal="adaptive",
                    uniform_fraction=args.uniform_fraction,
                    width_scale=width_scale,
                    min_half_width=min_half_width,
                    max_half_width=args.max_half_width,
                )
            )

    return variants


def momentum_options(momenta: list[float]) -> list[str]:
    options: list[str] = []
    for momentum in momenta:
        options.extend(["--momentum", f"{momentum:.17g}"])
    return options


def compute_references(args: argparse.Namespace, momenta: list[float]) -> dict[float, float]:
    references: dict[float, float] = {}
    if args.skip_reference:
        return references

    reference_dir = args.output_dir / "references"
    reference_dir.mkdir(parents=True, exist_ok=True)
    for momentum in momenta:
        grid_csv = reference_dir / f"grid_p{momentum:g}.csv"
        run_command(
            [
                str(args.grid_exe),
                str(grid_csv),
                "--momentum",
                f"{momentum:.17g}",
                "--eta-bins",
                str(args.eta_bins),
                "--phi-bins",
                str(args.phi_bins),
            ]
        )
        references[momentum] = integrate_grid(grid_csv, args.alpha)

    return references


def run_variant(
    args: argparse.Namespace,
    variant: Variant,
    momenta: list[float],
    seed: int,
) -> list[dict[str, str]]:
    output_csv = args.output_dir / f"{variant.label}.csv"
    command = [
        str(args.study_exe),
        str(output_csv),
        "--samples",
        str(args.samples),
        "--seed",
        str(seed),
        "--kt",
        f"{args.kt:.17g}",
        "--proposal",
        variant.proposal,
        "--adaptive-uniform-fraction",
        f"{variant.uniform_fraction:.17g}",
        "--adaptive-width-scale",
        f"{variant.width_scale:.17g}",
        "--adaptive-min-half-width",
        f"{variant.min_half_width:.17g}",
        "--adaptive-max-half-width",
        f"{variant.max_half_width:.17g}",
        *momentum_options(momenta),
    ]
    run_command(command)

    with output_csv.open(newline="") as handle:
        reader = csv.DictReader(handle)
        return list(reader)


def augment_row(
    row: dict[str, str],
    variant: Variant,
    references: dict[float, float],
) -> dict[str, str]:
    samples = float(row["samples"])
    momentum = float(row["momentum"])
    variance = float(row["varianceWeight"])
    mean = float(row["meanWeight"])
    rms = math.sqrt(max(variance, 0.0))
    standard_error = math.sqrt(variance / samples) if samples > 0.0 else math.nan
    reference = references.get(momentum, math.nan)
    difference = mean - reference if math.isfinite(reference) else math.nan
    pull = difference / standard_error if standard_error > 0.0 else math.nan

    output = dict(row)
    output["variant"] = variant.label
    output["weightRMS"] = f"{rms:.17g}"
    output["weightedStandardError"] = f"{standard_error:.17g}"
    output["referenceGridIntegral"] = f"{reference:.17g}"
    output["referenceDifference"] = f"{difference:.17g}"
    output["referencePull"] = f"{pull:.17g}"
    return output


def main() -> int:
    args = parse_args()
    if args.samples <= 0:
        print("--samples must be positive", file=sys.stderr)
        return 1
    if args.eta_bins < 2 or args.phi_bins < 2:
        print("--eta-bins and --phi-bins must be >= 2", file=sys.stderr)
        return 1

    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.summary_csv.parent.mkdir(parents=True, exist_ok=True)

    momenta = default_momenta(args)
    variants = scan_variants(args)
    references = compute_references(args, momenta)

    rows: list[dict[str, str]] = []
    for index, variant in enumerate(variants):
        variant_rows = run_variant(args, variant, momenta, args.seed + index)
        for row in variant_rows:
            augmented = augment_row(row, variant, references)
            rows.append(augmented)
            print(
                f"{augmented['variant']} p={float(augmented['momentum']):g}: "
                f"mean={float(augmented['meanWeight']):.6g} "
                f"RMS={float(augmented['weightRMS']):.6g} "
                f"ESS/N={float(augmented['effectiveSampleFraction']):.4g} "
                f"pull={float(augmented['referencePull']):.3g}"
            )

    fieldnames = [
        "variant",
        "proposal",
        "adaptiveUniformFraction",
        "adaptiveWidthScale",
        "adaptiveMinHalfWidth",
        "adaptiveMaxHalfWidth",
        "momentum",
        "kT",
        "samples",
        "sumWeight",
        "meanWeight",
        "minWeight",
        "maxWeight",
        "varianceWeight",
        "weightRMS",
        "weightedStandardError",
        "effectiveSampleSize",
        "effectiveSampleFraction",
        "referenceGridIntegral",
        "referenceDifference",
        "referencePull",
    ]
    with args.summary_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
