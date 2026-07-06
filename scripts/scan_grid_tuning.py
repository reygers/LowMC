#!/usr/bin/env python3
"""Scan grid-adaptive eta-phi proposal parameters for the p=100 toy case."""

from __future__ import annotations

import argparse
import csv
import math
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Variant:
    label: str
    eta_bins: int
    phi_bins: int
    refinement_levels: int
    refine_fraction: float
    uniform_fraction: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Scan grid-adaptive proposal parameters against a fixed grid reference."
    )
    parser.add_argument(
        "summary_csv",
        nargs="?",
        type=Path,
        default=Path("build/grid_tuning_summary.csv"),
        help="summary CSV to write",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("build/grid_tuning"),
        help="directory for intermediate sampler-study CSV files",
    )
    parser.add_argument(
        "--study-exe",
        type=Path,
        default=Path("build/lowmc_weighted_sampler_study"),
        help="path to lowmc_weighted_sampler_study",
    )
    parser.add_argument("--samples", type=int, default=100000, help="samples per variant")
    parser.add_argument("--seed", type=int, default=12345, help="base random seed")
    parser.add_argument("--momentum", type=float, default=100.0, help="toy pion momentum in GeV")
    parser.add_argument("--kt", type=float, default=0.02, help="fixed kT in GeV")
    parser.add_argument(
        "--reference",
        type=float,
        default=148.7805831520777,
        help="best available p=100 grid reference integral",
    )
    parser.add_argument(
        "--summary-md",
        type=Path,
        default=Path("build/grid_tuning_summary.md"),
        help="Markdown summary to write",
    )
    parser.add_argument(
        "--max-pull",
        type=float,
        default=3.0,
        help="maximum |pull| for stable-variant selection",
    )
    return parser.parse_args()


def default_variants() -> list[Variant]:
    variants: list[Variant] = []
    for eta_bins, phi_bins in [(16, 16), (24, 24), (32, 32)]:
        for refinement_levels in [2, 3]:
            for refine_fraction in [0.05, 0.1, 0.2]:
                for uniform_fraction in [0.1, 0.2, 0.4]:
                    variants.append(
                        Variant(
                            label=(
                                f"grid_{eta_bins}x{phi_bins}_l{refinement_levels}"
                                f"_r{refine_fraction:g}_u{uniform_fraction:g}"
                            ),
                            eta_bins=eta_bins,
                            phi_bins=phi_bins,
                            refinement_levels=refinement_levels,
                            refine_fraction=refine_fraction,
                            uniform_fraction=uniform_fraction,
                        )
                    )
    return variants


def run_command(command: list[str]) -> None:
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        if result.stdout:
            print(result.stdout, file=sys.stdout, end="")
        if result.stderr:
            print(result.stderr, file=sys.stderr, end="")
        raise subprocess.CalledProcessError(result.returncode, command)


def run_variant(args: argparse.Namespace, variant: Variant, seed: int) -> dict[str, str]:
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
        "--momentum",
        f"{args.momentum:.17g}",
        "--proposal",
        "grid",
        "--grid-eta-bins",
        str(variant.eta_bins),
        "--grid-phi-bins",
        str(variant.phi_bins),
        "--grid-refinement-levels",
        str(variant.refinement_levels),
        "--grid-refine-fraction",
        f"{variant.refine_fraction:.17g}",
        "--grid-uniform-fraction",
        f"{variant.uniform_fraction:.17g}",
    ]
    run_command(command)
    with output_csv.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 1:
        raise RuntimeError(f"expected one row in {output_csv}, got {len(rows)}")
    rows[0]["csv"] = str(output_csv)
    return rows[0]


def summarize_row(
    raw: dict[str, str],
    variant: Variant,
    reference: float,
) -> dict[str, str]:
    samples = float(raw["samples"])
    mean = float(raw["meanWeight"])
    variance = float(raw["varianceWeight"])
    rms = math.sqrt(max(variance, 0.0))
    error = math.sqrt(variance / samples) if samples > 0.0 else math.nan
    difference = mean - reference
    pull = difference / error if error > 0.0 else math.nan

    return {
        "variant": variant.label,
        "proposal": "grid",
        "gridEtaBins": str(variant.eta_bins),
        "gridPhiBins": str(variant.phi_bins),
        "gridRefinementLevels": str(variant.refinement_levels),
        "gridRefineFraction": f"{variant.refine_fraction:.17g}",
        "gridUniformFraction": f"{variant.uniform_fraction:.17g}",
        "momentum": raw["momentum"],
        "kT": raw["kT"],
        "samples": raw["samples"],
        "meanWeight": raw["meanWeight"],
        "weightRMS": f"{rms:.17g}",
        "weightedStandardError": f"{error:.17g}",
        "effectiveSampleSize": raw["effectiveSampleSize"],
        "effectiveSampleFraction": raw["effectiveSampleFraction"],
        "referenceGridIntegral": f"{reference:.17g}",
        "referenceDifference": f"{difference:.17g}",
        "referencePull": f"{pull:.17g}",
        "csv": raw["csv"],
    }


def write_summary_md(path: Path, rows: list[dict[str, str]], max_pull: float) -> None:
    sorted_rows = sorted(
        rows,
        key=lambda row: float(row["effectiveSampleFraction"]),
        reverse=True,
    )
    stable_rows = [
        row for row in sorted_rows
        if math.isfinite(float(row["referencePull"])) and abs(float(row["referencePull"])) <= max_pull
    ]

    lines = [
        "# Grid-Adaptive Tuning Summary",
        "",
        f"Stable variants use `|pull| <= {max_pull:g}` against the current grid reference.",
        "",
        "| variant | mean weight | weight RMS | ESS/N | pull |",
        "| --- | ---: | ---: | ---: | ---: |",
    ]
    for row in sorted_rows[:12]:
        lines.append(
            "| "
            f"{row['variant']} | "
            f"{float(row['meanWeight']):.6g} | "
            f"{float(row['weightRMS']):.6g} | "
            f"{float(row['effectiveSampleFraction']):.6g} | "
            f"{float(row['referencePull']):.3g} |"
        )

    lines.extend(["", "## Interpretation", ""])
    best = sorted_rows[0]
    lines.append(
        f"- Highest ESS/N tested: `{best['variant']}` with "
        f"`ESS/N = {float(best['effectiveSampleFraction']):.4g}` and "
        f"`pull = {float(best['referencePull']):.3g}`."
    )
    if stable_rows:
        stable = stable_rows[0]
        lines.append(
            f"- Highest ESS/N with acceptable pull: `{stable['variant']}` with "
            f"`ESS/N = {float(stable['effectiveSampleFraction']):.4g}` and "
            f"`pull = {float(stable['referencePull']):.3g}`."
        )
    else:
        lines.append("- No scanned grid-adaptive variant satisfies the pull cut.")
    lines.append(
        "- Treat this as diagnostic tuning: the p=100 reference still has visible "
        "grid-convergence uncertainty."
    )
    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    args = parse_args()
    if args.samples <= 0:
        print("--samples must be positive", file=sys.stderr)
        return 1

    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.summary_csv.parent.mkdir(parents=True, exist_ok=True)
    args.summary_md.parent.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, str]] = []
    for index, variant in enumerate(default_variants()):
        raw = run_variant(args, variant, args.seed + index)
        row = summarize_row(raw, variant, args.reference)
        rows.append(row)
        print(
            f"{row['variant']}: mean={float(row['meanWeight']):.6g} "
            f"RMS={float(row['weightRMS']):.6g} "
            f"ESS/N={float(row['effectiveSampleFraction']):.4g} "
            f"pull={float(row['referencePull']):.3g}"
        )

    fieldnames = [
        "variant",
        "proposal",
        "gridEtaBins",
        "gridPhiBins",
        "gridRefinementLevels",
        "gridRefineFraction",
        "gridUniformFraction",
        "momentum",
        "kT",
        "samples",
        "meanWeight",
        "weightRMS",
        "weightedStandardError",
        "effectiveSampleSize",
        "effectiveSampleFraction",
        "referenceGridIntegral",
        "referenceDifference",
        "referencePull",
        "csv",
    ]
    with args.summary_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    write_summary_md(args.summary_md, rows, args.max_pull)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
