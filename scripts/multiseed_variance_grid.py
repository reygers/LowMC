#!/usr/bin/env python3
"""Run a multi-seed p=100 comparison of mass and variance grid refinement."""

from __future__ import annotations

import argparse
import csv
import math
import statistics
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Variant:
    label: str
    refinement_score: str
    variance_probe_bins: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Multi-seed stability check for p=100 grid refinement modes."
    )
    parser.add_argument(
        "summary_csv",
        nargs="?",
        type=Path,
        default=Path("build/variance_grid_multiseed.csv"),
        help="per-seed summary CSV to write",
    )
    parser.add_argument(
        "--aggregate-csv",
        type=Path,
        default=Path("build/variance_grid_multiseed_aggregate.csv"),
        help="aggregate summary CSV to write",
    )
    parser.add_argument(
        "--summary-md",
        type=Path,
        default=Path("build/variance_grid_multiseed.md"),
        help="Markdown summary to write",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("build/variance_grid_multiseed"),
        help="directory for individual sampler-study CSV files",
    )
    parser.add_argument(
        "--study-exe",
        type=Path,
        default=Path("build/lowmc_weighted_sampler_study"),
        help="path to lowmc_weighted_sampler_study",
    )
    parser.add_argument("--samples", type=int, default=100000, help="samples per seed")
    parser.add_argument("--seeds", type=int, default=5, help="number of independent seeds")
    parser.add_argument("--seed", type=int, default=12345, help="first random seed")
    parser.add_argument("--momentum", type=float, default=100.0, help="toy pion momentum in GeV")
    parser.add_argument("--kt", type=float, default=0.02, help="fixed kT in GeV")
    parser.add_argument(
        "--reference",
        type=float,
        default=157.5401,
        help="high-stat weighted p=100 reference",
    )
    parser.add_argument(
        "--reference-error",
        type=float,
        default=0.319,
        help="standard error of the weighted reference",
    )
    return parser.parse_args()


def variants() -> list[Variant]:
    return [
        Variant("grid_mass_32x32_l3_r0.1_u0.1", "mass", 3),
        Variant("grid_variance_32x32_l3_r0.1_u0.1_p3", "variance", 3),
    ]


def run_command(command: list[str]) -> None:
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        if result.stdout:
            print(result.stdout, file=sys.stdout, end="")
        if result.stderr:
            print(result.stderr, file=sys.stderr, end="")
        raise subprocess.CalledProcessError(result.returncode, command)


def run_variant(args: argparse.Namespace, variant: Variant, seed: int) -> dict[str, str]:
    output_csv = args.output_dir / f"{variant.label}_seed{seed}.csv"
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
        "32",
        "--grid-phi-bins",
        "32",
        "--grid-refinement-levels",
        "3",
        "--grid-refine-fraction",
        "0.1",
        "--grid-uniform-fraction",
        "0.1",
        "--grid-refinement-score",
        variant.refinement_score,
        "--grid-variance-probe-bins",
        str(variant.variance_probe_bins),
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
    seed: int,
    reference: float,
    reference_error: float,
) -> dict[str, str]:
    samples = float(raw["samples"])
    mean = float(raw["meanWeight"])
    variance = float(raw["varianceWeight"])
    rms = math.sqrt(max(variance, 0.0))
    error = math.sqrt(variance / samples) if samples > 0.0 else math.nan
    total_error = math.sqrt(error * error + reference_error * reference_error)
    difference = mean - reference
    pull = difference / total_error if total_error > 0.0 else math.nan

    return {
        "variant": variant.label,
        "seed": str(seed),
        "gridRefinementScore": variant.refinement_score,
        "gridVarianceProbeBins": str(variant.variance_probe_bins),
        "samples": raw["samples"],
        "meanWeight": raw["meanWeight"],
        "weightRMS": f"{rms:.17g}",
        "weightedStandardError": f"{error:.17g}",
        "effectiveSampleSize": raw["effectiveSampleSize"],
        "effectiveSampleFraction": raw["effectiveSampleFraction"],
        "weightedReference": f"{reference:.17g}",
        "weightedReferenceError": f"{reference_error:.17g}",
        "referenceDifference": f"{difference:.17g}",
        "referencePull": f"{pull:.17g}",
        "csv": raw["csv"],
    }


def sample_stdev(values: list[float]) -> float:
    if len(values) < 2:
        return 0.0
    return statistics.stdev(values)


def aggregate_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    output: list[dict[str, str]] = []
    for variant in sorted({row["variant"] for row in rows}):
        group = [row for row in rows if row["variant"] == variant]
        mean_weights = [float(row["meanWeight"]) for row in group]
        rmss = [float(row["weightRMS"]) for row in group]
        ess_fractions = [float(row["effectiveSampleFraction"]) for row in group]
        pulls = [float(row["referencePull"]) for row in group]
        output.append(
            {
                "variant": variant,
                "gridRefinementScore": group[0]["gridRefinementScore"],
                "gridVarianceProbeBins": group[0]["gridVarianceProbeBins"],
                "seeds": str(len(group)),
                "meanWeightMean": f"{statistics.mean(mean_weights):.17g}",
                "meanWeightStd": f"{sample_stdev(mean_weights):.17g}",
                "weightRMSMean": f"{statistics.mean(rmss):.17g}",
                "weightRMSStd": f"{sample_stdev(rmss):.17g}",
                "effectiveSampleFractionMean": f"{statistics.mean(ess_fractions):.17g}",
                "effectiveSampleFractionStd": f"{sample_stdev(ess_fractions):.17g}",
                "referencePullMean": f"{statistics.mean(pulls):.17g}",
                "referencePullStd": f"{sample_stdev(pulls):.17g}",
            }
        )
    return output


def write_summary_md(path: Path, aggregate: list[dict[str, str]]) -> None:
    sorted_rows = sorted(
        aggregate,
        key=lambda row: float(row["effectiveSampleFractionMean"]),
        reverse=True,
    )
    lines = [
        "# Variance Grid Multi-Seed Summary",
        "",
        "| variant | seeds | mean weight | weight RMS | ESS/N | pull |",
        "| --- | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in sorted_rows:
        lines.append(
            "| "
            f"{row['variant']} | "
            f"{row['seeds']} | "
            f"{float(row['meanWeightMean']):.6g} +/- {float(row['meanWeightStd']):.3g} | "
            f"{float(row['weightRMSMean']):.6g} +/- {float(row['weightRMSStd']):.3g} | "
            f"{float(row['effectiveSampleFractionMean']):.6g} +/- "
            f"{float(row['effectiveSampleFractionStd']):.3g} | "
            f"{float(row['referencePullMean']):.3g} +/- "
            f"{float(row['referencePullStd']):.3g} |"
        )

    lines.extend(["", "## Interpretation", ""])
    best = sorted_rows[0]
    lines.append(
        f"- Highest mean ESS/N: `{best['variant']}` with "
        f"`ESS/N = {float(best['effectiveSampleFractionMean']):.4g}`."
    )
    if len(sorted_rows) >= 2:
        ratio = (
            float(sorted_rows[0]["effectiveSampleFractionMean"]) /
            float(sorted_rows[1]["effectiveSampleFractionMean"])
        )
        lines.append(f"- ESS/N ratio between the two variants is `{ratio:.4g}`.")
    lines.append(
        "- Treat this as a stochastic diagnostic; proposal-construction cost is not included."
    )
    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    args = parse_args()
    if args.samples <= 0:
        print("--samples must be positive", file=sys.stderr)
        return 1
    if args.seeds <= 0:
        print("--seeds must be positive", file=sys.stderr)
        return 1
    if args.reference_error < 0.0:
        print("--reference-error must be non-negative", file=sys.stderr)
        return 1

    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.summary_csv.parent.mkdir(parents=True, exist_ok=True)
    args.aggregate_csv.parent.mkdir(parents=True, exist_ok=True)
    args.summary_md.parent.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, str]] = []
    for variant in variants():
        for offset in range(args.seeds):
            seed = args.seed + offset
            raw = run_variant(args, variant, seed)
            row = summarize_row(raw, variant, seed, args.reference, args.reference_error)
            rows.append(row)
            print(
                f"{row['variant']} seed={seed}: "
                f"mean={float(row['meanWeight']):.6g} "
                f"RMS={float(row['weightRMS']):.6g} "
                f"ESS/N={float(row['effectiveSampleFraction']):.4g} "
                f"pull={float(row['referencePull']):.3g}"
            )

    fieldnames = [
        "variant",
        "seed",
        "gridRefinementScore",
        "gridVarianceProbeBins",
        "samples",
        "meanWeight",
        "weightRMS",
        "weightedStandardError",
        "effectiveSampleSize",
        "effectiveSampleFraction",
        "weightedReference",
        "weightedReferenceError",
        "referenceDifference",
        "referencePull",
        "csv",
    ]
    with args.summary_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    aggregate = aggregate_rows(rows)
    aggregate_fieldnames = [
        "variant",
        "gridRefinementScore",
        "gridVarianceProbeBins",
        "seeds",
        "meanWeightMean",
        "meanWeightStd",
        "weightRMSMean",
        "weightRMSStd",
        "effectiveSampleFractionMean",
        "effectiveSampleFractionStd",
        "referencePullMean",
        "referencePullStd",
    ]
    with args.aggregate_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=aggregate_fieldnames)
        writer.writeheader()
        writer.writerows(aggregate)
    write_summary_md(args.summary_md, aggregate)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
