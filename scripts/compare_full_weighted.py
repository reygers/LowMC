#!/usr/bin/env python3
"""Compare full phase-space uniform and grid-adaptive weighted sampling."""

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
    proposal: str
    options: tuple[str, ...]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare full log-uniform-kT weighted sampling paths."
    )
    parser.add_argument(
        "summary_csv",
        nargs="?",
        type=Path,
        default=Path("build/full_weighted_comparison.csv"),
        help="summary CSV to write",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("build/full_weighted_comparison"),
        help="directory for individual sampler-study CSV files",
    )
    parser.add_argument(
        "--study-exe",
        type=Path,
        default=Path("build/lowmc_full_weighted_sampler_study"),
        help="path to lowmc_full_weighted_sampler_study",
    )
    parser.add_argument("--samples", type=int, default=100000, help="samples per variant")
    parser.add_argument("--seed", type=int, default=12345, help="base random seed")
    parser.add_argument("--momentum", type=float, default=100.0, help="toy pion momentum in GeV")
    parser.add_argument("--kt-min", type=float, default=0.002, help="minimum kT in GeV")
    parser.add_argument("--kt-max", type=float, default=0.2, help="maximum kT in GeV")
    parser.add_argument(
        "--fixed-kt-angular-reference",
        type=float,
        default=157.5401,
        help="weighted angular reference at fixed kT=0.02",
    )
    parser.add_argument(
        "--fixed-kt-reference-kt",
        type=float,
        default=0.02,
        help="kT value of the fixed angular reference",
    )
    parser.add_argument(
        "--summary-md",
        type=Path,
        default=Path("build/full_weighted_comparison.md"),
        help="Markdown summary to write",
    )
    return parser.parse_args()


def variants() -> list[Variant]:
    return [
        Variant("full_uniform", "uniform", ()),
        Variant(
            "full_grid_32x32_l3_r0.1_u0.1",
            "grid",
            (
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
            ),
        ),
    ]


def run_command(command: list[str]) -> None:
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        if result.stdout:
            print(result.stdout, file=sys.stdout, end="")
        if result.stderr:
            print(result.stderr, file=sys.stderr, end="")
        raise subprocess.CalledProcessError(result.returncode, command)


def run_variant(args: argparse.Namespace, variant: Variant) -> dict[str, str]:
    output_csv = args.output_dir / f"{variant.label}.csv"
    command = [
        str(args.study_exe),
        str(output_csv),
        "--samples",
        str(args.samples),
        "--seed",
        str(args.seed),
        "--momentum",
        f"{args.momentum:.17g}",
        "--kt-min",
        f"{args.kt_min:.17g}",
        "--kt-max",
        f"{args.kt_max:.17g}",
        "--proposal",
        variant.proposal,
        *variant.options,
    ]
    run_command(command)
    with output_csv.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 1:
        raise RuntimeError(f"expected one row in {output_csv}, got {len(rows)}")
    rows[0]["csv"] = str(output_csv)
    return rows[0]


def expected_full_reference(args: argparse.Namespace) -> float:
    return (
        args.fixed_kt_angular_reference *
        args.fixed_kt_reference_kt *
        args.fixed_kt_reference_kt *
        math.log(args.kt_max / args.kt_min)
    )


def summarize_row(raw: dict[str, str], label: str, reference: float) -> dict[str, str]:
    samples = float(raw["samples"])
    mean = float(raw["meanWeight"])
    variance = float(raw["varianceWeight"])
    rms = math.sqrt(max(variance, 0.0))
    error = math.sqrt(variance / samples) if samples > 0.0 else math.nan
    difference = mean - reference
    pull = difference / error if error > 0.0 else math.nan

    return {
        "variant": label,
        "proposal": raw["proposal"],
        "momentum": raw["momentum"],
        "kTMin": raw["kTMin"],
        "kTMax": raw["kTMax"],
        "gridBuildKT": raw["gridBuildKT"],
        "samples": raw["samples"],
        "meanWeight": raw["meanWeight"],
        "weightRMS": f"{rms:.17g}",
        "weightedStandardError": f"{error:.17g}",
        "effectiveSampleSize": raw["effectiveSampleSize"],
        "effectiveSampleFraction": raw["effectiveSampleFraction"],
        "expectedFullReference": f"{reference:.17g}",
        "referenceDifference": f"{difference:.17g}",
        "referencePull": f"{pull:.17g}",
        "csv": raw["csv"],
    }


def write_summary_md(path: Path, rows: list[dict[str, str]]) -> None:
    sorted_rows = sorted(
        rows,
        key=lambda row: float(row["effectiveSampleFraction"]),
        reverse=True,
    )
    reference = float(rows[0]["expectedFullReference"])
    lines = [
        "# Full Weighted Sampling Comparison",
        "",
        f"Expected full-reference estimate from fixed-kT scaling: `{reference:.6g}`.",
        "",
        "| variant | mean weight | weight RMS | ESS/N | pull vs expected |",
        "| --- | ---: | ---: | ---: | ---: |",
    ]
    for row in sorted_rows:
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
        f"- Highest ESS/N: `{best['variant']}` with "
        f"`ESS/N = {float(best['effectiveSampleFraction']):.4g}`."
    )
    if len(sorted_rows) >= 2:
        ratio = (
            float(sorted_rows[0]["effectiveSampleFraction"]) /
            float(sorted_rows[1]["effectiveSampleFraction"])
        )
        lines.append(f"- ESS/N ratio between best and next variant is `{ratio:.4g}`.")
    lines.append(
        "- The grid proposal is built once at the geometric-mean kT and exact Low "
        "weights are evaluated at every sampled kT."
    )
    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    args = parse_args()
    if args.samples <= 0:
        print("--samples must be positive", file=sys.stderr)
        return 1
    if args.kt_min <= 0.0 or args.kt_max <= args.kt_min:
        print("kT bounds must satisfy 0 < kt-min < kt-max", file=sys.stderr)
        return 1

    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.summary_csv.parent.mkdir(parents=True, exist_ok=True)
    args.summary_md.parent.mkdir(parents=True, exist_ok=True)

    reference = expected_full_reference(args)
    rows: list[dict[str, str]] = []
    for variant in variants():
        raw = run_variant(args, variant)
        row = summarize_row(raw, variant.label, reference)
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
        "momentum",
        "kTMin",
        "kTMax",
        "gridBuildKT",
        "samples",
        "meanWeight",
        "weightRMS",
        "weightedStandardError",
        "effectiveSampleSize",
        "effectiveSampleFraction",
        "expectedFullReference",
        "referenceDifference",
        "referencePull",
        "csv",
    ]
    with args.summary_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    write_summary_md(args.summary_md, rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
