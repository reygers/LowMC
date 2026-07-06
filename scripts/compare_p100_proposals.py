#!/usr/bin/env python3
"""Compare p=100 fixed-kT eta-phi proposals against the weighted reference."""

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
        description="Compare tuned grid, local-box, and smooth-kernel p=100 proposals."
    )
    parser.add_argument(
        "summary_csv",
        nargs="?",
        type=Path,
        default=Path("build/p100_proposal_comparison.csv"),
        help="summary CSV to write",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("build/p100_proposal_comparison"),
        help="directory for individual sampler-study CSV files",
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
        default=157.5401,
        help="high-stat weighted p=100 reference",
    )
    parser.add_argument(
        "--reference-error",
        type=float,
        default=0.319,
        help="standard error of the weighted reference",
    )
    parser.add_argument(
        "--fixed-grid-reference",
        type=float,
        default=148.7805831520777,
        help="finest fixed-grid lower cross-check",
    )
    parser.add_argument(
        "--summary-md",
        type=Path,
        default=Path("build/p100_proposal_comparison.md"),
        help="Markdown summary to write",
    )
    return parser.parse_args()


def default_variants() -> list[Variant]:
    return [
        Variant(
            "grid_32x32_l3_r0.1_u0.1",
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
        Variant(
            "local_box_w40_min0.1_u0.2",
            "adaptive",
            (
                "--adaptive-uniform-fraction",
                "0.2",
                "--adaptive-width-scale",
                "40",
                "--adaptive-min-half-width",
                "0.1",
                "--adaptive-max-half-width",
                "0.8",
            ),
        ),
        Variant(
            "smooth_triangular_w40_min0.1_u0.2",
            "smooth",
            (
                "--smooth-uniform-fraction",
                "0.2",
                "--smooth-width-scale",
                "40",
                "--smooth-min-half-width",
                "0.1",
                "--smooth-max-half-width",
                "0.8",
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


def summarize_row(
    raw: dict[str, str],
    variant: Variant,
    reference: float,
    reference_error: float,
    fixed_grid_reference: float,
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
        "proposal": variant.proposal,
        "momentum": raw["momentum"],
        "kT": raw["kT"],
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
        "fixedGridReference": f"{fixed_grid_reference:.17g}",
        "fixedGridDifference": f"{mean - fixed_grid_reference:.17g}",
        "csv": raw["csv"],
    }


def write_summary_md(path: Path, rows: list[dict[str, str]]) -> None:
    sorted_rows = sorted(
        rows,
        key=lambda row: float(row["effectiveSampleFraction"]),
        reverse=True,
    )
    reference = float(rows[0]["weightedReference"])
    reference_error = float(rows[0]["weightedReferenceError"])
    fixed_grid = float(rows[0]["fixedGridReference"])

    lines = [
        "# p=100 Proposal Comparison",
        "",
        f"Weighted reference: `{reference:.6g} +/- {reference_error:.3g}`.",
        f"Fixed-grid lower cross-check: `{fixed_grid:.6g}`.",
        "",
        "| variant | mean weight | weight RMS | ESS/N | pull vs weighted ref |",
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
        f"- Highest ESS/N in this comparison: `{best['variant']}` with "
        f"`ESS/N = {float(best['effectiveSampleFraction']):.4g}`."
    )
    smooth_rows = [row for row in rows if row["proposal"] == "smooth"]
    if smooth_rows:
        smooth = smooth_rows[0]
        lines.append(
            f"- Smooth triangular kernel result: `ESS/N = "
            f"{float(smooth['effectiveSampleFraction']):.4g}`, "
            f"`pull = {float(smooth['referencePull']):.3g}`."
        )
    lines.append(
        "- This is a diagnostic comparison only; exact Low weights are still "
        "evaluated at every sampled point."
    )
    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    args = parse_args()
    if args.samples <= 0:
        print("--samples must be positive", file=sys.stderr)
        return 1
    if args.reference_error < 0.0:
        print("--reference-error must be non-negative", file=sys.stderr)
        return 1

    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.summary_csv.parent.mkdir(parents=True, exist_ok=True)
    args.summary_md.parent.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, str]] = []
    for index, variant in enumerate(default_variants()):
        raw = run_variant(args, variant, args.seed + index)
        row = summarize_row(
            raw,
            variant,
            args.reference,
            args.reference_error,
            args.fixed_grid_reference,
        )
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
        "kT",
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
        "fixedGridReference",
        "fixedGridDifference",
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
