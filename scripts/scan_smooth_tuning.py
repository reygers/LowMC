#!/usr/bin/env python3
"""Scan smooth local-kernel eta-phi proposal parameters for the p=100 toy case."""

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
    role: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Scan triangular smooth-kernel proposal parameters against the "
            "p=100 high-stat weighted reference."
        )
    )
    parser.add_argument(
        "summary_csv",
        nargs="?",
        type=Path,
        default=Path("build/smooth_tuning_summary.csv"),
        help="summary CSV to write",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("build/smooth_tuning"),
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
        default=Path("build/smooth_tuning_summary.md"),
        help="Markdown summary to write",
    )
    parser.add_argument(
        "--max-pull",
        type=float,
        default=3.0,
        help="maximum |pull| for stable-variant selection",
    )
    return parser.parse_args()


def fixed_variants() -> list[Variant]:
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
            "benchmark",
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
            "baseline",
        ),
    ]


def smooth_variants() -> list[Variant]:
    variants: list[Variant] = []
    for width_scale in [10.0, 20.0, 40.0, 80.0]:
        for min_half_width in [0.02, 0.05, 0.1, 0.2]:
            for uniform_fraction in [0.1, 0.2, 0.4]:
                variants.append(
                    Variant(
                        (
                            f"smooth_triangular_w{width_scale:g}"
                            f"_min{min_half_width:g}_u{uniform_fraction:g}"
                        ),
                        "smooth",
                        (
                            "--smooth-uniform-fraction",
                            f"{uniform_fraction:.17g}",
                            "--smooth-width-scale",
                            f"{width_scale:.17g}",
                            "--smooth-min-half-width",
                            f"{min_half_width:.17g}",
                            "--smooth-max-half-width",
                            "0.8",
                        ),
                        "scan",
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
        "role": variant.role,
        "proposal": variant.proposal,
        "smoothUniformFraction": raw.get("smoothUniformFraction", ""),
        "smoothWidthScale": raw.get("smoothWidthScale", ""),
        "smoothMinHalfWidth": raw.get("smoothMinHalfWidth", ""),
        "smoothMaxHalfWidth": raw.get("smoothMaxHalfWidth", ""),
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


def write_summary_md(path: Path, rows: list[dict[str, str]], max_pull: float) -> None:
    sorted_rows = sorted(
        rows,
        key=lambda row: float(row["effectiveSampleFraction"]),
        reverse=True,
    )
    smooth_rows = [row for row in rows if row["proposal"] == "smooth"]
    smooth_sorted = sorted(
        smooth_rows,
        key=lambda row: float(row["effectiveSampleFraction"]),
        reverse=True,
    )
    stable_smooth_rows = [
        row for row in smooth_sorted
        if math.isfinite(float(row["referencePull"])) and
        abs(float(row["referencePull"])) <= max_pull
    ]
    reference = float(rows[0]["weightedReference"])
    reference_error = float(rows[0]["weightedReferenceError"])
    fixed_grid = float(rows[0]["fixedGridReference"])

    lines = [
        "# Smooth-Kernel Tuning Summary",
        "",
        f"Weighted reference: `{reference:.6g} +/- {reference_error:.3g}`.",
        f"Fixed-grid lower cross-check: `{fixed_grid:.6g}`.",
        f"Stable variants use `|pull| <= {max_pull:g}` against the weighted reference.",
        "",
        "| variant | role | mean weight | weight RMS | ESS/N | pull vs weighted ref |",
        "| --- | --- | ---: | ---: | ---: | ---: |",
    ]
    for row in sorted_rows[:14]:
        lines.append(
            "| "
            f"{row['variant']} | "
            f"{row['role']} | "
            f"{float(row['meanWeight']):.6g} | "
            f"{float(row['weightRMS']):.6g} | "
            f"{float(row['effectiveSampleFraction']):.6g} | "
            f"{float(row['referencePull']):.3g} |"
        )

    lines.extend(["", "## Interpretation", ""])
    best_overall = sorted_rows[0]
    lines.append(
        f"- Highest ESS/N overall: `{best_overall['variant']}` with "
        f"`ESS/N = {float(best_overall['effectiveSampleFraction']):.4g}`."
    )
    if smooth_sorted:
        best_smooth = smooth_sorted[0]
        lines.append(
            f"- Highest ESS/N smooth variant: `{best_smooth['variant']}` with "
            f"`ESS/N = {float(best_smooth['effectiveSampleFraction']):.4g}` and "
            f"`pull = {float(best_smooth['referencePull']):.3g}`."
        )
    if stable_smooth_rows:
        stable = stable_smooth_rows[0]
        lines.append(
            f"- Highest ESS/N stable smooth variant: `{stable['variant']}` with "
            f"`ESS/N = {float(stable['effectiveSampleFraction']):.4g}` and "
            f"`pull = {float(stable['referencePull']):.3g}`."
        )
    else:
        lines.append("- No scanned smooth variant satisfies the pull cut.")

    baseline_rows = [row for row in rows if row["role"] == "baseline"]
    benchmark_rows = [row for row in rows if row["role"] == "benchmark"]
    if baseline_rows and smooth_sorted:
        baseline = baseline_rows[0]
        best_smooth = smooth_sorted[0]
        ratio = (
            float(best_smooth["effectiveSampleFraction"]) /
            float(baseline["effectiveSampleFraction"])
        )
        lines.append(
            f"- Best smooth ESS/N is `{ratio:.3g}` times the local-box baseline."
        )
    if benchmark_rows and smooth_sorted:
        benchmark = benchmark_rows[0]
        best_smooth = smooth_sorted[0]
        ratio = (
            float(best_smooth["effectiveSampleFraction"]) /
            float(benchmark["effectiveSampleFraction"])
        )
        lines.append(
            f"- Best smooth ESS/N is `{ratio:.3g}` times the tuned grid benchmark."
        )

    lines.append(
        "- This remains a diagnostic comparison; exact Low weights are evaluated "
        "at every sampled point."
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
    variants = fixed_variants() + smooth_variants()
    for index, variant in enumerate(variants):
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
        "role",
        "proposal",
        "smoothUniformFraction",
        "smoothWidthScale",
        "smoothMinHalfWidth",
        "smoothMaxHalfWidth",
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
    write_summary_md(args.summary_md, rows, args.max_pull)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
