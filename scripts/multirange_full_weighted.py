#!/usr/bin/env python3
"""Check full grid-adaptive sampling across seeds and kT ranges."""

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
class KTRange:
    label: str
    kt_min: float
    kt_max: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Multi-seed and multi-range checks for full grid-adaptive sampling."
    )
    parser.add_argument(
        "summary_csv",
        nargs="?",
        type=Path,
        default=Path("build/full_weighted_multirange.csv"),
        help="per-seed summary CSV to write",
    )
    parser.add_argument(
        "--aggregate-csv",
        type=Path,
        default=Path("build/full_weighted_multirange_aggregate.csv"),
        help="aggregate summary CSV to write",
    )
    parser.add_argument(
        "--summary-md",
        type=Path,
        default=Path("build/full_weighted_multirange.md"),
        help="Markdown summary to write",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("build/full_weighted_multirange"),
        help="directory for individual sampler-study CSV files",
    )
    parser.add_argument(
        "--study-exe",
        type=Path,
        default=Path("build/lowmc_full_weighted_sampler_study"),
        help="path to lowmc_full_weighted_sampler_study",
    )
    parser.add_argument("--samples", type=int, default=100000, help="samples per seed")
    parser.add_argument("--seeds", type=int, default=5, help="number of seeds per kT range")
    parser.add_argument("--seed", type=int, default=12345, help="first random seed")
    parser.add_argument("--momentum", type=float, default=100.0, help="toy pion momentum in GeV")
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
    return parser.parse_args()


def kt_ranges() -> list[KTRange]:
    return [
        KTRange("kt_0.002_0.2", 0.002, 0.2),
        KTRange("kt_0.001_0.1", 0.001, 0.1),
        KTRange("kt_0.02_0.2", 0.02, 0.2),
    ]


def expected_full_reference(args: argparse.Namespace, kt_range: KTRange) -> float:
    return (
        args.fixed_kt_angular_reference *
        args.fixed_kt_reference_kt *
        args.fixed_kt_reference_kt *
        math.log(kt_range.kt_max / kt_range.kt_min)
    )


def run_command(command: list[str]) -> None:
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        if result.stdout:
            print(result.stdout, file=sys.stdout, end="")
        if result.stderr:
            print(result.stderr, file=sys.stderr, end="")
        raise subprocess.CalledProcessError(result.returncode, command)


def run_point(args: argparse.Namespace, kt_range: KTRange, seed: int) -> dict[str, str]:
    output_csv = args.output_dir / f"full_grid_{kt_range.label}_seed{seed}.csv"
    command = [
        str(args.study_exe),
        str(output_csv),
        "--samples",
        str(args.samples),
        "--seed",
        str(seed),
        "--momentum",
        f"{args.momentum:.17g}",
        "--kt-min",
        f"{kt_range.kt_min:.17g}",
        "--kt-max",
        f"{kt_range.kt_max:.17g}",
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
    kt_range: KTRange,
    seed: int,
    reference: float,
) -> dict[str, str]:
    samples = float(raw["samples"])
    mean = float(raw["meanWeight"])
    variance = float(raw["varianceWeight"])
    rms = math.sqrt(max(variance, 0.0))
    error = math.sqrt(variance / samples) if samples > 0.0 else math.nan
    difference = mean - reference
    pull = difference / error if error > 0.0 else math.nan
    ratio = mean / reference if reference > 0.0 else math.nan

    return {
        "range": kt_range.label,
        "seed": str(seed),
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
        "referenceRatio": f"{ratio:.17g}",
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
    for label in sorted({row["range"] for row in rows}):
        group = [row for row in rows if row["range"] == label]
        means = [float(row["meanWeight"]) for row in group]
        rmss = [float(row["weightRMS"]) for row in group]
        ess = [float(row["effectiveSampleFraction"]) for row in group]
        ratios = [float(row["referenceRatio"]) for row in group]
        pulls = [float(row["referencePull"]) for row in group]
        output.append(
            {
                "range": label,
                "kTMin": group[0]["kTMin"],
                "kTMax": group[0]["kTMax"],
                "gridBuildKT": group[0]["gridBuildKT"],
                "seeds": str(len(group)),
                "expectedFullReference": group[0]["expectedFullReference"],
                "meanWeightMean": f"{statistics.mean(means):.17g}",
                "meanWeightStd": f"{sample_stdev(means):.17g}",
                "weightRMSMean": f"{statistics.mean(rmss):.17g}",
                "weightRMSStd": f"{sample_stdev(rmss):.17g}",
                "effectiveSampleFractionMean": f"{statistics.mean(ess):.17g}",
                "effectiveSampleFractionStd": f"{sample_stdev(ess):.17g}",
                "referenceRatioMean": f"{statistics.mean(ratios):.17g}",
                "referenceRatioStd": f"{sample_stdev(ratios):.17g}",
                "referencePullMean": f"{statistics.mean(pulls):.17g}",
                "referencePullStd": f"{sample_stdev(pulls):.17g}",
            }
        )
    return output


def write_summary_md(path: Path, aggregate: list[dict[str, str]]) -> None:
    lines = [
        "# Full Weighted Multi-Range Summary",
        "",
        "| range | seeds | expected | mean weight | weight RMS | ESS/N | mean/expected | pull |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in aggregate:
        lines.append(
            "| "
            f"{row['range']} | "
            f"{row['seeds']} | "
            f"{float(row['expectedFullReference']):.6g} | "
            f"{float(row['meanWeightMean']):.6g} +/- {float(row['meanWeightStd']):.3g} | "
            f"{float(row['weightRMSMean']):.6g} +/- {float(row['weightRMSStd']):.3g} | "
            f"{float(row['effectiveSampleFractionMean']):.6g} +/- "
            f"{float(row['effectiveSampleFractionStd']):.3g} | "
            f"{float(row['referenceRatioMean']):.6g} +/- "
            f"{float(row['referenceRatioStd']):.3g} | "
            f"{float(row['referencePullMean']):.3g} +/- "
            f"{float(row['referencePullStd']):.3g} |"
        )

    lines.extend(["", "## Interpretation", ""])
    lines.append(
        "- The expected value uses the fixed-kT angular reference scaled by "
        "`kTRef^2 ln(kTMax/kTMin)`."
    )
    lines.append(
        "- Stable `mean/expected` values near one indicate the full path respects "
        "the expected logarithmic kT scaling."
    )
    lines.append(
        "- The grid proposal is built at the geometric-mean kT for each range, "
        "while exact Low weights are evaluated at sampled kT."
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

    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.summary_csv.parent.mkdir(parents=True, exist_ok=True)
    args.aggregate_csv.parent.mkdir(parents=True, exist_ok=True)
    args.summary_md.parent.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, str]] = []
    for kt_range in kt_ranges():
        reference = expected_full_reference(args, kt_range)
        for offset in range(args.seeds):
            seed = args.seed + offset
            raw = run_point(args, kt_range, seed)
            row = summarize_row(raw, kt_range, seed, reference)
            rows.append(row)
            print(
                f"{row['range']} seed={seed}: mean={float(row['meanWeight']):.6g} "
                f"RMS={float(row['weightRMS']):.6g} "
                f"ESS/N={float(row['effectiveSampleFraction']):.4g} "
                f"ratio={float(row['referenceRatio']):.4g} "
                f"pull={float(row['referencePull']):.3g}"
            )

    fieldnames = [
        "range",
        "seed",
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
        "referenceRatio",
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
        "range",
        "kTMin",
        "kTMax",
        "gridBuildKT",
        "seeds",
        "expectedFullReference",
        "meanWeightMean",
        "meanWeightStd",
        "weightRMSMean",
        "weightRMSStd",
        "effectiveSampleFractionMean",
        "effectiveSampleFractionStd",
        "referenceRatioMean",
        "referenceRatioStd",
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
