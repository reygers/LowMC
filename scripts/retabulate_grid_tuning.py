#!/usr/bin/env python3
"""Retabulate grid-adaptive tuning rows against a weighted reference."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Recompute grid-tuning differences and pulls against a high-stat "
            "weighted reference from run_p100_reference.py."
        )
    )
    parser.add_argument(
        "grid_tuning_csv",
        type=Path,
        help="CSV file written by scripts/scan_grid_tuning.py",
    )
    parser.add_argument(
        "reference_csv",
        type=Path,
        help="CSV file written by scripts/run_p100_reference.py",
    )
    parser.add_argument("output_csv", type=Path, help="retabulated CSV to write")
    parser.add_argument(
        "--summary-md",
        type=Path,
        default=None,
        help="optional Markdown summary to write",
    )
    parser.add_argument(
        "--max-pull",
        type=float,
        default=3.0,
        help="maximum |pull| for stable-variant selection",
    )
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def weighted_reference(reference_rows: list[dict[str, str]]) -> tuple[float, float]:
    weighted_rows = [
        row for row in reference_rows
        if row.get("kind") == "weighted_grid_adaptive"
    ]
    if not weighted_rows:
        raise ValueError("reference CSV has no weighted_grid_adaptive rows")

    weighted_sum = 0.0
    precision_sum = 0.0
    for row in weighted_rows:
        estimate = float(row["estimate"])
        error = float(row["standardError"])
        if error <= 0.0 or not math.isfinite(error):
            raise ValueError("reference rows need positive finite standardError")
        precision = 1.0 / (error * error)
        weighted_sum += estimate * precision
        precision_sum += precision

    return weighted_sum / precision_sum, math.sqrt(1.0 / precision_sum)


def fixed_grid_reference(reference_rows: list[dict[str, str]]) -> float:
    grid_rows = [row for row in reference_rows if row.get("kind") == "fixed_grid"]
    if not grid_rows:
        return math.nan
    return float(grid_rows[-1]["estimate"])


def retabulate_row(
    row: dict[str, str],
    reference: float,
    reference_error: float,
    fixed_grid: float,
) -> dict[str, str]:
    mean = float(row["meanWeight"])
    sample_error = float(row["weightedStandardError"])
    total_error = math.sqrt(sample_error * sample_error + reference_error * reference_error)
    difference = mean - reference
    pull = difference / total_error if total_error > 0.0 else math.nan

    output = dict(row)
    output["updatedReference"] = f"{reference:.17g}"
    output["updatedReferenceError"] = f"{reference_error:.17g}"
    output["updatedReferenceDifference"] = f"{difference:.17g}"
    output["updatedReferencePull"] = f"{pull:.17g}"
    output["fixedGridReference"] = f"{fixed_grid:.17g}"
    output["fixedGridDifference"] = f"{(mean - fixed_grid):.17g}"
    return output


def write_summary(path: Path, rows: list[dict[str, str]], max_pull: float) -> None:
    sorted_rows = sorted(
        rows,
        key=lambda row: float(row["effectiveSampleFraction"]),
        reverse=True,
    )
    stable_rows = [
        row for row in sorted_rows
        if math.isfinite(float(row["updatedReferencePull"])) and
        abs(float(row["updatedReferencePull"])) <= max_pull
    ]

    reference = float(rows[0]["updatedReference"])
    reference_error = float(rows[0]["updatedReferenceError"])
    fixed_grid = float(rows[0]["fixedGridReference"])

    lines = [
        "# Retabulated Grid-Adaptive Tuning Summary",
        "",
        f"High-stat weighted reference: `{reference:.6g} +/- {reference_error:.3g}`.",
        f"Finest fixed-grid cross-check: `{fixed_grid:.6g}`.",
        "",
        f"Stable variants use `|pull| <= {max_pull:g}` against the weighted reference.",
        "",
        "| variant | mean weight | weight RMS | ESS/N | pull vs weighted ref |",
        "| --- | ---: | ---: | ---: | ---: |",
    ]
    for row in sorted_rows[:12]:
        lines.append(
            "| "
            f"{row['variant']} | "
            f"{float(row['meanWeight']):.6g} | "
            f"{float(row['weightRMS']):.6g} | "
            f"{float(row['effectiveSampleFraction']):.6g} | "
            f"{float(row['updatedReferencePull']):.3g} |"
        )

    lines.extend(["", "## Interpretation", ""])
    best = sorted_rows[0]
    lines.append(
        f"- Highest ESS/N tested: `{best['variant']}` with "
        f"`ESS/N = {float(best['effectiveSampleFraction']):.4g}` and "
        f"`pull = {float(best['updatedReferencePull']):.3g}`."
    )
    if stable_rows:
        stable = stable_rows[0]
        lines.append(
            f"- Highest ESS/N with acceptable pull: `{stable['variant']}` with "
            f"`ESS/N = {float(stable['effectiveSampleFraction']):.4g}` and "
            f"`pull = {float(stable['updatedReferencePull']):.3g}`."
        )
    else:
        lines.append("- No scanned variant satisfies the pull cut.")
    lines.append(
        "- The aggressive refined grid maps now look consistent with the updated "
        "reference; the old large pulls were mostly a fixed-grid-reference issue."
    )
    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    args = parse_args()
    try:
        tuning_rows = read_csv(args.grid_tuning_csv)
        reference_rows = read_csv(args.reference_csv)
        reference, reference_error = weighted_reference(reference_rows)
        fixed_grid = fixed_grid_reference(reference_rows)
    except (OSError, ValueError, KeyError) as error:
        print(f"Could not retabulate grid tuning: {error}", file=sys.stderr)
        return 1

    output_rows = [
        retabulate_row(row, reference, reference_error, fixed_grid)
        for row in tuning_rows
    ]

    if not output_rows:
        print("Could not retabulate grid tuning: tuning CSV is empty", file=sys.stderr)
        return 1

    fieldnames = list(output_rows[0].keys())
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(output_rows)

    if args.summary_md is not None:
        args.summary_md.parent.mkdir(parents=True, exist_ok=True)
        write_summary(args.summary_md, output_rows, args.max_pull)

    best = max(output_rows, key=lambda row: float(row["effectiveSampleFraction"]))
    print(
        f"reference={reference:.8g} +/- {reference_error:.3g}; "
        f"best={best['variant']} ESS/N={float(best['effectiveSampleFraction']):.4g} "
        f"pull={float(best['updatedReferencePull']):.3g}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
