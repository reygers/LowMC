#!/usr/bin/env python3
"""Plot and summarize adaptive eta-phi tuning scan output."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from collections import defaultdict
from pathlib import Path


REQUIRED_COLUMNS = {
    "variant",
    "proposal",
    "momentum",
    "meanWeight",
    "weightRMS",
    "effectiveSampleFraction",
    "referencePull",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create a PNG and Markdown summary from an adaptive tuning scan CSV."
    )
    parser.add_argument(
        "input_csv",
        type=Path,
        help="CSV file written by scripts/scan_adaptive_tuning.py",
    )
    parser.add_argument("output_png", type=Path, help="PNG file to write")
    parser.add_argument(
        "--summary-md",
        type=Path,
        default=None,
        help="optional Markdown summary file to write",
    )
    parser.add_argument(
        "--max-pull",
        type=float,
        default=3.0,
        help="maximum |pull| for the preferred-variant table",
    )
    return parser.parse_args()


def parse_float(row: dict[str, str], name: str) -> float:
    value = float(row[name])
    if not math.isfinite(value):
        return math.nan
    return value


def read_rows(path: Path) -> list[dict[str, float | str]]:
    rows: list[dict[str, float | str]] = []

    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = set(reader.fieldnames or [])
        missing = REQUIRED_COLUMNS - fieldnames
        if missing:
            names = ", ".join(sorted(missing))
            raise ValueError(f"input CSV is missing required column(s): {names}")

        for row in reader:
            parsed = {
                "variant": row["variant"],
                "proposal": row["proposal"],
                "momentum": parse_float(row, "momentum"),
                "meanWeight": parse_float(row, "meanWeight"),
                "weightRMS": parse_float(row, "weightRMS"),
                "effectiveSampleFraction": parse_float(row, "effectiveSampleFraction"),
                "referencePull": parse_float(row, "referencePull"),
            }
            rows.append(parsed)

    if not rows:
        raise ValueError("input CSV is empty")

    return rows


def is_usable_pull(row: dict[str, float | str], max_pull: float) -> bool:
    pull = row["referencePull"]
    return isinstance(pull, float) and math.isfinite(pull) and abs(pull) <= max_pull


def best_by_momentum(
    rows: list[dict[str, float | str]],
    max_pull: float,
) -> dict[float, dict[str, float | str]]:
    best: dict[float, dict[str, float | str]] = {}
    fallback: dict[float, dict[str, float | str]] = {}

    for row in rows:
        momentum = row["momentum"]
        if not isinstance(momentum, float):
            continue
        ess = row["effectiveSampleFraction"]
        if not isinstance(ess, float):
            continue

        if momentum not in fallback or ess > float(fallback[momentum]["effectiveSampleFraction"]):
            fallback[momentum] = row
        if is_usable_pull(row, max_pull):
            if momentum not in best or ess > float(best[momentum]["effectiveSampleFraction"]):
                best[momentum] = row

    for momentum, row in fallback.items():
        best.setdefault(momentum, row)

    return best


def grouped_by_variant(rows: list[dict[str, float | str]]) -> dict[str, list[dict[str, float | str]]]:
    grouped: dict[str, list[dict[str, float | str]]] = defaultdict(list)
    for row in rows:
        grouped[str(row["variant"])].append(row)
    for variant_rows in grouped.values():
        variant_rows.sort(key=lambda row: float(row["momentum"]))
    return dict(sorted(grouped.items()))


def write_summary(
    path: Path,
    rows: list[dict[str, float | str]],
    max_pull: float,
) -> None:
    best = best_by_momentum(rows, max_pull)

    lines = [
        "# Adaptive Tuning Summary",
        "",
        "Preferred variants choose the highest `ESS/N` among rows with acceptable",
        f"`|pull| <= {max_pull:g}` when available; otherwise they fall back to the",
        "highest `ESS/N` row for that momentum.",
        "",
        "| momentum [GeV] | variant | mean weight | weight RMS | ESS/N | pull |",
        "| ---: | --- | ---: | ---: | ---: | ---: |",
    ]

    for momentum in sorted(best):
        row = best[momentum]
        lines.append(
            "| "
            f"{momentum:g} | {row['variant']} | "
            f"{float(row['meanWeight']):.6g} | "
            f"{float(row['weightRMS']):.6g} | "
            f"{float(row['effectiveSampleFraction']):.6g} | "
            f"{float(row['referencePull']):.3g} |"
        )

    unstable = [
        row for row in rows
        if str(row["proposal"]) == "adaptive" and not is_usable_pull(row, max_pull)
    ]
    if unstable:
        lines.extend(
            [
                "",
                "Rows with large reference pulls remain unstable or reference-limited.",
                "They should not be used to declare a production default.",
            ]
        )

    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    args = parse_args()

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print(
            "matplotlib is required for plotting. Install it or inspect the CSV directly.",
            file=sys.stderr,
        )
        return 2

    try:
        rows = read_rows(args.input_csv)
    except (OSError, ValueError) as error:
        print(f"Could not read adaptive tuning CSV: {error}", file=sys.stderr)
        return 1

    grouped = grouped_by_variant(rows)
    fig, axes = plt.subplots(3, 1, figsize=(9, 9), sharex=True, constrained_layout=True)

    for variant, variant_rows in grouped.items():
        momenta = [float(row["momentum"]) for row in variant_rows]
        rms = [float(row["weightRMS"]) for row in variant_rows]
        ess = [float(row["effectiveSampleFraction"]) for row in variant_rows]
        pulls = [abs(float(row["referencePull"])) for row in variant_rows]
        linewidth = 2.2 if variant == "uniform" else 1.0
        alpha = 1.0 if variant == "uniform" else 0.65
        axes[0].plot(momenta, rms, marker="o", linewidth=linewidth, alpha=alpha, label=variant)
        axes[1].plot(momenta, ess, marker="o", linewidth=linewidth, alpha=alpha)
        axes[2].plot(momenta, pulls, marker="o", linewidth=linewidth, alpha=alpha)

    axes[0].set_ylabel("weight RMS")
    axes[0].set_yscale("log")
    axes[0].set_title("Adaptive eta-phi tuning scan")
    axes[0].legend(fontsize="x-small", ncols=2)

    axes[1].set_ylabel("ESS / samples")
    axes[1].set_yscale("log")

    axes[2].set_ylabel("|pull| vs grid")
    axes[2].set_xlabel("charged-pion momentum [GeV]")
    axes[2].set_yscale("log")

    for ax in axes:
        ax.set_xscale("log")
        ax.grid(True, which="both", alpha=0.25)

    fig.savefig(args.output_png, dpi=160)
    plt.close(fig)

    if args.summary_md is not None:
        write_summary(args.summary_md, rows, args.max_pull)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
