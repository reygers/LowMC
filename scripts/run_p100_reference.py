#!/usr/bin/env python3
"""Check p=100 GeV reference quality with fixed grids and high-stat weighted runs."""

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
class GridSize:
    eta_bins: int
    phi_bins: int


@dataclass(frozen=True)
class WeightedReferenceVariant:
    label: str
    eta_bins: int
    phi_bins: int
    refinement_levels: int
    refine_fraction: float
    uniform_fraction: float


DEFAULT_GRIDS = [
    GridSize(801, 961),
    GridSize(1201, 1441),
    GridSize(1601, 1921),
]

DEFAULT_WEIGHTED_VARIANTS = [
    WeightedReferenceVariant("stable_grid_16x16_l3_r0.2_u0.1", 16, 16, 3, 0.2, 0.1),
    WeightedReferenceVariant("high_ess_grid_32x32_l3_r0.1_u0.1", 32, 32, 3, 0.1, 0.1),
]


def parse_grid(text: str) -> GridSize:
    try:
        eta_text, phi_text = text.lower().split("x", 1)
        grid = GridSize(int(eta_text), int(phi_text))
    except ValueError as error:
        raise argparse.ArgumentTypeError("grid must have the form ETABINSxPHIBINS") from error
    if grid.eta_bins < 2 or grid.phi_bins < 2:
        raise argparse.ArgumentTypeError("grid bin counts must be >= 2")
    return grid


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compare p=100 fixed-grid references with high-stat weighted "
            "grid-adaptive estimates."
        )
    )
    parser.add_argument(
        "summary_csv",
        nargs="?",
        type=Path,
        default=Path("build/p100_reference_summary.csv"),
        help="summary CSV to write",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("build/p100_reference"),
        help="directory for intermediate CSV files",
    )
    parser.add_argument(
        "--grid-exe",
        type=Path,
        default=Path("build/lowmc_two_particle_diagnostic"),
        help="path to lowmc_two_particle_diagnostic",
    )
    parser.add_argument(
        "--study-exe",
        type=Path,
        default=Path("build/lowmc_weighted_sampler_study"),
        help="path to lowmc_weighted_sampler_study",
    )
    parser.add_argument("--samples", type=int, default=1000000, help="samples per weighted variant")
    parser.add_argument("--seed", type=int, default=12345, help="base random seed")
    parser.add_argument("--momentum", type=float, default=100.0, help="toy pion momentum in GeV")
    parser.add_argument("--kt", type=float, default=0.02, help="fixed kT in GeV")
    parser.add_argument(
        "--grid",
        action="append",
        type=parse_grid,
        help="fixed grid size ETABINSxPHIBINS; may be repeated",
    )
    parser.add_argument(
        "--summary-md",
        type=Path,
        default=Path("build/p100_reference_summary.md"),
        help="Markdown summary to write",
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


def compute_grid_rows(args: argparse.Namespace, grids: list[GridSize]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    grid_dir = args.output_dir / "grids"
    grid_dir.mkdir(parents=True, exist_ok=True)
    previous = math.nan

    for grid in grids:
        grid_csv = grid_dir / f"p{args.momentum:g}_{grid.eta_bins}x{grid.phi_bins}.csv"
        if not grid_csv.exists():
            run_command(
                [
                    str(args.grid_exe),
                    str(grid_csv),
                    "--momentum",
                    f"{args.momentum:.17g}",
                    "--eta-bins",
                    str(grid.eta_bins),
                    "--phi-bins",
                    str(grid.phi_bins),
                ]
            )
        integral = integrate_grid(grid_csv, args.alpha)
        delta = integral - previous if math.isfinite(previous) else math.nan
        rows.append(
            {
                "kind": "fixed_grid",
                "label": f"grid_{grid.eta_bins}x{grid.phi_bins}",
                "etaBins": str(grid.eta_bins),
                "phiBins": str(grid.phi_bins),
                "samples": "",
                "estimate": f"{integral:.17g}",
                "standardError": "",
                "differenceFromFinestGrid": "",
                "pullVsFinestGrid": "",
                "deltaFromPreviousGrid": f"{delta:.17g}",
                "relativeDeltaFromPreviousGrid": (
                    f"{(delta / integral):.17g}" if math.isfinite(delta) and integral != 0.0 else "nan"
                ),
                "csv": str(grid_csv),
            }
        )
        previous = integral
    return rows


def run_weighted_variant(
    args: argparse.Namespace,
    variant: WeightedReferenceVariant,
    seed: int,
    fixed_grid_reference: float,
) -> dict[str, str]:
    output_csv = args.output_dir / f"{variant.label}.csv"
    run_command(
        [
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
    )
    with output_csv.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 1:
        raise RuntimeError(f"expected one row in {output_csv}, got {len(rows)}")

    row = rows[0]
    mean = float(row["meanWeight"])
    samples = float(row["samples"])
    variance = float(row["varianceWeight"])
    error = math.sqrt(variance / samples)
    difference = mean - fixed_grid_reference
    pull = difference / error if error > 0.0 else math.nan
    return {
        "kind": "weighted_grid_adaptive",
        "label": variant.label,
        "etaBins": str(variant.eta_bins),
        "phiBins": str(variant.phi_bins),
        "samples": row["samples"],
        "estimate": f"{mean:.17g}",
        "standardError": f"{error:.17g}",
        "differenceFromFinestGrid": f"{difference:.17g}",
        "pullVsFinestGrid": f"{pull:.17g}",
        "deltaFromPreviousGrid": "",
        "relativeDeltaFromPreviousGrid": "",
        "csv": str(output_csv),
    }


def write_markdown(path: Path, rows: list[dict[str, str]]) -> None:
    grid_rows = [row for row in rows if row["kind"] == "fixed_grid"]
    weighted_rows = [row for row in rows if row["kind"] == "weighted_grid_adaptive"]
    finest = float(grid_rows[-1]["estimate"])

    lines = [
        "# p = 100 GeV Reference Summary",
        "",
        "## Fixed Grids",
        "",
        "| grid | integral | relative change from previous |",
        "| --- | ---: | ---: |",
    ]
    for row in grid_rows:
        lines.append(
            f"| {row['label']} | {float(row['estimate']):.8g} | "
            f"{float(row['relativeDeltaFromPreviousGrid']):.4g} |"
        )

    lines.extend(
        [
            "",
            "## Weighted Grid-Adaptive References",
            "",
            "| variant | estimate | stat. error | difference from finest grid | pull |",
            "| --- | ---: | ---: | ---: | ---: |",
        ]
    )
    for row in weighted_rows:
        lines.append(
            f"| {row['label']} | {float(row['estimate']):.8g} | "
            f"{float(row['standardError']):.4g} | "
            f"{float(row['differenceFromFinestGrid']):.4g} | "
            f"{float(row['pullVsFinestGrid']):.3g} |"
        )

    weighted_estimates = [float(row["estimate"]) for row in weighted_rows]
    min_weighted = min(weighted_estimates)
    lines.extend(["", "## Interpretation", ""])
    if min_weighted > finest:
        lines.append(
            "- The high-stat weighted estimates are above the finest fixed grid, "
            "so the fixed-grid reference is likely still systematically low."
        )
    else:
        lines.append(
            "- The high-stat weighted estimates do not all sit above the finest fixed grid; "
            "more cross-checks are needed before replacing the reference."
        )
    lines.append(
        "- The weighted estimates are cheap enough for proposal tuning when using "
        "the current grid-adaptive maps, but they remain stochastic references."
    )
    lines.append(
        "- For regular tuning, use both the finest fixed grid and at least one "
        "high-stat weighted grid-adaptive estimate."
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

    grids = args.grid if args.grid else DEFAULT_GRIDS
    grid_rows = compute_grid_rows(args, grids)
    fixed_grid_reference = float(grid_rows[-1]["estimate"])

    rows = list(grid_rows)
    for index, variant in enumerate(DEFAULT_WEIGHTED_VARIANTS):
        row = run_weighted_variant(args, variant, args.seed + index, fixed_grid_reference)
        rows.append(row)
        print(
            f"{row['label']}: estimate={float(row['estimate']):.8g} "
            f"+/- {float(row['standardError']):.4g}, "
            f"pull_vs_grid={float(row['pullVsFinestGrid']):.3g}"
        )

    fieldnames = [
        "kind",
        "label",
        "etaBins",
        "phiBins",
        "samples",
        "estimate",
        "standardError",
        "differenceFromFinestGrid",
        "pullVsFinestGrid",
        "deltaFromPreviousGrid",
        "relativeDeltaFromPreviousGrid",
        "csv",
    ]
    with args.summary_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    write_markdown(args.summary_md, rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
