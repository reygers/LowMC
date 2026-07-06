#!/usr/bin/env python3
"""Focused p=100 GeV grid-reference and proposal-shape checks."""

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
class Variant:
    label: str
    proposal: str
    uniform_fraction: float = 0.2
    width_scale: float = 5.0
    min_half_width: float = 0.02
    max_half_width: float = 0.8
    grid_eta_bins: int = 16
    grid_phi_bins: int = 16
    grid_refinement_levels: int = 2
    grid_refine_fraction: float = 0.1
    grid_uniform_fraction: float = 0.2


DEFAULT_GRIDS = [
    GridSize(401, 481),
    GridSize(801, 961),
    GridSize(1201, 1441),
    GridSize(1601, 1921),
]

DEFAULT_VARIANTS = [
    Variant(label="uniform", proposal="uniform"),
    Variant(label="adaptive_w10_min0.02", proposal="adaptive", width_scale=10.0, min_half_width=0.02),
    Variant(label="adaptive_w20_min0.02", proposal="adaptive", width_scale=20.0, min_half_width=0.02),
    Variant(label="adaptive_w40_min0.02", proposal="adaptive", width_scale=40.0, min_half_width=0.02),
    Variant(label="adaptive_w80_min0.02", proposal="adaptive", width_scale=80.0, min_half_width=0.02),
    Variant(label="adaptive_w20_min0.05", proposal="adaptive", width_scale=20.0, min_half_width=0.05),
    Variant(label="adaptive_w40_min0.05", proposal="adaptive", width_scale=40.0, min_half_width=0.05),
    Variant(label="adaptive_w80_min0.05", proposal="adaptive", width_scale=80.0, min_half_width=0.05),
    Variant(label="adaptive_w40_min0.1", proposal="adaptive", width_scale=40.0, min_half_width=0.1),
    Variant(label="adaptive_w80_min0.1", proposal="adaptive", width_scale=80.0, min_half_width=0.1),
    Variant(label="grid_16x16_l2_f0.1", proposal="grid"),
    Variant(label="grid_24x24_l2_f0.1", proposal="grid", grid_eta_bins=24, grid_phi_bins=24),
    Variant(label="grid_32x32_l2_f0.1", proposal="grid", grid_eta_bins=32, grid_phi_bins=32),
    Variant(label="grid_24x24_l3_f0.1", proposal="grid", grid_eta_bins=24, grid_phi_bins=24, grid_refinement_levels=3),
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
            "Run p=100 GeV midrapidity grid-convergence and local-box proposal checks."
        )
    )
    parser.add_argument(
        "summary_csv",
        nargs="?",
        type=Path,
        default=Path("build/p100_focus_summary.csv"),
        help="summary CSV to write",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("build/p100_focus"),
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
    parser.add_argument("--samples", type=int, default=100000, help="samples per proposal variant")
    parser.add_argument("--seed", type=int, default=12345, help="base random seed")
    parser.add_argument("--momentum", type=float, default=100.0, help="toy pion momentum in GeV")
    parser.add_argument("--kt", type=float, default=0.02, help="fixed kT in GeV")
    parser.add_argument(
        "--grid",
        action="append",
        type=parse_grid,
        help="grid size ETABINSxPHIBINS; may be repeated",
    )
    parser.add_argument(
        "--summary-md",
        type=Path,
        default=Path("build/p100_focus_summary.md"),
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


def compute_grid_references(args: argparse.Namespace, grids: list[GridSize]) -> list[dict[str, float | str]]:
    rows: list[dict[str, float | str]] = []
    grid_dir = args.output_dir / "grids"
    grid_dir.mkdir(parents=True, exist_ok=True)

    previous_integral = math.nan
    for grid in grids:
        grid_csv = grid_dir / f"p{args.momentum:g}_{grid.eta_bins}x{grid.phi_bins}.csv"
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
        delta = integral - previous_integral if math.isfinite(previous_integral) else math.nan
        relative_delta = delta / integral if math.isfinite(delta) and integral != 0.0 else math.nan
        rows.append(
            {
                "kind": "grid",
                "variant": f"grid_{grid.eta_bins}x{grid.phi_bins}",
                "etaBins": grid.eta_bins,
                "phiBins": grid.phi_bins,
                "gridIntegral": integral,
                "gridDeltaFromPrevious": delta,
                "gridRelativeDeltaFromPrevious": relative_delta,
                "csv": str(grid_csv),
            }
        )
        previous_integral = integral

    return rows


def run_variant(args: argparse.Namespace, variant: Variant, seed: int) -> dict[str, str]:
    variant_csv = args.output_dir / f"{variant.label}.csv"
    command = [
        str(args.study_exe),
        str(variant_csv),
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
        "--adaptive-uniform-fraction",
        f"{variant.uniform_fraction:.17g}",
        "--adaptive-width-scale",
        f"{variant.width_scale:.17g}",
        "--adaptive-min-half-width",
        f"{variant.min_half_width:.17g}",
        "--adaptive-max-half-width",
        f"{variant.max_half_width:.17g}",
    ]
    if variant.proposal == "grid":
        command.extend(
            [
                "--grid-eta-bins",
                str(variant.grid_eta_bins),
                "--grid-phi-bins",
                str(variant.grid_phi_bins),
                "--grid-refinement-levels",
                str(variant.grid_refinement_levels),
                "--grid-refine-fraction",
                f"{variant.grid_refine_fraction:.17g}",
                "--grid-uniform-fraction",
                f"{variant.grid_uniform_fraction:.17g}",
            ]
        )
    run_command(command)

    with variant_csv.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 1:
        raise RuntimeError(f"expected one row in {variant_csv}, got {len(rows)}")
    rows[0]["csv"] = str(variant_csv)
    return rows[0]


def summarize_variant(
    row: dict[str, str],
    variant: Variant,
    reference: float,
) -> dict[str, float | str]:
    samples = float(row["samples"])
    mean = float(row["meanWeight"])
    variance = float(row["varianceWeight"])
    rms = math.sqrt(max(variance, 0.0))
    standard_error = math.sqrt(variance / samples) if samples > 0.0 else math.nan
    difference = mean - reference
    pull = difference / standard_error if standard_error > 0.0 else math.nan

    return {
        "kind": "proposal",
        "variant": variant.label,
        "proposal": variant.proposal,
        "adaptiveUniformFraction": variant.uniform_fraction,
        "adaptiveWidthScale": variant.width_scale,
        "adaptiveMinHalfWidth": variant.min_half_width,
        "adaptiveMaxHalfWidth": variant.max_half_width,
        "samples": int(samples),
        "meanWeight": mean,
        "weightRMS": rms,
        "weightedStandardError": standard_error,
        "effectiveSampleFraction": float(row["effectiveSampleFraction"]),
        "referenceGridIntegral": reference,
        "referenceDifference": difference,
        "referencePull": pull,
        "csv": row["csv"],
    }


def proposal_rows(args: argparse.Namespace, reference: float) -> list[dict[str, float | str]]:
    rows: list[dict[str, float | str]] = []
    for index, variant in enumerate(DEFAULT_VARIANTS):
        raw_row = run_variant(args, variant, args.seed + index)
        rows.append(summarize_variant(raw_row, variant, reference))
    return rows


def write_summary_csv(
    path: Path,
    grid_rows: list[dict[str, float | str]],
    proposal_summary_rows: list[dict[str, float | str]],
) -> None:
    fieldnames = [
        "kind",
        "variant",
        "proposal",
        "etaBins",
        "phiBins",
        "gridIntegral",
        "gridDeltaFromPrevious",
        "gridRelativeDeltaFromPrevious",
        "adaptiveUniformFraction",
        "adaptiveWidthScale",
        "adaptiveMinHalfWidth",
        "adaptiveMaxHalfWidth",
        "samples",
        "meanWeight",
        "weightRMS",
        "weightedStandardError",
        "effectiveSampleFraction",
        "referenceGridIntegral",
        "referenceDifference",
        "referencePull",
        "csv",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(grid_rows)
        writer.writerows(proposal_summary_rows)


def write_summary_md(
    path: Path,
    grid_rows: list[dict[str, float | str]],
    proposal_summary_rows: list[dict[str, float | str]],
) -> None:
    best_ess = max(proposal_summary_rows, key=lambda row: float(row["effectiveSampleFraction"]))
    stable_rows = [
        row for row in proposal_summary_rows
        if math.isfinite(float(row["referencePull"])) and abs(float(row["referencePull"])) <= 3.0
    ]
    best_stable = (
        max(stable_rows, key=lambda row: float(row["effectiveSampleFraction"]))
        if stable_rows else None
    )

    lines = [
        "# p = 100 GeV Focus Summary",
        "",
        "## Grid Reference",
        "",
        "| eta bins | phi bins | integral | relative change from previous |",
        "| ---: | ---: | ---: | ---: |",
    ]
    for row in grid_rows:
        lines.append(
            "| "
            f"{int(row['etaBins'])} | {int(row['phiBins'])} | "
            f"{float(row['gridIntegral']):.8g} | "
            f"{float(row['gridRelativeDeltaFromPrevious']):.4g} |"
        )

    lines.extend(
        [
            "",
            "## Proposal Variants",
            "",
            "| variant | mean weight | weight RMS | ESS/N | pull |",
            "| --- | ---: | ---: | ---: | ---: |",
        ]
    )
    for row in sorted(
        proposal_summary_rows,
        key=lambda item: float(item["effectiveSampleFraction"]),
        reverse=True,
    ):
        lines.append(
            "| "
            f"{row['variant']} | "
            f"{float(row['meanWeight']):.6g} | "
            f"{float(row['weightRMS']):.6g} | "
            f"{float(row['effectiveSampleFraction']):.6g} | "
            f"{float(row['referencePull']):.3g} |"
        )

    lines.extend(["", "## Interpretation", ""])
    lines.append(
        f"- Highest ESS/N tested: `{best_ess['variant']}` with "
        f"`ESS/N = {float(best_ess['effectiveSampleFraction']):.4g}`."
    )
    if best_stable is not None:
        lines.append(
            f"- Highest ESS/N with `|pull| <= 3`: `{best_stable['variant']}` with "
            f"`ESS/N = {float(best_stable['effectiveSampleFraction']):.4g}`."
        )
    else:
        lines.append("- No tested proposal variant satisfies `|pull| <= 3`.")
    lines.append(
        "- `p = 100 GeV` is a useful high-energy tuning benchmark; "
        "`p = 1000 GeV` should remain an extreme stress case."
    )
    lines.append(
        "- Grid-adaptive variants with higher refinement can strongly improve ESS/N, "
        "but large pulls should be treated as a warning to tune the map and/or improve "
        "the grid reference before declaring a default."
    )

    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    args = parse_args()
    if args.samples <= 0:
        print("--samples must be positive", file=sys.stderr)
        return 1
    if args.momentum <= 0.0:
        print("--momentum must be positive", file=sys.stderr)
        return 1

    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.summary_csv.parent.mkdir(parents=True, exist_ok=True)
    args.summary_md.parent.mkdir(parents=True, exist_ok=True)

    grids = args.grid if args.grid else DEFAULT_GRIDS
    grid_rows = compute_grid_references(args, grids)
    reference = float(grid_rows[-1]["gridIntegral"])
    proposals = proposal_rows(args, reference)

    write_summary_csv(args.summary_csv, grid_rows, proposals)
    write_summary_md(args.summary_md, grid_rows, proposals)

    print(f"reference={reference:.8g} from {grid_rows[-1]['variant']}")
    best = max(proposals, key=lambda row: float(row["effectiveSampleFraction"]))
    print(
        f"best_ess={best['variant']} ESS/N={float(best['effectiveSampleFraction']):.4g} "
        f"pull={float(best['referencePull']):.3g}"
    )

    stable = [
        row for row in proposals
        if math.isfinite(float(row["referencePull"])) and abs(float(row["referencePull"])) <= 3.0
    ]
    if stable:
        best_stable = max(stable, key=lambda row: float(row["effectiveSampleFraction"]))
        print(
            f"best_stable={best_stable['variant']} "
            f"ESS/N={float(best_stable['effectiveSampleFraction']):.4g} "
            f"pull={float(best_stable['referencePull']):.3g}"
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
