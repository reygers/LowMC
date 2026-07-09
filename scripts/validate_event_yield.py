#!/usr/bin/env python3
"""Validate deterministic event-yield estimates against weighted sampling."""

from __future__ import annotations

import argparse
import csv
import math
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Case:
    label: str
    momentum: float | None = None
    particles: Path | None = None


@dataclass(frozen=True)
class Grid:
    label: str
    method: str
    log_kt_bins: int
    eta_bins: int
    phi_bins: int
    refinement_levels: int = 0
    refine_fraction: float = 0.1


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Scan event-yield convergence and compare with weighted samples."
    )
    parser.add_argument(
        "summary_csv",
        nargs="?",
        type=Path,
        default=Path("build/event_yield_validation.csv"),
        help="summary CSV to write",
    )
    parser.add_argument(
        "--summary-md",
        type=Path,
        default=Path("build/event_yield_validation.md"),
        help="Markdown summary to write",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("build/event_yield_validation"),
        help="directory for individual diagnostic CSV files",
    )
    parser.add_argument(
        "--yield-exe",
        type=Path,
        default=Path("build/lowmc_event_yield_diagnostic"),
        help="path to lowmc_event_yield_diagnostic",
    )
    parser.add_argument(
        "--weighted-exe",
        type=Path,
        default=Path("build/lowmc_full_weighted_sampler_study"),
        help="path to lowmc_full_weighted_sampler_study",
    )
    parser.add_argument("--samples", type=int, default=20000, help="weighted samples per case")
    parser.add_argument("--seed", type=int, default=12345, help="base random seed")
    parser.add_argument("--kt-min", type=float, default=0.002, help="minimum kT in GeV")
    parser.add_argument("--kt-max", type=float, default=0.2, help="maximum kT in GeV")
    parser.add_argument("--eta-min", type=float, default=-5.0, help="minimum photon eta")
    parser.add_argument("--eta-max", type=float, default=5.0, help="maximum photon eta")
    parser.add_argument("--phi-min", type=float, default=0.0, help="minimum photon phi")
    parser.add_argument(
        "--phi-max",
        type=float,
        default=2.0 * math.pi,
        help="maximum photon phi",
    )
    parser.add_argument(
        "--particles",
        type=Path,
        default=Path("fixtures/event4_particles.lowmc"),
        help="text fixture to include in the validation",
    )
    return parser.parse_args()


def grids() -> list[Grid]:
    return [
        Grid("midpoint_coarse", "midpoint", 8, 20, 24),
        Grid("midpoint_medium", "midpoint", 16, 40, 48),
        Grid("midpoint_fine", "midpoint", 32, 80, 96),
        Grid("adaptive_l3", "adaptive", 16, 2, 2, 3, 0.1),
        Grid("adaptive_l5", "adaptive", 16, 2, 2, 5, 0.1),
        Grid("adaptive_l7", "adaptive", 16, 2, 2, 7, 0.1),
        Grid("factorized_l7", "factorized", 16, 2, 2, 7, 0.1),
    ]


def cases(args: argparse.Namespace) -> list[Case]:
    return [
        Case("toy_p1", momentum=1.0),
        Case("toy_p10", momentum=10.0),
        Case("toy_p100", momentum=100.0),
        Case("toy_p1000", momentum=1000.0),
        Case("text_event4", particles=args.particles),
    ]


def run_command(command: list[str]) -> None:
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        if result.stdout:
            print(result.stdout, file=sys.stdout, end="")
        if result.stderr:
            print(result.stderr, file=sys.stderr, end="")
        raise subprocess.CalledProcessError(result.returncode, command)


def read_one_row(path: Path) -> dict[str, str]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 1:
        raise RuntimeError(f"expected one row in {path}, got {len(rows)}")
    return rows[0]


def common_phase_space_args(args: argparse.Namespace) -> list[str]:
    return [
        "--kt-min",
        f"{args.kt_min:.17g}",
        "--kt-max",
        f"{args.kt_max:.17g}",
        "--eta-min",
        f"{args.eta_min:.17g}",
        "--eta-max",
        f"{args.eta_max:.17g}",
        "--phi-min",
        f"{args.phi_min:.17g}",
        "--phi-max",
        f"{args.phi_max:.17g}",
    ]


def case_args(case: Case) -> list[str]:
    if case.particles is not None:
        return ["--particles", str(case.particles)]
    if case.momentum is None:
        raise ValueError(f"case {case.label} has neither momentum nor particles")
    return ["--momentum", f"{case.momentum:.17g}"]


def run_yield(args: argparse.Namespace, case: Case, grid: Grid) -> dict[str, str]:
    output_csv = args.output_dir / f"{case.label}_{grid.label}_yield.csv"
    command = [
        str(args.yield_exe),
        str(output_csv),
        *case_args(case),
        *common_phase_space_args(args),
        "--method",
        grid.method,
        "--log-kt-bins",
        str(grid.log_kt_bins),
        "--seed",
        str(args.seed),
    ]
    if grid.method in {"adaptive", "factorized"}:
        command.extend(
            [
                "--adaptive-grid-eta-bins",
                "32",
                "--adaptive-grid-phi-bins",
                "32",
                "--adaptive-refinement-levels",
                str(grid.refinement_levels),
                "--adaptive-refine-fraction",
                f"{grid.refine_fraction:.17g}",
                "--adaptive-eta-sub-bins",
                str(grid.eta_bins),
                "--adaptive-phi-sub-bins",
                str(grid.phi_bins),
            ]
        )
    else:
        command.extend(
            [
                "--eta-bins",
                str(grid.eta_bins),
                "--phi-bins",
                str(grid.phi_bins),
            ]
        )
    run_command(command)
    row = read_one_row(output_csv)
    row["csv"] = str(output_csv)
    return row


def run_weighted(args: argparse.Namespace, case: Case, index: int) -> dict[str, str]:
    output_csv = args.output_dir / f"{case.label}_weighted.csv"
    command = [
        str(args.weighted_exe),
        str(output_csv),
        *case_args(case),
        *common_phase_space_args(args),
        "--samples",
        str(args.samples),
        "--seed",
        str(args.seed + 1000 + index),
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
    row = read_one_row(output_csv)
    row["csv"] = str(output_csv)
    return row


def weighted_error(weighted: dict[str, str]) -> float:
    samples = float(weighted["samples"])
    variance = float(weighted["varianceWeight"])
    return math.sqrt(max(variance, 0.0) / samples) if samples > 0.0 else math.nan


def summarize_case(
    case: Case,
    grid: Grid,
    yield_row: dict[str, str],
    weighted: dict[str, str],
    previous_lambda: float | None,
) -> dict[str, str]:
    lambda_event = float(yield_row["lambda_event"])
    weighted_mean = float(weighted["meanWeight"])
    weighted_se = weighted_error(weighted)
    difference = lambda_event - weighted_mean
    pull = difference / weighted_se if weighted_se > 0.0 else math.nan
    step_difference = (
        lambda_event - previous_lambda if previous_lambda is not None else math.nan
    )
    step_relative = (
        abs(step_difference) / abs(lambda_event)
        if previous_lambda is not None and lambda_event != 0.0
        else math.nan
    )
    ratio = lambda_event / weighted_mean if weighted_mean != 0.0 else math.nan
    input_type = "particles" if case.particles is not None else "toy"
    momentum = "nan" if case.momentum is None else f"{case.momentum:.17g}"

    return {
        "case": case.label,
        "inputType": input_type,
        "momentum": momentum,
        "particles": "" if case.particles is None else str(case.particles),
        "method": yield_row["method"],
        "grid": grid.label,
        "logKTBins": yield_row["log_kt_bins"],
        "etaBins": yield_row["eta_bins"],
        "phiBins": yield_row["phi_bins"],
        "evaluations": yield_row["evaluations"],
        "gridCells": yield_row["grid_cells"],
        "lambdaEvent": yield_row["lambda_event"],
        "coarseLambdaEvent": yield_row["coarse_lambda_event"],
        "coarseRelativeDifference": yield_row["relative_difference"],
        "stepRelativeDifference": f"{step_relative:.17g}",
        "weightedSamples": weighted["samples"],
        "weightedMean": weighted["meanWeight"],
        "weightedStandardError": f"{weighted_se:.17g}",
        "weightedRMS": f"{math.sqrt(max(float(weighted['varianceWeight']), 0.0)):.17g}",
        "effectiveSampleFraction": weighted["effectiveSampleFraction"],
        "lambdaOverWeightedMean": f"{ratio:.17g}",
        "differenceVsWeighted": f"{difference:.17g}",
        "pullVsWeighted": f"{pull:.17g}",
        "yieldCsv": yield_row["csv"],
        "weightedCsv": weighted["csv"],
    }


def write_markdown(path: Path, rows: list[dict[str, str]]) -> None:
    lines = [
        "# Event Yield Validation",
        "",
        "| case | grid | method | lambda_event | step diff | weighted mean | MC stderr | pull | ESS/N |",
        "| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in rows:
        lines.append(
            "| "
            f"{row['case']} | "
            f"{row['grid']} | "
            f"{row['method']} | "
            f"{float(row['lambdaEvent']):.6g} | "
            f"{float(row['stepRelativeDifference']):.3g} | "
            f"{float(row['weightedMean']):.6g} | "
            f"{float(row['weightedStandardError']):.3g} | "
            f"{float(row['pullVsWeighted']):.3g} | "
            f"{float(row['effectiveSampleFraction']):.3g} |"
        )

    lines.extend(["", "## Interpretation", ""])
    lines.append(
        "- `lambda_event` is the deterministic midpoint estimate in "
        "`(log k_T, eta, phi)` for `method=log-kt-eta-phi-midpoint`, or the "
        "adaptive-cell deterministic estimate for "
        "`method=adaptive-log-kt-grid-eta-phi-midpoint`. The factorized "
        "adaptive method uses the same angular grid once at the geometric-mean "
        "`k_T` and multiplies by `ln(kTMax / kTMin)`."
    )
    lines.append(
        "- `weighted mean` is the Monte Carlo mean photon weight from the existing "
        "grid-adaptive weighted sampler on the same phase space."
    )
    lines.append(
        "- The pull uses only the weighted-sampler standard error; deterministic "
        "grid error is tracked separately by the step and coarse/fine differences."
    )
    lines.append(
        "- The text-fixture row verifies the same comparison path for direct "
        "charged-particle input before any HepMC3 event-writing step."
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
    if args.eta_max <= args.eta_min:
        print("eta bounds must satisfy eta-min < eta-max", file=sys.stderr)
        return 1
    if args.phi_max <= args.phi_min:
        print("phi bounds must satisfy phi-min < phi-max", file=sys.stderr)
        return 1

    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.summary_csv.parent.mkdir(parents=True, exist_ok=True)
    args.summary_md.parent.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, str]] = []
    for case_index, case in enumerate(cases(args)):
        weighted = run_weighted(args, case, case_index)
        previous_lambda_by_method: dict[str, float] = {}
        for grid in grids():
            yield_row = run_yield(args, case, grid)
            previous_lambda = previous_lambda_by_method.get(grid.method)
            row = summarize_case(case, grid, yield_row, weighted, previous_lambda)
            rows.append(row)
            previous_lambda_by_method[grid.method] = float(row["lambdaEvent"])
            print(
                f"{case.label} {grid.label}: lambda={float(row['lambdaEvent']):.6g} "
                f"weighted={float(row['weightedMean']):.6g} "
                f"pull={float(row['pullVsWeighted']):.3g} "
                f"step={float(row['stepRelativeDifference']):.3g}"
            )

    fieldnames = [
        "case",
        "inputType",
        "momentum",
        "particles",
        "method",
        "grid",
        "logKTBins",
        "etaBins",
        "phiBins",
        "evaluations",
        "gridCells",
        "lambdaEvent",
        "coarseLambdaEvent",
        "coarseRelativeDifference",
        "stepRelativeDifference",
        "weightedSamples",
        "weightedMean",
        "weightedStandardError",
        "weightedRMS",
        "effectiveSampleFraction",
        "lambdaOverWeightedMean",
        "differenceVsWeighted",
        "pullVsWeighted",
        "yieldCsv",
        "weightedCsv",
    ]
    with args.summary_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    write_markdown(args.summary_md, rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
