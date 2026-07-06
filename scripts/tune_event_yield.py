#!/usr/bin/env python3
"""Tune adaptive event-yield integration settings."""

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
class Config:
    label: str
    refinement_levels: int
    refine_fraction: float
    eta_sub_bins: int
    phi_sub_bins: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Scan adaptive event-yield settings and provisional gates."
    )
    parser.add_argument(
        "summary_csv",
        nargs="?",
        type=Path,
        default=Path("build/event_yield_tuning.csv"),
        help="summary CSV to write",
    )
    parser.add_argument(
        "--summary-md",
        type=Path,
        default=Path("build/event_yield_tuning.md"),
        help="Markdown summary to write",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("build/event_yield_tuning"),
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
    parser.add_argument("--weighted-samples", type=int, default=50000)
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--kt-min", type=float, default=0.002)
    parser.add_argument("--kt-max", type=float, default=0.2)
    parser.add_argument("--eta-min", type=float, default=-5.0)
    parser.add_argument("--eta-max", type=float, default=5.0)
    parser.add_argument("--phi-min", type=float, default=0.0)
    parser.add_argument("--phi-max", type=float, default=2.0 * math.pi)
    parser.add_argument(
        "--particles",
        type=Path,
        default=Path("fixtures/event4_particles.lowmc"),
        help="text fixture to include",
    )
    return parser.parse_args()


def cases(args: argparse.Namespace) -> list[Case]:
    return [
        Case("toy_p1", momentum=1.0),
        Case("toy_p10", momentum=10.0),
        Case("toy_p100", momentum=100.0),
        Case("toy_p1000", momentum=1000.0),
        Case("text_event4", particles=args.particles),
    ]


def configs() -> list[Config]:
    return [
        Config("l5_r0.1_sub2", 5, 0.1, 2, 2),
        Config("l7_r0.1_sub2", 7, 0.1, 2, 2),
        Config("l8_r0.1_sub2", 8, 0.1, 2, 2),
        Config("l9_r0.1_sub2", 9, 0.1, 2, 2),
        Config("l7_r0.05_sub2", 7, 0.05, 2, 2),
        Config("l7_r0.2_sub2", 7, 0.2, 2, 2),
        Config("l7_r0.1_sub3", 7, 0.1, 3, 3),
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


def phase_space_args(args: argparse.Namespace) -> list[str]:
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


def run_yield(args: argparse.Namespace, case: Case, config: Config) -> dict[str, str]:
    output = args.output_dir / f"{case.label}_{config.label}.csv"
    command = [
        str(args.yield_exe),
        str(output),
        *case_args(case),
        *phase_space_args(args),
        "--method",
        "adaptive",
        "--log-kt-bins",
        "16",
        "--adaptive-grid-eta-bins",
        "32",
        "--adaptive-grid-phi-bins",
        "32",
        "--adaptive-refinement-levels",
        str(config.refinement_levels),
        "--adaptive-refine-fraction",
        f"{config.refine_fraction:.17g}",
        "--adaptive-eta-sub-bins",
        str(config.eta_sub_bins),
        "--adaptive-phi-sub-bins",
        str(config.phi_sub_bins),
        "--seed",
        str(args.seed),
    ]
    run_command(command)
    row = read_one_row(output)
    row["yieldCsv"] = str(output)
    return row


def run_weighted(args: argparse.Namespace, case: Case, case_index: int) -> dict[str, str]:
    output = args.output_dir / f"{case.label}_weighted.csv"
    command = [
        str(args.weighted_exe),
        str(output),
        *case_args(case),
        *phase_space_args(args),
        "--samples",
        str(args.weighted_samples),
        "--seed",
        str(args.seed + 5000 + case_index),
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
    row = read_one_row(output)
    row["weightedCsv"] = str(output)
    return row


def weighted_standard_error(weighted: dict[str, str]) -> float:
    samples = float(weighted["samples"])
    variance = float(weighted["varianceWeight"])
    return math.sqrt(max(variance, 0.0) / samples) if samples > 0.0 else math.nan


def summarize_rows(
    case: Case,
    yield_rows: list[tuple[Config, dict[str, str]]],
    weighted: dict[str, str],
) -> list[dict[str, str]]:
    reference = next(
        float(row["lambda_event"])
        for config, row in yield_rows
        if config.label == "l9_r0.1_sub2"
    )
    weighted_mean = float(weighted["meanWeight"])
    weighted_se = weighted_standard_error(weighted)
    weighted_ess_fraction = float(weighted["effectiveSampleFraction"])
    output: list[dict[str, str]] = []
    for config, row in yield_rows:
        lambda_event = float(row["lambda_event"])
        deterministic_difference = lambda_event - reference
        deterministic_relative = (
            abs(deterministic_difference) / abs(reference)
            if reference != 0.0
            else math.nan
        )
        weighted_difference = lambda_event - weighted_mean
        weighted_pull = weighted_difference / weighted_se if weighted_se > 0.0 else math.nan
        input_type = "particles" if case.particles is not None else "toy"
        output.append(
            {
                "case": case.label,
                "inputType": input_type,
                "momentum": "nan" if case.momentum is None else f"{case.momentum:.17g}",
                "particles": "" if case.particles is None else str(case.particles),
                "config": config.label,
                "refinementLevels": str(config.refinement_levels),
                "refineFraction": f"{config.refine_fraction:.17g}",
                "etaSubBins": str(config.eta_sub_bins),
                "phiSubBins": str(config.phi_sub_bins),
                "lambdaEvent": row["lambda_event"],
                "referenceLambda": f"{reference:.17g}",
                "deterministicRelativeDifference": f"{deterministic_relative:.17g}",
                "weightedMean": weighted["meanWeight"],
                "weightedStandardError": f"{weighted_se:.17g}",
                "weightedPull": f"{weighted_pull:.17g}",
                "weightedEffectiveSampleFraction": f"{weighted_ess_fraction:.17g}",
                "evaluations": row["evaluations"],
                "gridCells": row["grid_cells"],
                "yieldCsv": row["yieldCsv"],
                "weightedCsv": weighted["weightedCsv"],
            }
        )
    return output


def gate_label(row: dict[str, str]) -> str:
    rel = float(row["deterministicRelativeDifference"])
    ess = float(row["weightedEffectiveSampleFraction"])
    case = row["case"]
    if case == "toy_p1000":
        if rel < 0.02:
            return "stress-ok"
        return "stress-watch"
    if rel < 0.01 and ess >= 0.05:
        return "ordinary-ok"
    return "watch"


def write_markdown(path: Path, rows: list[dict[str, str]]) -> None:
    selected = [
        row
        for row in rows
        if row["config"] in {"l7_r0.1_sub2", "l8_r0.1_sub2", "l9_r0.1_sub2"}
    ]
    lines = [
        "# Event Yield Tuning",
        "",
        "| case | config | lambda_event | rel vs L9 | weighted mean | pull | ESS/N | gate |",
        "| --- | --- | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    for row in selected:
        lines.append(
            "| "
            f"{row['case']} | "
            f"{row['config']} | "
            f"{float(row['lambdaEvent']):.6g} | "
            f"{float(row['deterministicRelativeDifference']):.3g} | "
            f"{float(row['weightedMean']):.6g} | "
            f"{float(row['weightedPull']):.3g} | "
            f"{float(row['weightedEffectiveSampleFraction']):.3g} | "
            f"{gate_label(row)} |"
        )

    lines.extend(["", "## Provisional Gate", ""])
    lines.append(
        "- Ordinary current diagnostic gate: adaptive L7 or better, deterministic "
        "relative difference below `1%` versus the L9 tuning reference, and "
        "weighted-reference `ESS/N >= 0.05` when using the current weighted "
        "comparison."
    )
    lines.append(
        "- Extreme high-boost stress gate: treat `p = 1000 GeV` separately for now; "
        "L9 is the current reference row, while L7/L8 are convergence probes. "
        "The current weighted reference remains secondary because its ESS/N is low."
    )
    lines.append(
        "- The tuning reference is not a publication-quality proof; it is a "
        "diagnostic gate for deciding when to connect `lambda_event` to the "
        "next count-plus-shape prototype."
    )
    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    args = parse_args()
    if args.weighted_samples <= 0:
        print("--weighted-samples must be positive", file=sys.stderr)
        return 1
    if args.kt_min <= 0.0 or args.kt_max <= args.kt_min:
        print("kT bounds must satisfy 0 < kt-min < kt-max", file=sys.stderr)
        return 1

    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.summary_csv.parent.mkdir(parents=True, exist_ok=True)
    args.summary_md.parent.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, str]] = []
    for case_index, case in enumerate(cases(args)):
        weighted = run_weighted(args, case, case_index)
        case_yields = [(config, run_yield(args, case, config)) for config in configs()]
        case_rows = summarize_rows(case, case_yields, weighted)
        rows.extend(case_rows)
        l7 = next(row for row in case_rows if row["config"] == "l7_r0.1_sub2")
        l9 = next(row for row in case_rows if row["config"] == "l9_r0.1_sub2")
        print(
            f"{case.label}: L7={float(l7['lambdaEvent']):.6g} "
            f"L9={float(l9['lambdaEvent']):.6g} "
            f"rel={float(l7['deterministicRelativeDifference']):.3g} "
            f"weighted={float(l9['weightedMean']):.6g} "
            f"ESS/N={float(l9['weightedEffectiveSampleFraction']):.3g}"
        )

    fieldnames = [
        "case",
        "inputType",
        "momentum",
        "particles",
        "config",
        "refinementLevels",
        "refineFraction",
        "etaSubBins",
        "phiSubBins",
        "lambdaEvent",
        "referenceLambda",
        "deterministicRelativeDifference",
        "weightedMean",
        "weightedStandardError",
        "weightedPull",
        "weightedEffectiveSampleFraction",
        "evaluations",
        "gridCells",
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
