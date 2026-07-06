#!/usr/bin/env python3
"""Stress-test the shape-unweighting API on selected two-pion toy cases."""

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
    name: str
    momentum: float
    particle_eta: float = 0.0
    particle_phi: float = math.pi / 2.0


DEFAULT_CASES = [
    Case("mid_p100", 100.0),
    Case("mid_p1000", 1000.0),
    Case("forward_eta4_p100", 100.0, 4.0),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Stress-test shape-unweighting API diagnostics on selected two-pion toy cases."
    )
    parser.add_argument(
        "summary_csv",
        nargs="?",
        type=Path,
        default=Path("build/shape_unweighting_stress_summary.csv"),
        help="summary CSV to write",
    )
    parser.add_argument(
        "--diagnostic-exe",
        type=Path,
        default=Path("build/lowmc_shape_unweighting_diagnostic"),
    )
    parser.add_argument(
        "--summary-md",
        type=Path,
        default=Path("build/shape_unweighting_stress_summary.md"),
    )
    parser.add_argument("--calibration-trials", type=int, default=1000)
    parser.add_argument("--trials", type=int, default=2000)
    parser.add_argument("--seed", type=int, action="append", dest="seeds")
    parser.add_argument("--safety-factor", type=float, default=3.0)
    parser.add_argument("--probe-bins", type=int, default=11)
    parser.add_argument("--probe-rule", choices=("interior", "edge"), default="interior")
    parser.add_argument(
        "--failure-policy",
        choices=("diagnostic", "strict"),
        default="diagnostic",
        help="shape-unweighting envelope failure policy",
    )
    return parser.parse_args()


def run_command(command: list[str]) -> str:
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        if result.stdout:
            print(result.stdout, file=sys.stdout, end="")
        if result.stderr:
            print(result.stderr, file=sys.stderr, end="")
        raise subprocess.CalledProcessError(result.returncode, command)
    return result.stdout.strip()


def parse_key_value_line(line: str) -> dict[str, str]:
    row: dict[str, str] = {}
    for item in line.split():
        if "=" not in item:
            continue
        key, value = item.split("=", 1)
        row[key] = value
    return row


def run_case(args: argparse.Namespace, case: Case, seed: int) -> dict[str, str]:
    command = [
        str(args.diagnostic_exe),
        "--calibration-trials",
        str(args.calibration_trials),
        "--trials",
        str(args.trials),
        "--momentum",
        f"{case.momentum:.17g}",
        "--particle-eta",
        f"{case.particle_eta:.17g}",
        "--particle-phi",
        f"{case.particle_phi:.17g}",
        "--seed",
        str(seed),
        "--safety-factor",
        f"{args.safety_factor:.17g}",
        "--probe-bins",
        str(args.probe_bins),
        "--probe-rule",
        args.probe_rule,
        "--failure-policy",
        args.failure_policy,
    ]
    row = parse_key_value_line(run_command(command))
    row["case"] = case.name
    row["seed"] = str(seed)
    return row


def mean(values: list[float]) -> float:
    return sum(values) / len(values) if values else math.nan


def aggregate(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    groups: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        groups.setdefault(row["case"], []).append(row)

    summaries: list[dict[str, str]] = []
    for case, group_rows in groups.items():
        failures = sum(int(row["failed_envelope_count"]) for row in group_rows)
        failed_runs = sum(1 for row in group_rows if int(row["failed_envelope_count"]) > 0)
        non_ok_runs = sum(1 for row in group_rows if row["status"] != "ok")
        summaries.append(
            {
                "case": case,
                "runs": str(len(group_rows)),
                "non_ok_runs": str(non_ok_runs),
                "failed_runs": str(failed_runs),
                "total_failed_envelope_count": str(failures),
                "total_accepted_candidates": str(
                    sum(int(row["accepted_candidates"]) for row in group_rows)
                ),
                "mean_observed_acceptance_efficiency": f"{mean([float(row['observed_acceptance_efficiency']) for row in group_rows]):.17g}",
                "mean_expected_acceptance_efficiency": f"{mean([float(row['expected_acceptance_efficiency']) for row in group_rows]):.17g}",
                "mean_envelope_integral": f"{mean([float(row['envelope_integral']) for row in group_rows]):.17g}",
                "max_trial_bound_ratio": f"{max(float(row['max_trial_bound_ratio']) for row in group_rows):.17g}",
                "mean_suspect_cell_count": f"{mean([float(row['suspect_cell_count']) for row in group_rows]):.17g}",
                "mean_probe_evaluations": f"{mean([float(row['probe_evaluations']) for row in group_rows]):.17g}",
                "total_fallback_envelope_count": str(
                    sum(int(row["fallback_count"]) for row in group_rows)
                ),
                "momentum": group_rows[0]["momentum"],
                "particle_eta": group_rows[0]["particle_eta"],
            }
        )
    return summaries


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise RuntimeError("no rows to write")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def write_summary_md(path: Path, summaries: list[dict[str, str]], args: argparse.Namespace) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Unweighting Bound Stress Summary",
        "",
        f"`lowmc_shape_unweighting_diagnostic`, `{args.probe_bins} x {args.probe_bins}` `{args.probe_rule}` probes, safety factor `{args.safety_factor:g}`, failure policy `{args.failure_policy}`.",
        "",
        "| case | runs | non-ok runs | failed runs | failed samples | accepted candidates | mean observed eff. | mean expected eff. | mean envelope integral | max trial/bound | suspect cells | probe evals | fallback samples |",
        "| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in summaries:
        lines.append(
            "| "
            f"{row['case']} | "
            f"{row['runs']} | "
            f"{row['non_ok_runs']} | "
            f"{row['failed_runs']} | "
            f"{row['total_failed_envelope_count']} | "
            f"{row['total_accepted_candidates']} | "
            f"{float(row['mean_observed_acceptance_efficiency']):.6g} | "
            f"{float(row['mean_expected_acceptance_efficiency']):.6g} | "
            f"{float(row['mean_envelope_integral']):.6g} | "
            f"{float(row['max_trial_bound_ratio']):.6g} | "
            f"{float(row['mean_suspect_cell_count']):.6g} | "
            f"{float(row['mean_probe_evaluations']):.6g} | "
            f"{row['total_fallback_envelope_count']} |"
        )
    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    args = parse_args()
    seeds = args.seeds if args.seeds else [12345, 23456, 34567, 45678, 56789]

    rows: list[dict[str, str]] = []
    for case in DEFAULT_CASES:
        for seed in seeds:
            rows.append(run_case(args, case, int(seed)))

    summaries = aggregate(rows)
    write_csv(args.summary_csv, rows)
    write_summary_md(args.summary_md, summaries, args)
    print(f"wrote {args.summary_csv}")
    print(f"wrote {args.summary_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
