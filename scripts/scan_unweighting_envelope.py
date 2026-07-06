#!/usr/bin/env python3
"""Scan empirical unweighting envelopes for the two-pion toy diagnostic."""

from __future__ import annotations

import argparse
import csv
import math
import subprocess
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Scan lowmc_unweighting_diagnostic over toy momenta, seeds, and safety factors."
    )
    parser.add_argument(
        "summary_csv",
        nargs="?",
        type=Path,
        default=Path("build/unweighting_envelope_scan.csv"),
        help="summary CSV to write",
    )
    parser.add_argument(
        "--diagnostic-exe",
        type=Path,
        default=Path("build/lowmc_unweighting_diagnostic"),
        help="path to lowmc_unweighting_diagnostic",
    )
    parser.add_argument(
        "--summary-md",
        type=Path,
        default=Path("build/unweighting_envelope_scan.md"),
        help="Markdown summary to write",
    )
    parser.add_argument("--calibration-samples", type=int, default=1000)
    parser.add_argument("--trial-samples", type=int, default=1000)
    parser.add_argument("--momentum", type=float, action="append", dest="momenta")
    parser.add_argument("--particle-eta", type=float, default=0.0)
    parser.add_argument("--particle-phi", type=float, default=math.pi / 2.0)
    parser.add_argument("--seed", type=int, action="append", dest="seeds")
    parser.add_argument("--safety-factor", type=float, action="append", dest="safety_factors")
    parser.add_argument(
        "--envelope-mode",
        choices=("global", "cell", "probed-cell"),
        action="append",
        dest="envelope_modes",
        help="envelope mode to scan; may be repeated",
    )
    parser.add_argument("--probe-bins", type=int, default=11)
    parser.add_argument(
        "--probe-rule",
        choices=("interior", "edge"),
        action="append",
        dest="probe_rules",
        help="probe rule to scan; may be repeated",
    )
    return parser.parse_args()


def values_or_default(values: list[float] | list[int] | None, default: list[float] | list[int]):
    return values if values else default


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


def run_point(
    args: argparse.Namespace,
    momentum: float,
    seed: int,
    safety_factor: float,
    envelope_mode: str,
    probe_rule: str,
) -> dict[str, str]:
    command = [
        str(args.diagnostic_exe),
        "--calibration-samples",
        str(args.calibration_samples),
        "--trial-samples",
        str(args.trial_samples),
        "--momentum",
        f"{momentum:.17g}",
        "--particle-eta",
        f"{args.particle_eta:.17g}",
        "--particle-phi",
        f"{args.particle_phi:.17g}",
        "--seed",
        str(seed),
        "--safety-factor",
        f"{safety_factor:.17g}",
        "--envelope-mode",
        envelope_mode,
        "--probe-bins",
        str(args.probe_bins),
        "--probe-rule",
        probe_rule,
    ]
    row = parse_key_value_line(run_command(command))
    row["seed"] = str(seed)
    return row


def aggregate(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    groups: dict[tuple[str, str, str, str, str, str], list[dict[str, str]]] = {}
    for row in rows:
        groups.setdefault(
            (
                row["momentum"],
                row["safety_factor"],
                row["envelope_mode"],
                row["probe_rule"],
                row["particle_eta"],
                row["particle_phi"],
            ),
            [],
        ).append(row)

    summaries: list[dict[str, str]] = []
    for (
        momentum,
        safety_factor,
        envelope_mode,
        probe_rule,
        particle_eta,
        particle_phi,
    ), group_rows in sorted(
        groups.items(),
        key=lambda item: (
            float(item[0][0]),
            float(item[0][4]),
            float(item[0][5]),
            float(item[0][1]),
            item[0][2],
            item[0][3],
        ),
    ):
        failures = sum(int(row["failed_envelope_count"]) for row in group_rows)
        failed_runs = sum(1 for row in group_rows if int(row["failed_envelope_count"]) > 0)
        fallback_count = sum(int(row["fallback_envelope_count"]) for row in group_rows)
        observed_acceptances = [
            float(row["observed_acceptance_efficiency"]) for row in group_rows
        ]
        expected_acceptances = [
            float(row["expected_acceptance_efficiency"]) for row in group_rows
        ]
        max_ratios = [float(row["trial_calibration_max_ratio"]) for row in group_rows]
        local_max_ratios = [
            float(row["max_local_trial_calibration_ratio"]) for row in group_rows
        ]
        empty_cell_counts = [float(row["calibration_empty_cells"]) for row in group_rows]
        suspect_cell_counts = [float(row["suspect_cell_count"]) for row in group_rows]
        probe_evaluations = [float(row["probe_evaluations"]) for row in group_rows]
        summaries.append(
            {
                "momentum": momentum,
                "safety_factor": safety_factor,
                "envelope_mode": envelope_mode,
                "probe_rule": probe_rule,
                "particle_eta": particle_eta,
                "particle_phi": particle_phi,
                "runs": str(len(group_rows)),
                "failed_runs": str(failed_runs),
                "total_failed_envelope_count": str(failures),
                "total_fallback_envelope_count": str(fallback_count),
                "mean_calibration_empty_cells": f"{mean(empty_cell_counts):.17g}",
                "mean_suspect_cell_count": f"{mean(suspect_cell_counts):.17g}",
                "mean_probe_evaluations": f"{mean(probe_evaluations):.17g}",
                "mean_observed_acceptance_efficiency": f"{mean(observed_acceptances):.17g}",
                "mean_expected_acceptance_efficiency": f"{mean(expected_acceptances):.17g}",
                "max_trial_calibration_max_ratio": f"{max(max_ratios):.17g}",
                "max_local_trial_calibration_ratio": f"{max(local_max_ratios):.17g}",
            }
        )
    return summaries


def mean(values: list[float]) -> float:
    return sum(values) / float(len(values)) if values else math.nan


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise RuntimeError("no rows to write")
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_summary_md(path: Path, summaries: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Unweighting Envelope Scan",
        "",
        "| momentum | particle eta | safety factor | mode | probe rule | runs | failed runs | failed samples | fallback samples | empty cells | suspect cells | probe evals | mean observed eff. | mean expected eff. | max global ratio | max local ratio |",
        "| ---: | ---: | ---: | :--- | :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in summaries:
        lines.append(
            "| "
            f"{float(row['momentum']):.6g} | "
            f"{float(row['particle_eta']):.6g} | "
            f"{float(row['safety_factor']):.6g} | "
            f"{row['envelope_mode']} | "
            f"{row['probe_rule']} | "
            f"{row['runs']} | "
            f"{row['failed_runs']} | "
            f"{row['total_failed_envelope_count']} | "
            f"{row['total_fallback_envelope_count']} | "
            f"{float(row['mean_calibration_empty_cells']):.6g} | "
            f"{float(row['mean_suspect_cell_count']):.6g} | "
            f"{float(row['mean_probe_evaluations']):.6g} | "
            f"{float(row['mean_observed_acceptance_efficiency']):.6g} | "
            f"{float(row['mean_expected_acceptance_efficiency']):.6g} | "
            f"{float(row['max_trial_calibration_max_ratio']):.6g} | "
            f"{float(row['max_local_trial_calibration_ratio']):.6g} |"
        )
    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    args = parse_args()
    momenta = values_or_default(args.momenta, [1.0, 10.0, 100.0])
    seeds = values_or_default(args.seeds, [12345, 23456, 34567])
    safety_factors = values_or_default(args.safety_factors, [1.5, 3.0, 10.0])
    envelope_modes = values_or_default(args.envelope_modes, ["global", "cell", "probed-cell"])
    probe_rules = values_or_default(args.probe_rules, ["interior"])

    rows: list[dict[str, str]] = []
    for momentum in momenta:
        for safety_factor in safety_factors:
            for envelope_mode in envelope_modes:
                for probe_rule in probe_rules:
                    for seed in seeds:
                        rows.append(
                            run_point(
                                args,
                                float(momentum),
                                int(seed),
                                float(safety_factor),
                                str(envelope_mode),
                                str(probe_rule),
                            )
                        )

    summaries = aggregate(rows)
    write_csv(args.summary_csv, rows)
    write_summary_md(args.summary_md, summaries)

    print(f"wrote {args.summary_csv}")
    print(f"wrote {args.summary_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
