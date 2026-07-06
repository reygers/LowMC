#!/usr/bin/env python3
"""Compare API and legacy shape-unweighting diagnostics on toy cases."""

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
        description="Compare corrected API and legacy shape-unweighting diagnostics."
    )
    parser.add_argument(
        "summary_csv",
        nargs="?",
        type=Path,
        default=Path("build/shape_unweighting_diagnostic_comparison.csv"),
        help="per-run comparison CSV to write",
    )
    parser.add_argument(
        "--summary-md",
        type=Path,
        default=Path("build/shape_unweighting_diagnostic_comparison.md"),
    )
    parser.add_argument(
        "--api-exe",
        type=Path,
        default=Path("build/lowmc_shape_unweighting_diagnostic"),
    )
    parser.add_argument(
        "--legacy-exe",
        type=Path,
        default=Path("build/lowmc_unweighting_diagnostic"),
    )
    parser.add_argument("--calibration-trials", type=int, default=1000)
    parser.add_argument("--trials", type=int, default=2000)
    parser.add_argument("--seed", type=int, action="append", dest="seeds")
    parser.add_argument("--safety-factor", type=float, default=3.0)
    parser.add_argument("--probe-bins", type=int, default=11)
    parser.add_argument("--probe-rule", choices=("interior", "edge"), default="interior")
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


def run_api(args: argparse.Namespace, case: Case, seed: int) -> dict[str, str]:
    command = [
        str(args.api_exe),
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
        "diagnostic",
    ]
    return parse_key_value_line(run_command(command))


def run_legacy(args: argparse.Namespace, case: Case, seed: int) -> dict[str, str]:
    command = [
        str(args.legacy_exe),
        "--calibration-samples",
        str(args.calibration_trials),
        "--trial-samples",
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
        "--envelope-mode",
        "probed-cell",
        "--probe-bins",
        str(args.probe_bins),
        "--probe-rule",
        args.probe_rule,
    ]
    return parse_key_value_line(run_command(command))


def int_value(row: dict[str, str], key: str) -> int:
    return int(row[key])


def float_value(row: dict[str, str], key: str) -> float:
    return float(row[key])


def compare_run(
    args: argparse.Namespace, case: Case, seed: int, api: dict[str, str], legacy: dict[str, str]
) -> dict[str, str]:
    api_ratio = float_value(api, "max_trial_bound_ratio")
    legacy_local_ratio = float_value(legacy, "max_local_trial_calibration_ratio")
    legacy_ratio_scaled = legacy_local_ratio / args.safety_factor
    return {
        "case": case.name,
        "seed": str(seed),
        "momentum": api["momentum"],
        "particle_eta": api["particle_eta"],
        "api_status": api["status"],
        "legacy_status": legacy["envelope_status"],
        "api_accepted": api["accepted_candidates"],
        "legacy_accepted": legacy["accepted_photons"],
        "accepted_delta": str(
            int_value(api, "accepted_candidates") - int_value(legacy, "accepted_photons")
        ),
        "api_failed_envelope_count": api["failed_envelope_count"],
        "legacy_failed_envelope_count": legacy["failed_envelope_count"],
        "failed_envelope_delta": str(
            int_value(api, "failed_envelope_count")
            - int_value(legacy, "failed_envelope_count")
        ),
        "api_observed_acceptance_efficiency": api["observed_acceptance_efficiency"],
        "legacy_observed_acceptance_efficiency": legacy[
            "observed_acceptance_efficiency"
        ],
        "observed_acceptance_efficiency_delta": f"{float_value(api, 'observed_acceptance_efficiency') - float_value(legacy, 'observed_acceptance_efficiency'):.17g}",
        "api_expected_acceptance_efficiency": api["expected_acceptance_efficiency"],
        "legacy_expected_acceptance_efficiency": legacy[
            "expected_acceptance_efficiency"
        ],
        "expected_acceptance_efficiency_delta": f"{float_value(api, 'expected_acceptance_efficiency') - float_value(legacy, 'expected_acceptance_efficiency'):.17g}",
        "api_max_trial_bound_ratio": api["max_trial_bound_ratio"],
        "legacy_max_local_trial_calibration_ratio": legacy[
            "max_local_trial_calibration_ratio"
        ],
        "legacy_max_local_trial_calibration_ratio_scaled": f"{legacy_ratio_scaled:.17g}",
        "scaled_ratio_delta": f"{api_ratio - legacy_ratio_scaled:.17g}",
        "api_suspect_cell_count": api["suspect_cell_count"],
        "legacy_suspect_cell_count": legacy["suspect_cell_count"],
        "suspect_cell_delta": str(
            int_value(api, "suspect_cell_count") - int_value(legacy, "suspect_cell_count")
        ),
        "api_fallback_count": api["fallback_count"],
        "legacy_fallback_envelope_count": legacy["fallback_envelope_count"],
        "fallback_delta": str(
            int_value(api, "fallback_count")
            - int_value(legacy, "fallback_envelope_count")
        ),
        "api_probe_evaluations": api["probe_evaluations"],
        "legacy_probe_evaluations": legacy["probe_evaluations"],
    }


def mean(values: list[float]) -> float:
    return sum(values) / len(values) if values else math.nan


def aggregate(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    groups: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        groups.setdefault(row["case"], []).append(row)

    summaries: list[dict[str, str]] = []
    for case, group_rows in groups.items():
        summaries.append(
            {
                "case": case,
                "runs": str(len(group_rows)),
                "status_mismatches": str(
                    sum(
                        1
                        for row in group_rows
                        if (row["api_status"] == "ok") != (row["legacy_status"] == "held")
                    )
                ),
                "accepted_delta_sum": str(
                    sum(int_value(row, "accepted_delta") for row in group_rows)
                ),
                "failed_envelope_delta_sum": str(
                    sum(int_value(row, "failed_envelope_delta") for row in group_rows)
                ),
                "max_abs_observed_eff_delta": f"{max(abs(float_value(row, 'observed_acceptance_efficiency_delta')) for row in group_rows):.17g}",
                "max_abs_expected_eff_delta": f"{max(abs(float_value(row, 'expected_acceptance_efficiency_delta')) for row in group_rows):.17g}",
                "max_abs_scaled_ratio_delta": f"{max(abs(float_value(row, 'scaled_ratio_delta')) for row in group_rows):.17g}",
                "api_max_trial_bound_ratio": f"{max(float_value(row, 'api_max_trial_bound_ratio') for row in group_rows):.17g}",
                "legacy_max_local_ratio_scaled": f"{max(float_value(row, 'legacy_max_local_trial_calibration_ratio_scaled') for row in group_rows):.17g}",
                "mean_api_suspect_cells": f"{mean([float_value(row, 'api_suspect_cell_count') for row in group_rows]):.17g}",
                "mean_legacy_suspect_cells": f"{mean([float_value(row, 'legacy_suspect_cell_count') for row in group_rows]):.17g}",
                "fallback_delta_sum": str(
                    sum(int_value(row, "fallback_delta") for row in group_rows)
                ),
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
        "# Shape-Unweighting Diagnostic Comparison",
        "",
        f"`lowmc_shape_unweighting_diagnostic` versus `lowmc_unweighting_diagnostic`, `{args.probe_bins} x {args.probe_bins}` `{args.probe_rule}` probes, safety factor `{args.safety_factor:g}`.",
        "",
        "The API diagnostic uses the corrected cell-local envelope sampler. The legacy diagnostic samples from the original proposal and is retained as a bound/efficiency development comparison, so accepted counts and efficiencies are no longer expected to agree.",
        "",
        "The API diagnostic reports `max(w / M_c)` after the safety factor is applied. The legacy local ratio below is `max_local_trial_calibration_ratio / safetyFactor` so that the two maximum-ratio columns use the same scale.",
        "",
        "| case | runs | status mismatches | accepted delta | failed-envelope delta | max abs obs eff delta | max abs exp eff delta | max abs scaled ratio delta | API max trial/bound | legacy max scaled ratio | mean API suspect cells | mean legacy suspect cells | fallback delta |",
        "| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in summaries:
        lines.append(
            "| "
            f"{row['case']} | "
            f"{row['runs']} | "
            f"{row['status_mismatches']} | "
            f"{row['accepted_delta_sum']} | "
            f"{row['failed_envelope_delta_sum']} | "
            f"{float(row['max_abs_observed_eff_delta']):.6g} | "
            f"{float(row['max_abs_expected_eff_delta']):.6g} | "
            f"{float(row['max_abs_scaled_ratio_delta']):.6g} | "
            f"{float(row['api_max_trial_bound_ratio']):.6g} | "
            f"{float(row['legacy_max_local_ratio_scaled']):.6g} | "
            f"{float(row['mean_api_suspect_cells']):.6g} | "
            f"{float(row['mean_legacy_suspect_cells']):.6g} | "
            f"{row['fallback_delta_sum']} |"
        )
    lines.extend(
        [
            "",
            "Interpretation: maximum-ratio agreement is still useful for checking the same cell bounds. Accepted-count and efficiency differences are expected because the corrected API sampler draws from the normalized local envelope, while the legacy diagnostic draws from the original proposal.",
        ]
    )
    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    args = parse_args()
    seeds = args.seeds if args.seeds else [12345, 23456, 34567, 45678, 56789]

    rows: list[dict[str, str]] = []
    for case in DEFAULT_CASES:
        for seed in seeds:
            api = run_api(args, case, int(seed))
            legacy = run_legacy(args, case, int(seed))
            rows.append(compare_run(args, case, int(seed), api, legacy))

    summaries = aggregate(rows)
    write_csv(args.summary_csv, rows)
    write_summary_md(args.summary_md, summaries, args)
    print(f"wrote {args.summary_csv}")
    print(f"wrote {args.summary_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
