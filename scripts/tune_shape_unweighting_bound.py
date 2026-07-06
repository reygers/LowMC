#!/usr/bin/env python3
"""Tune corrected shape-unweighting bounds for the boosted two-pion toy case."""

from __future__ import annotations

import argparse
import csv
import math
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Config:
    name: str
    probe_bins: int
    probe_rule: str
    safety_factor: float
    probe_policy: str = "fixed"
    boosted_probe_bins: int = 21
    boosted_safety_factor: float = 20.0
    boosted_neighbor_scale: float = 1.0
    grid_eta_bins: int = 32
    grid_phi_bins: int = 32
    grid_refinement_levels: int = 3
    grid_refine_fraction: float = 0.1
    grid_uniform_fraction: float = 0.1
    targeted_refinement_levels: int = 0
    targeted_refinement_neighbor_scale: float = 0.0


DEFAULT_CONFIGS = [
    Config("p11_interior_s3", 11, "interior", 3.0),
    Config("p11_interior_s10", 11, "interior", 10.0),
    Config("adaptive_p11_to_p21_s3_local15_n0", 11, "interior", 3.0, "boosted-cells", 21, 15.0, 0.0),
    Config("adaptive_p11_to_p21_s3_local15_n0.5", 11, "interior", 3.0, "boosted-cells", 21, 15.0, 0.5),
    Config("targeted_l1_p11_to_p21_s3_local15_n0", 11, "interior", 3.0, "boosted-cells", 21, 15.0, 0.0, 32, 32, 3, 0.1, 0.1, 1, 0.0),
    Config("targeted_l2_p11_to_p21_s3_local15_n0", 11, "interior", 3.0, "boosted-cells", 21, 15.0, 0.0, 32, 32, 3, 0.1, 0.1, 2, 0.0),
    Config("targeted_l3_p11_to_p21_s3_local15_n0", 11, "interior", 3.0, "boosted-cells", 21, 15.0, 0.0, 32, 32, 3, 0.1, 0.1, 3, 0.0),
    Config("targeted_l4_p11_to_p21_s3_local5_n0", 11, "interior", 3.0, "boosted-cells", 21, 5.0, 0.0, 32, 32, 3, 0.1, 0.1, 4, 0.0),
    Config("targeted_l4_p11_to_p21_s3_local15_n0", 11, "interior", 3.0, "boosted-cells", 21, 15.0, 0.0, 32, 32, 3, 0.1, 0.1, 4, 0.0),
    Config("adaptive_l4_p11_to_p21_s3_local15_n0", 11, "interior", 3.0, "boosted-cells", 21, 15.0, 0.0, 32, 32, 4, 0.1, 0.1),
    Config("adaptive_48_l4_p11_to_p21_s3_local15_n0", 11, "interior", 3.0, "boosted-cells", 21, 15.0, 0.0, 48, 48, 4, 0.1, 0.1),
    Config("adaptive_p11_to_p21_s3_local10", 11, "interior", 3.0, "boosted-cells", 21, 10.0, 1.0),
    Config("adaptive_p11_to_p21_s3_local15", 11, "interior", 3.0, "boosted-cells", 21, 15.0, 1.0),
    Config("adaptive_p11_to_p21_s3_local20", 11, "interior", 3.0, "boosted-cells", 21, 20.0),
    Config("adaptive_p11_to_p31_s3_local20", 11, "interior", 3.0, "boosted-cells", 31, 20.0),
    Config("p21_interior_s3", 21, "interior", 3.0),
    Config("p21_interior_s5", 21, "interior", 5.0),
    Config("p21_interior_s10", 21, "interior", 10.0),
    Config("p21_interior_s20", 21, "interior", 20.0),
    Config("p21_edge_s5", 21, "edge", 5.0),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Scan corrected shape-unweighting bounds for transverse p=1000 GeV."
    )
    parser.add_argument(
        "summary_csv",
        nargs="?",
        type=Path,
        default=Path("build/shape_unweighting_bound_tuning.csv"),
    )
    parser.add_argument(
        "--summary-md",
        type=Path,
        default=Path("build/shape_unweighting_bound_tuning.md"),
    )
    parser.add_argument(
        "--diagnostic-exe",
        type=Path,
        default=Path("build/lowmc_shape_unweighting_diagnostic"),
    )
    parser.add_argument("--trials", type=int, default=2000)
    parser.add_argument("--calibration-trials", type=int, default=1000)
    parser.add_argument("--seed", type=int, action="append", dest="seeds")
    parser.add_argument("--momentum", type=float, default=1000.0)
    parser.add_argument("--particle-eta", type=float, default=0.0)
    parser.add_argument("--particle-phi", type=float, default=math.pi / 2.0)
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


def run_config(args: argparse.Namespace, config: Config, seed: int) -> dict[str, str]:
    command = [
        str(args.diagnostic_exe),
        "--trials",
        str(args.trials),
        "--calibration-trials",
        str(args.calibration_trials),
        "--seed",
        str(seed),
        "--momentum",
        f"{args.momentum:.17g}",
        "--particle-eta",
        f"{args.particle_eta:.17g}",
        "--particle-phi",
        f"{args.particle_phi:.17g}",
        "--safety-factor",
        f"{config.safety_factor:.17g}",
        "--probe-bins",
        str(config.probe_bins),
        "--probe-rule",
        config.probe_rule,
        "--probe-policy",
        config.probe_policy,
        "--boosted-probe-bins",
        str(config.boosted_probe_bins),
        "--boosted-safety-factor",
        f"{config.boosted_safety_factor:.17g}",
        "--boosted-neighbor-scale",
        f"{config.boosted_neighbor_scale:.17g}",
        "--grid-eta-bins",
        str(config.grid_eta_bins),
        "--grid-phi-bins",
        str(config.grid_phi_bins),
        "--grid-refinement-levels",
        str(config.grid_refinement_levels),
        "--grid-refine-fraction",
        f"{config.grid_refine_fraction:.17g}",
        "--grid-uniform-fraction",
        f"{config.grid_uniform_fraction:.17g}",
        "--targeted-refinement-levels",
        str(config.targeted_refinement_levels),
        "--targeted-refinement-neighbor-scale",
        f"{config.targeted_refinement_neighbor_scale:.17g}",
        "--failure-policy",
        "diagnostic",
    ]
    row = parse_key_value_line(run_command(command))
    row["config"] = config.name
    row["seed"] = str(seed)
    return row


def mean(values: list[float]) -> float:
    return sum(values) / len(values) if values else math.nan


def aggregate(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    groups: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        groups.setdefault(row["config"], []).append(row)

    summaries: list[dict[str, str]] = []
    for name, group_rows in groups.items():
        failures = sum(int(row["failed_envelope_count"]) for row in group_rows)
        failed_runs = sum(1 for row in group_rows if int(row["failed_envelope_count"]) > 0)
        non_ok_runs = sum(1 for row in group_rows if row["status"] != "ok")
        summaries.append(
            {
                "config": name,
                "runs": str(len(group_rows)),
                "non_ok_runs": str(non_ok_runs),
                "failed_runs": str(failed_runs),
                "failed_envelope_count": str(failures),
                "accepted_candidates": str(
                    sum(int(row["accepted_candidates"]) for row in group_rows)
                ),
                "mean_observed_acceptance_efficiency": f"{mean([float(row['observed_acceptance_efficiency']) for row in group_rows]):.17g}",
                "mean_expected_acceptance_efficiency": f"{mean([float(row['expected_acceptance_efficiency']) for row in group_rows]):.17g}",
                "mean_envelope_integral": f"{mean([float(row['envelope_integral']) for row in group_rows]):.17g}",
                "max_trial_bound_ratio": f"{max(float(row['max_trial_bound_ratio']) for row in group_rows):.17g}",
                "mean_suspect_cell_count": f"{mean([float(row['suspect_cell_count']) for row in group_rows]):.17g}",
                "mean_promoted_cell_count": f"{mean([float(row['promoted_cell_count']) for row in group_rows]):.17g}",
                "mean_grid_cells": f"{mean([float(row['grid_cells']) for row in group_rows]):.17g}",
                "mean_probe_evaluations": f"{mean([float(row['probe_evaluations']) for row in group_rows]):.17g}",
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


def write_markdown(path: Path, summaries: list[dict[str, str]], args: argparse.Namespace) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Shape-Unweighting Bound Tuning",
        "",
        f"Transverse two-pion toy, `p = {args.momentum:g} GeV`, `{args.calibration_trials}` calibration trials and `{args.trials}` corrected-shape trials per seed.",
        "",
        "| config | runs | non-ok runs | failed runs | failed samples | accepted | obs. eff. | exp. eff. | envelope integral | max trial/bound | suspect cells | promoted cells | grid cells | probe evals |",
        "| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in summaries:
        lines.append(
            "| "
            f"{row['config']} | "
            f"{row['runs']} | "
            f"{row['non_ok_runs']} | "
            f"{row['failed_runs']} | "
            f"{row['failed_envelope_count']} | "
            f"{row['accepted_candidates']} | "
            f"{float(row['mean_observed_acceptance_efficiency']):.6g} | "
            f"{float(row['mean_expected_acceptance_efficiency']):.6g} | "
            f"{float(row['mean_envelope_integral']):.6g} | "
            f"{float(row['max_trial_bound_ratio']):.6g} | "
            f"{float(row['mean_suspect_cell_count']):.6g} | "
            f"{float(row['mean_promoted_cell_count']):.6g} | "
            f"{float(row['mean_grid_cells']):.6g} | "
            f"{float(row['mean_probe_evaluations']):.6g} |"
        )
    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    args = parse_args()
    seeds = args.seeds if args.seeds else [12345, 23456, 34567, 45678, 56789]
    rows: list[dict[str, str]] = []
    for config in DEFAULT_CONFIGS:
        for seed in seeds:
            rows.append(run_config(args, config, int(seed)))
    summaries = aggregate(rows)
    write_csv(args.summary_csv, rows)
    write_markdown(args.summary_md, summaries, args)
    print(f"wrote {args.summary_csv}")
    print(f"wrote {args.summary_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
