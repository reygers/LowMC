#!/usr/bin/env python3
"""Compare shape-unweighting candidate configurations across cases and seeds."""

from __future__ import annotations

import argparse
import csv
import math
import statistics
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

import validate_shape_unweighting as shape_validation


@dataclass(frozen=True)
class Case:
    name: str
    momentum: float
    particle_eta: float = 0.0
    particle_phi: float = math.pi / 2.0
    second_momentum: float | None = None
    second_particle_eta: float | None = None
    second_particle_phi: float | None = None
    topology: str = "two-pion"
    particle_file: Path | None = None


@dataclass(frozen=True)
class Candidate:
    name: str
    probe_policy: str
    boosted_safety_factor: float
    boosted_neighbor_scale: float
    grid_eta_bins: int
    grid_phi_bins: int
    grid_refinement_levels: int
    targeted_refinement_levels: int = 0
    targeted_refinement_neighbor_scale: float = 0.0


DEFAULT_TARGETED_CANDIDATE = Candidate(
    "targeted_l5_halo0.5_p11_to_p21_s3_local5",
    "boosted-cells",
    5.0,
    0.0,
    32,
    32,
    3,
    5,
    0.5,
)

COMPARISON_CANDIDATES = [
    Candidate("adaptive_p11_to_p21_s3_local15_n0", "boosted-cells", 15.0, 0.0, 32, 32, 3),
    Candidate("adaptive_48_l4_p11_to_p21_s3_local15_n0", "boosted-cells", 15.0, 0.0, 48, 48, 4),
    Candidate(
        "targeted_l4_p11_to_p21_s3_local5_n0",
        "boosted-cells",
        5.0,
        0.0,
        32,
        32,
        3,
        4,
        0.0,
    ),
    DEFAULT_TARGETED_CANDIDATE,
]

TARGETED_CANDIDATES = [DEFAULT_TARGETED_CANDIDATE]

TOPOLOGY_CANDIDATES = [
    Candidate("targeted_l4_p11_to_p21_s3_local5_n0", "boosted-cells", 5.0, 0.0, 32, 32, 3, 4, 0.0),
    Candidate("targeted_l4_halo0.5_p11_to_p21_s3_local5", "boosted-cells", 5.0, 0.0, 32, 32, 3, 4, 0.5),
    Candidate("targeted_l5_p11_to_p21_s3_local5_n0", "boosted-cells", 5.0, 0.0, 32, 32, 3, 5, 0.0),
    DEFAULT_TARGETED_CANDIDATE,
    Candidate("targeted_grid40_l3_l4_p11_to_p21_s3_local5_n0", "boosted-cells", 5.0, 0.0, 40, 40, 3, 4, 0.0),
    Candidate("adaptive_48_l4_p11_to_p21_s3_local15_n0", "boosted-cells", 15.0, 0.0, 48, 48, 4),
]

DEFAULT_CASES = [
    Case("mid_p1", 1.0),
    Case("mid_p1000", 1000.0),
    Case("forward_eta4_p100", 100.0, 4.0),
    Case(
        "asymmetric_mid_p100",
        100.0,
        0.35,
        1.1,
        35.0,
        -0.75,
        4.0,
    ),
    Case("event4_inspectable", 0.0, topology="event4"),
    Case("event4_text_inspectable", 0.0, particle_file=Path("fixtures/event4_particles.lowmc")),
    Case("forward4_inspectable", 0.0, topology="forward4"),
    Case("multi8_inspectable", 0.0, topology="multi8"),
    Case("clustered8_text_inspectable", 0.0, particle_file=Path("fixtures/clustered8_particles.lowmc")),
    Case(
        "incoming_outgoing4_text_inspectable",
        0.0,
        particle_file=Path("fixtures/incoming_outgoing4_particles.lowmc"),
    ),
    Case(
        "hepmc3_validation_text_inspectable",
        0.0,
        particle_file=Path("fixtures/hepmc3_validation_event_particles.lowmc"),
    ),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Validate selected shape-unweighting candidates across toy cases."
    )
    parser.add_argument(
        "summary_csv",
        nargs="?",
        type=Path,
        default=Path("build/shape_unweighting_candidate_validation.csv"),
    )
    parser.add_argument(
        "--summary-md",
        type=Path,
        default=Path("build/shape_unweighting_candidate_validation.md"),
    )
    parser.add_argument(
        "--truth-exe",
        type=Path,
        default=Path("build/lowmc_two_particle_diagnostic"),
    )
    parser.add_argument(
        "--shape-exe",
        type=Path,
        default=Path("build/lowmc_shape_unweighting_diagnostic"),
    )
    parser.add_argument("--work-dir", type=Path, default=Path("build/shape_candidate_validation"))
    parser.add_argument("--eta-bins", type=int, default=101)
    parser.add_argument("--phi-bins", type=int, default=121)
    parser.add_argument("--trials", type=int, default=50000)
    parser.add_argument("--calibration-trials", type=int, default=1000)
    parser.add_argument("--seed", type=int, action="append", dest="seeds")
    parser.add_argument(
        "--case",
        choices=tuple(case.name for case in DEFAULT_CASES),
        action="append",
        dest="cases",
        help="case to include; may be repeated; default includes all cases",
    )
    parser.add_argument(
        "--candidate-set",
        choices=("targeted", "comparison", "topology"),
        default="targeted",
        help="named candidate set to validate; comparison includes older baselines",
    )
    parser.add_argument("--safety-factor", type=float, default=3.0)
    parser.add_argument("--probe-bins", type=int, default=11)
    parser.add_argument("--boosted-probe-bins", type=int, default=21)
    parser.add_argument(
        "--failure-policy",
        choices=("diagnostic", "strict"),
        default="diagnostic",
    )
    parser.add_argument(
        "--gate-policy",
        choices=("off", "diagnostic", "strict"),
        default="off",
        help="optional toy-suite validation gate: off, mark rows, or exit nonzero on failure",
    )
    parser.add_argument(
        "--gate-max-trial-bound-ratio",
        type=float,
        help="maximum allowed trial/bound ratio for enabled validation gates; must be below 1",
    )
    args = parser.parse_args()
    if args.gate_policy != "off":
        if args.gate_max_trial_bound_ratio is None:
            parser.error("--gate-max-trial-bound-ratio is required when --gate-policy is enabled")
        if not (0.0 < args.gate_max_trial_bound_ratio < 1.0):
            parser.error("--gate-max-trial-bound-ratio must satisfy 0 < value < 1")
    return args


def candidates_for_set(name: str) -> list[Candidate]:
    if name == "targeted":
        return TARGETED_CANDIDATES
    if name == "comparison":
        return COMPARISON_CANDIDATES
    if name == "topology":
        return TOPOLOGY_CANDIDATES
    raise ValueError(f"unknown candidate set: {name}")


def cases_for_args(args: argparse.Namespace) -> list[Case]:
    if not args.cases:
        return DEFAULT_CASES
    requested = set(args.cases)
    return [case for case in DEFAULT_CASES if case.name in requested]


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


def add_case_arguments(command: list[str], case: Case) -> None:
    if case.particle_file is not None:
        command.extend(["--particles", str(case.particle_file)])
        return

    if case.topology != "two-pion":
        command.extend(["--topology", case.topology])
        return

    command.extend(
        [
            "--momentum",
            f"{case.momentum:.17g}",
            "--particle-eta",
            f"{case.particle_eta:.17g}",
            "--particle-phi",
            f"{case.particle_phi:.17g}",
        ]
    )
    if case.second_momentum is not None:
        command.extend(["--second-momentum", f"{case.second_momentum:.17g}"])
    if case.second_particle_eta is not None:
        command.extend(["--second-particle-eta", f"{case.second_particle_eta:.17g}"])
    if case.second_particle_phi is not None:
        command.extend(["--second-particle-phi", f"{case.second_particle_phi:.17g}"])


def run_truth(args: argparse.Namespace, case: Case, path: Path) -> None:
    command = [
        str(args.truth_exe),
        str(path),
        "--eta-bins",
        str(args.eta_bins),
        "--phi-bins",
        str(args.phi_bins),
    ]
    add_case_arguments(command, case)
    run_command(command)


def run_shape(
    args: argparse.Namespace,
    candidate: Candidate,
    case: Case,
    seed: int,
    path: Path,
) -> dict[str, str]:
    command = [
        str(args.shape_exe),
        "--trials",
        str(args.trials),
        "--calibration-trials",
        str(args.calibration_trials),
        "--seed",
        str(seed),
        "--safety-factor",
        f"{args.safety_factor:.17g}",
        "--probe-bins",
        str(args.probe_bins),
        "--probe-rule",
        "interior",
        "--probe-policy",
        candidate.probe_policy,
        "--boosted-probe-bins",
        str(args.boosted_probe_bins),
        "--boosted-safety-factor",
        f"{candidate.boosted_safety_factor:.17g}",
        "--boosted-neighbor-scale",
        f"{candidate.boosted_neighbor_scale:.17g}",
        "--grid-eta-bins",
        str(candidate.grid_eta_bins),
        "--grid-phi-bins",
        str(candidate.grid_phi_bins),
        "--grid-refinement-levels",
        str(candidate.grid_refinement_levels),
        "--grid-refine-fraction",
        "0.10000000000000001",
        "--grid-uniform-fraction",
        "0.10000000000000001",
        "--targeted-refinement-levels",
        str(candidate.targeted_refinement_levels),
        "--targeted-refinement-neighbor-scale",
        f"{candidate.targeted_refinement_neighbor_scale:.17g}",
        "--failure-policy",
        args.failure_policy,
        "--samples-csv",
        str(path),
    ]
    add_case_arguments(command, case)
    return parse_key_value_line(run_command(command))


def run_candidate_case_seed(
    args: argparse.Namespace,
    candidate: Candidate,
    case: Case,
    seed: int,
    truth_cache: dict[str, Path],
) -> dict[str, str]:
    args.work_dir.mkdir(parents=True, exist_ok=True)
    truth_csv = truth_cache.get(case.name)
    if truth_csv is None:
        truth_csv = args.work_dir / f"{case.name}_truth.csv"
        run_truth(args, case, truth_csv)
        truth_cache[case.name] = truth_csv

    sample_csv = args.work_dir / f"{candidate.name}_{case.name}_seed{seed}_samples.csv"
    diagnostics = run_shape(args, candidate, case, seed, sample_csv)
    etas, phis, truth = shape_validation.read_truth(truth_csv)
    samples = shape_validation.read_samples(sample_csv)
    metrics = shape_validation.shape_metrics(etas, phis, truth, samples)

    row = {
        "candidate": candidate.name,
        "case": case.name,
        "seed": str(seed),
        "momentum": diagnostics["momentum"],
        "particle_eta": diagnostics["particle_eta"],
        "particle_phi": diagnostics["particle_phi"],
        "second_momentum": diagnostics["second_momentum"],
        "second_particle_eta": diagnostics["second_particle_eta"],
        "second_particle_phi": diagnostics["second_particle_phi"],
        "topology": diagnostics.get("topology", case.topology),
        "particle_file": diagnostics.get(
            "particle_file",
            str(case.particle_file) if case.particle_file is not None else "none",
        ),
        "particle_count": diagnostics.get("particle_count", "2"),
        "status": diagnostics["status"],
        "trials": diagnostics["trials"],
        "accepted_candidates": diagnostics["accepted_candidates"],
        "failed_envelope_count": diagnostics["failed_envelope_count"],
        "observed_acceptance_efficiency": diagnostics["observed_acceptance_efficiency"],
        "expected_acceptance_efficiency": diagnostics["expected_acceptance_efficiency"],
        "envelope_integral": diagnostics["envelope_integral"],
        "max_trial_bound_ratio": diagnostics["max_trial_bound_ratio"],
        "suspect_cell_count": diagnostics["suspect_cell_count"],
        "promoted_cell_count": diagnostics["promoted_cell_count"],
        "fallback_count": diagnostics["fallback_count"],
        "grid_cells": diagnostics["grid_cells"],
        "probe_evaluations": diagnostics["probe_evaluations"],
        "truth_csv": str(truth_csv),
        "sample_csv": str(sample_csv),
    }
    for key, value in metrics.items():
        row[key] = f"{value:.17g}"
    add_gate_columns(row, args)
    return row


def mean_float(rows: list[dict[str, str]], key: str) -> float:
    return statistics.mean(float(row[key]) for row in rows)


def row_gate_reasons(row: dict[str, str], args: argparse.Namespace) -> list[str]:
    if args.gate_policy == "off":
        return []

    reasons: list[str] = []
    if row["status"] != "ok":
        reasons.append("non_ok_status")
    if int(row["failed_envelope_count"]) != 0:
        reasons.append("envelope_failures")
    if int(row["fallback_count"]) != 0:
        reasons.append("fallback_used")
    if args.gate_max_trial_bound_ratio is None:
        reasons.append("missing_ratio_threshold")
    elif float(row["max_trial_bound_ratio"]) >= args.gate_max_trial_bound_ratio:
        reasons.append("max_trial_bound_ratio")
    return reasons


def add_gate_columns(row: dict[str, str], args: argparse.Namespace) -> None:
    if args.gate_policy == "off":
        row["gate_policy"] = "off"
        row["gate_max_trial_bound_ratio"] = ""
        row["gate_status"] = "not_run"
        row["gate_reasons"] = ""
        return

    reasons = row_gate_reasons(row, args)
    row["gate_policy"] = args.gate_policy
    row["gate_max_trial_bound_ratio"] = f"{args.gate_max_trial_bound_ratio:.17g}"
    row["gate_status"] = "pass" if not reasons else "fail"
    row["gate_reasons"] = ",".join(reasons)


def summary_gate_reasons(row: dict[str, str], args: argparse.Namespace) -> list[str]:
    if args.gate_policy == "off":
        return []

    reasons: list[str] = []
    if int(row["non_ok_runs"]) != 0:
        reasons.append("non_ok_status")
    if int(row["failed_envelope_count"]) != 0:
        reasons.append("envelope_failures")
    if int(row["fallback_count"]) != 0:
        reasons.append("fallback_used")
    if args.gate_max_trial_bound_ratio is None:
        reasons.append("missing_ratio_threshold")
    elif float(row["max_trial_bound_ratio"]) >= args.gate_max_trial_bound_ratio:
        reasons.append("max_trial_bound_ratio")
    return reasons


def add_summary_gate_columns(row: dict[str, str], args: argparse.Namespace) -> None:
    if args.gate_policy == "off":
        row["gate_policy"] = "off"
        row["gate_max_trial_bound_ratio"] = ""
        row["gate_status"] = "not_run"
        row["gate_reasons"] = ""
        return

    reasons = summary_gate_reasons(row, args)
    row["gate_policy"] = args.gate_policy
    row["gate_max_trial_bound_ratio"] = f"{args.gate_max_trial_bound_ratio:.17g}"
    row["gate_status"] = "pass" if not reasons else "fail"
    row["gate_reasons"] = ",".join(reasons)


def aggregate(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    groups: dict[tuple[str, str], list[dict[str, str]]] = {}
    for row in rows:
        groups.setdefault((row["candidate"], row["case"]), []).append(row)

    summaries: list[dict[str, str]] = []
    for (candidate, case), group_rows in groups.items():
        summary = {
            "candidate": candidate,
            "case": case,
            "runs": str(len(group_rows)),
            "non_ok_runs": str(sum(1 for row in group_rows if row["status"] != "ok")),
            "failed_envelope_count": str(
                sum(int(row["failed_envelope_count"]) for row in group_rows)
            ),
            "accepted_candidates": str(sum(int(row["accepted_candidates"]) for row in group_rows)),
            "mean_observed_acceptance_efficiency": f"{mean_float(group_rows, 'observed_acceptance_efficiency'):.17g}",
            "mean_expected_acceptance_efficiency": f"{mean_float(group_rows, 'expected_acceptance_efficiency'):.17g}",
            "mean_envelope_integral": f"{mean_float(group_rows, 'envelope_integral'):.17g}",
            "max_trial_bound_ratio": f"{max(float(row['max_trial_bound_ratio']) for row in group_rows):.17g}",
            "mean_suspect_cell_count": f"{mean_float(group_rows, 'suspect_cell_count'):.17g}",
            "mean_promoted_cell_count": f"{mean_float(group_rows, 'promoted_cell_count'):.17g}",
            "fallback_count": str(sum(int(row["fallback_count"]) for row in group_rows)),
            "mean_fallback_count": f"{mean_float(group_rows, 'fallback_count'):.17g}",
            "mean_grid_cells": f"{mean_float(group_rows, 'grid_cells'):.17g}",
            "mean_probe_evaluations": f"{mean_float(group_rows, 'probe_evaluations'):.17g}",
            "mean_shape_l1": f"{mean_float(group_rows, 'shape_l1'):.17g}",
            "mean_eta_projection_rms": f"{mean_float(group_rows, 'eta_projection_rms'):.17g}",
            "mean_phi_projection_rms": f"{mean_float(group_rows, 'phi_projection_rms'):.17g}",
        }
        summaries.append(summary)
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
        "# Shape-Unweighting Candidate Validation",
        "",
        f"`{args.eta_bins} x {args.phi_bins}` truth grid, `{args.trials}` corrected shape trials per seed, `{args.calibration_trials}` calibration trials, seeds `{', '.join(str(seed) for seed in (args.seeds or [12345, 23456, 34567]))}`.",
        "",
        f"Candidate set: `{args.candidate_set}`.",
        "",
        f"Gate policy: `{args.gate_policy}`"
        + (
            f", max trial/bound threshold `{args.gate_max_trial_bound_ratio:.6g}`."
            if args.gate_max_trial_bound_ratio is not None
            else "."
        ),
        "",
        "| candidate | case | runs | gate | reasons | non-ok | failures | accepted | obs. eff. | exp. eff. | envelope | max trial/bound | suspect | promoted | fallback | grid cells | probe evals | L1 | eta RMS | phi RMS |",
        "| :--- | :--- | ---: | :--- | :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in summaries:
        lines.append(
            "| "
            f"{row['candidate']} | "
            f"{row['case']} | "
            f"{row['runs']} | "
            f"{row['gate_status']} | "
            f"{row['gate_reasons'] or '-'} | "
            f"{row['non_ok_runs']} | "
            f"{row['failed_envelope_count']} | "
            f"{row['accepted_candidates']} | "
            f"{float(row['mean_observed_acceptance_efficiency']):.6g} | "
            f"{float(row['mean_expected_acceptance_efficiency']):.6g} | "
            f"{float(row['mean_envelope_integral']):.6g} | "
            f"{float(row['max_trial_bound_ratio']):.6g} | "
            f"{float(row['mean_suspect_cell_count']):.6g} | "
            f"{float(row['mean_promoted_cell_count']):.6g} | "
            f"{row['fallback_count']} | "
            f"{float(row['mean_grid_cells']):.6g} | "
            f"{float(row['mean_probe_evaluations']):.6g} | "
            f"{float(row['mean_shape_l1']):.6g} | "
            f"{float(row['mean_eta_projection_rms']):.6g} | "
            f"{float(row['mean_phi_projection_rms']):.6g} |"
        )
    lines.extend(
        [
            "",
            "Shape metrics compare normalized accepted-candidate histograms with fixed-grid Low-intensity references. They are diagnostic statistics, not release thresholds.",
            "The gate, when enabled, checks only envelope diagnostics for this toy suite: zero non-ok runs, zero envelope failures, zero fallback uses, and a maximum trial/bound ratio below the configured threshold.",
        ]
    )
    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    args = parse_args()
    seeds = args.seeds if args.seeds else [12345, 23456, 34567]
    rows: list[dict[str, str]] = []
    truth_cache: dict[str, Path] = {}
    for candidate in candidates_for_set(args.candidate_set):
        for case in cases_for_args(args):
            for seed in seeds:
                rows.append(run_candidate_case_seed(args, candidate, case, int(seed), truth_cache))
    summaries = aggregate(rows)
    for summary in summaries:
        add_summary_gate_columns(summary, args)
    write_csv(args.summary_csv, rows)
    write_markdown(args.summary_md, summaries, args)
    print(f"wrote {args.summary_csv}")
    print(f"wrote {args.summary_md}")
    if args.gate_policy == "strict" and any(row["gate_status"] == "fail" for row in summaries):
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
