#!/usr/bin/env python3
"""Compare corrected shape-unweighted samples with fixed-grid eta-phi truth."""

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
        description="Validate corrected shape-unweighted eta-phi samples against grid truth."
    )
    parser.add_argument(
        "summary_csv",
        nargs="?",
        type=Path,
        default=Path("build/shape_unweighting_validation_summary.csv"),
        help="summary CSV to write",
    )
    parser.add_argument(
        "--summary-md",
        type=Path,
        default=Path("build/shape_unweighting_validation_summary.md"),
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
    parser.add_argument("--work-dir", type=Path, default=Path("build/shape_validation"))
    parser.add_argument("--eta-bins", type=int, default=101)
    parser.add_argument("--phi-bins", type=int, default=121)
    parser.add_argument("--trials", type=int, default=50000)
    parser.add_argument("--calibration-trials", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--safety-factor", type=float, default=3.0)
    parser.add_argument("--probe-bins", type=int, default=11)
    parser.add_argument("--probe-rule", choices=("interior", "edge"), default="interior")
    parser.add_argument("--probe-policy", choices=("fixed", "boosted-cells"), default="fixed")
    parser.add_argument("--boosted-probe-bins", type=int, default=21)
    parser.add_argument("--boosted-safety-factor", type=float, default=20.0)
    parser.add_argument("--boosted-neighbor-scale", type=float, default=1.0)
    parser.add_argument("--grid-eta-bins", type=int, default=32)
    parser.add_argument("--grid-phi-bins", type=int, default=32)
    parser.add_argument("--grid-refinement-levels", type=int, default=3)
    parser.add_argument("--grid-refine-fraction", type=float, default=0.1)
    parser.add_argument("--grid-uniform-fraction", type=float, default=0.1)
    parser.add_argument("--targeted-refinement-levels", type=int, default=0)
    parser.add_argument("--targeted-refinement-neighbor-scale", type=float, default=0.0)
    parser.add_argument(
        "--failure-policy",
        choices=("diagnostic", "strict"),
        default="diagnostic",
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


def run_truth(args: argparse.Namespace, case: Case, path: Path) -> None:
    command = [
        str(args.truth_exe),
        str(path),
        "--momentum",
        f"{case.momentum:.17g}",
        "--particle-eta",
        f"{case.particle_eta:.17g}",
        "--particle-phi",
        f"{case.particle_phi:.17g}",
        "--eta-bins",
        str(args.eta_bins),
        "--phi-bins",
        str(args.phi_bins),
    ]
    run_command(command)


def run_shape(args: argparse.Namespace, case: Case, path: Path) -> dict[str, str]:
    command = [
        str(args.shape_exe),
        "--trials",
        str(args.trials),
        "--calibration-trials",
        str(args.calibration_trials),
        "--seed",
        str(args.seed),
        "--momentum",
        f"{case.momentum:.17g}",
        "--particle-eta",
        f"{case.particle_eta:.17g}",
        "--particle-phi",
        f"{case.particle_phi:.17g}",
        "--safety-factor",
        f"{args.safety_factor:.17g}",
        "--probe-bins",
        str(args.probe_bins),
        "--probe-rule",
        args.probe_rule,
        "--probe-policy",
        args.probe_policy,
        "--boosted-probe-bins",
        str(args.boosted_probe_bins),
        "--boosted-safety-factor",
        f"{args.boosted_safety_factor:.17g}",
        "--boosted-neighbor-scale",
        f"{args.boosted_neighbor_scale:.17g}",
        "--grid-eta-bins",
        str(args.grid_eta_bins),
        "--grid-phi-bins",
        str(args.grid_phi_bins),
        "--grid-refinement-levels",
        str(args.grid_refinement_levels),
        "--grid-refine-fraction",
        f"{args.grid_refine_fraction:.17g}",
        "--grid-uniform-fraction",
        f"{args.grid_uniform_fraction:.17g}",
        "--targeted-refinement-levels",
        str(args.targeted_refinement_levels),
        "--targeted-refinement-neighbor-scale",
        f"{args.targeted_refinement_neighbor_scale:.17g}",
        "--failure-policy",
        args.failure_policy,
        "--samples-csv",
        str(path),
    ]
    return parse_key_value_line(run_command(command))


def centers_to_edges(centers: list[float]) -> list[float]:
    if len(centers) < 2:
        raise ValueError("need at least two centers")
    edges = [centers[0] - 0.5 * (centers[1] - centers[0])]
    for index in range(len(centers) - 1):
        edges.append(0.5 * (centers[index] + centers[index + 1]))
    edges.append(centers[-1] + 0.5 * (centers[-1] - centers[-2]))
    return edges


def read_truth(path: Path) -> tuple[list[float], list[float], list[list[float]]]:
    rows: list[tuple[float, float, float]] = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            rows.append((float(row["eta"]), float(row["phi"]), float(row["S"])))
    if not rows:
        raise ValueError("empty truth CSV")

    etas = sorted({row[0] for row in rows})
    phis = sorted({row[1] for row in rows})
    eta_index = {eta: index for index, eta in enumerate(etas)}
    phi_index = {phi: index for index, phi in enumerate(phis)}
    values = [[0.0 for _ in phis] for _ in etas]
    for eta, phi, value in rows:
        values[eta_index[eta]][phi_index[phi]] = value
    return etas, phis, values


def read_samples(path: Path) -> list[tuple[float, float]]:
    samples: list[tuple[float, float]] = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            samples.append((float(row["eta"]), float(row["phi"])))
    return samples


def find_bin(edges: list[float], value: float) -> int | None:
    if value < edges[0] or value > edges[-1]:
        return None
    if value == edges[-1]:
        return len(edges) - 2
    lo = 0
    hi = len(edges) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if value < edges[mid]:
            hi = mid
        else:
            lo = mid
    return lo


def normalize(values: list[float]) -> list[float]:
    total = sum(values)
    if total <= 0.0:
        return [0.0 for _ in values]
    return [value / total for value in values]


def shape_metrics(
    etas: list[float],
    phis: list[float],
    truth: list[list[float]],
    samples: list[tuple[float, float]],
) -> dict[str, float]:
    eta_edges = centers_to_edges(etas)
    phi_edges = centers_to_edges(phis)
    truth_flat = [value for row in truth for value in row]
    truth_shape = normalize(truth_flat)

    counts = [[0.0 for _ in phis] for _ in etas]
    outside = 0
    for eta, phi in samples:
        eta_bin = find_bin(eta_edges, eta)
        phi_bin = find_bin(phi_edges, phi)
        if eta_bin is None or phi_bin is None:
            outside += 1
            continue
        counts[eta_bin][phi_bin] += 1.0

    sample_flat = [value for row in counts for value in row]
    sample_shape = normalize(sample_flat)
    l1 = sum(abs(sample - truth) for sample, truth in zip(sample_shape, truth_shape))

    truth_eta = normalize([sum(row) for row in truth])
    sample_eta = normalize([sum(row) for row in counts])
    truth_phi = normalize([sum(truth[eta_i][phi_i] for eta_i in range(len(etas))) for phi_i in range(len(phis))])
    sample_phi = normalize([sum(counts[eta_i][phi_i] for eta_i in range(len(etas))) for phi_i in range(len(phis))])

    eta_rms = math.sqrt(
        sum((sample - truth) ** 2 for sample, truth in zip(sample_eta, truth_eta)) /
        float(len(truth_eta))
    )
    phi_rms = math.sqrt(
        sum((sample - truth) ** 2 for sample, truth in zip(sample_phi, truth_phi)) /
        float(len(truth_phi))
    )
    eta_max = max(abs(sample - truth) for sample, truth in zip(sample_eta, truth_eta))
    phi_max = max(abs(sample - truth) for sample, truth in zip(sample_phi, truth_phi))

    return {
        "shape_l1": l1,
        "eta_projection_rms": eta_rms,
        "phi_projection_rms": phi_rms,
        "eta_projection_max_abs": eta_max,
        "phi_projection_max_abs": phi_max,
        "outside_samples": float(outside),
    }


def run_case(args: argparse.Namespace, case: Case) -> dict[str, str]:
    args.work_dir.mkdir(parents=True, exist_ok=True)
    truth_csv = args.work_dir / f"{case.name}_truth.csv"
    sample_csv = args.work_dir / f"{case.name}_shape_samples.csv"
    run_truth(args, case, truth_csv)
    diagnostics = run_shape(args, case, sample_csv)
    etas, phis, truth = read_truth(truth_csv)
    samples = read_samples(sample_csv)
    metrics = shape_metrics(etas, phis, truth, samples)

    row = {
        "case": case.name,
        "momentum": diagnostics["momentum"],
        "particle_eta": diagnostics["particle_eta"],
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
        "probe_evaluations": diagnostics["probe_evaluations"],
        "truth_csv": str(truth_csv),
        "sample_csv": str(sample_csv),
    }
    for key, value in metrics.items():
        row[key] = f"{value:.17g}"
    return row


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise RuntimeError("no rows to write")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def write_markdown(path: Path, rows: list[dict[str, str]], args: argparse.Namespace) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Shape-Unweighting Validation Summary",
        "",
        f"`{args.eta_bins} x {args.phi_bins}` truth grid, `{args.trials}` corrected shape trials, `{args.probe_bins} x {args.probe_bins}` `{args.probe_rule}` probes, policy `{args.probe_policy}`, safety factor `{args.safety_factor:g}`, boosted neighbor scale `{args.boosted_neighbor_scale:g}`.",
        "",
        "| case | status | accepted | failures | obs. eff. | envelope integral | max trial/bound | promoted cells | L1 shape | eta RMS | phi RMS |",
        "| :--- | :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in rows:
        lines.append(
            "| "
            f"{row['case']} | "
            f"{row['status']} | "
            f"{row['accepted_candidates']} | "
            f"{row['failed_envelope_count']} | "
            f"{float(row['observed_acceptance_efficiency']):.6g} | "
            f"{float(row['envelope_integral']):.6g} | "
            f"{float(row['max_trial_bound_ratio']):.6g} | "
            f"{row['promoted_cell_count']} | "
            f"{float(row['shape_l1']):.6g} | "
            f"{float(row['eta_projection_rms']):.6g} | "
            f"{float(row['phi_projection_rms']):.6g} |"
        )
    lines.extend(
        [
            "",
            "The L1 and projection RMS values compare normalized accepted-candidate histograms with the fixed-grid Low intensity shape. They are diagnostic statistics, not final release thresholds.",
        ]
    )
    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    args = parse_args()
    rows = [run_case(args, case) for case in DEFAULT_CASES]
    write_csv(args.summary_csv, rows)
    write_markdown(args.summary_md, rows, args)
    print(f"wrote {args.summary_csv}")
    print(f"wrote {args.summary_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
