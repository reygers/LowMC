#!/usr/bin/env python3
"""Validate count-plus-shape diagnostics across controlled fixtures."""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Case:
    label: str
    validation_class: str
    momentum: float | None = None
    particles: Path | None = None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run lowmc_count_shape_diagnostic across controlled fixtures."
    )
    parser.add_argument(
        "summary_csv",
        nargs="?",
        type=Path,
        default=Path("build/count_shape_validation.csv"),
        help="per-run summary CSV to write",
    )
    parser.add_argument(
        "--aggregate-csv",
        type=Path,
        default=Path("build/count_shape_validation_aggregate.csv"),
        help="aggregate summary CSV to write",
    )
    parser.add_argument(
        "--summary-md",
        type=Path,
        default=Path("build/count_shape_validation.md"),
        help="Markdown summary to write",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("build/count_shape_validation"),
        help="directory for individual diagnostic CSV files",
    )
    parser.add_argument(
        "--diagnostic-exe",
        type=Path,
        default=Path("build/lowmc_count_shape_diagnostic"),
        help="path to lowmc_count_shape_diagnostic",
    )
    parser.add_argument("--seeds", type=int, default=3, help="number of seeds per case")
    parser.add_argument("--seed", type=int, default=12345, help="first seed")
    parser.add_argument("--oversampling-factor", type=float, default=100.0)
    parser.add_argument("--min-shape-trials", type=int, default=10000)
    parser.add_argument("--shape-trials-per-photon", type=int, default=1000)
    parser.add_argument(
        "--event4-particles",
        type=Path,
        default=Path("fixtures/event4_particles.lowmc"),
    )
    parser.add_argument(
        "--hepmc3-particles",
        type=Path,
        default=Path("fixtures/hepmc3_validation_event_particles.lowmc"),
    )
    return parser.parse_args()


def cases(args: argparse.Namespace) -> list[Case]:
    return [
        Case("toy_p1", "ordinary", momentum=1.0),
        Case("toy_p10", "ordinary", momentum=10.0),
        Case("toy_p100", "ordinary", momentum=100.0),
        Case("toy_p1000", "stress", momentum=1000.0),
        Case("event4_text", "ordinary", particles=args.event4_particles),
        Case("hepmc3_validation_text", "ordinary", particles=args.hepmc3_particles),
    ]


def case_args(case: Case) -> list[str]:
    if case.particles is not None:
        return ["--particles", str(case.particles)]
    if case.momentum is None:
        raise ValueError(f"case {case.label} has neither momentum nor particles")
    return ["--momentum", f"{case.momentum:.17g}"]


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


def run_case_seed(args: argparse.Namespace, case: Case, seed_index: int) -> dict[str, str]:
    multiplicity_seed = args.seed + 1000 * seed_index
    shape_seed = args.seed + 500 + 1000 * seed_index
    output_csv = args.output_dir / f"{case.label}_seed{seed_index}.csv"
    candidates_csv = args.output_dir / f"{case.label}_seed{seed_index}_candidates.csv"
    command = [
        str(args.diagnostic_exe),
        str(output_csv),
        *case_args(case),
        "--validation-class",
        case.validation_class,
        "--oversampling-factor",
        f"{args.oversampling_factor:.17g}",
        "--multiplicity-seed",
        str(multiplicity_seed),
        "--shape-seed",
        str(shape_seed),
        "--min-shape-trials",
        str(args.min_shape_trials),
        "--shape-trials-per-photon",
        str(args.shape_trials_per_photon),
        "--candidates-csv",
        str(candidates_csv),
    ]
    run_command(command)
    row = read_one_row(output_csv)
    row["case"] = case.label
    row["seed_index"] = str(seed_index)
    row["summary_csv"] = str(output_csv)
    row["candidates_csv"] = str(candidates_csv)
    return row


def to_float(row: dict[str, str], key: str) -> float:
    return float(row[key])


def to_int(row: dict[str, str], key: str) -> int:
    return int(float(row[key]))


def aggregate(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    output: list[dict[str, str]] = []
    for label in sorted({row["case"] for row in rows}):
        group = [row for row in rows if row["case"] == label]
        requested = [to_int(row, "requested_photons") for row in group]
        retained = [to_int(row, "retained_photons") for row in group]
        enough_count = sum(to_int(row, "enough_candidates") for row in group)
        envelope_count = sum(to_int(row, "envelope_ok") for row in group)
        status_count = sum(1 for row in group if row["status"] == "ok")
        max_ratio = max(to_float(row, "max_trial_bound_ratio") for row in group)
        total_requested = sum(requested)
        total_retained = sum(retained)
        output.append(
            {
                "case": label,
                "runs": str(len(group)),
                "validationClass": group[0]["validation_class"],
                "lambdaEvent": group[0]["lambda_event"],
                "oversamplingFactor": group[0]["oversampling_factor"],
                "poissonMean": group[0]["poisson_mean"],
                "requestedTotal": str(total_requested),
                "retainedTotal": str(total_retained),
                "requestedMin": str(min(requested)),
                "requestedMax": str(max(requested)),
                "allEnoughCandidates": "1" if enough_count == len(group) else "0",
                "allEnvelopeOk": "1" if envelope_count == len(group) else "0",
                "allStatusOk": "1" if status_count == len(group) else "0",
                "maxTrialBoundRatioMax": f"{max_ratio:.17g}",
            }
        )
    return output


def write_markdown(path: Path, aggregate_rows: list[dict[str, str]]) -> None:
    lines = [
        "# Count Plus Shape Validation",
        "",
        "| case | class | runs | lambda_event | poisson mean | requested | retained | envelope | enough | status | max ratio |",
        "| --- | --- | ---: | ---: | ---: | ---: | ---: | --- | --- | --- | ---: |",
    ]
    for row in aggregate_rows:
        lines.append(
            "| "
            f"{row['case']} | "
            f"{row['validationClass']} | "
            f"{row['runs']} | "
            f"{float(row['lambdaEvent']):.6g} | "
            f"{float(row['poissonMean']):.6g} | "
            f"{row['requestedTotal']} | "
            f"{row['retainedTotal']} | "
            f"{'ok' if row['allEnvelopeOk'] == '1' else 'watch'} | "
            f"{'ok' if row['allEnoughCandidates'] == '1' else 'watch'} | "
            f"{'ok' if row['allStatusOk'] == '1' else 'watch'} | "
            f"{float(row['maxTrialBoundRatioMax']):.3g} |"
        )
    lines.extend(["", "## Gate", ""])
    lines.append(
        "- Provisional pass: all runs have `status=ok`, `envelope_ok=1`, "
        "`enough_candidates=1`, and `max_trial_bound_ratio < 0.8`."
    )
    lines.append(
        "- The counts are diagnostic counts; with oversampling, each retained "
        "candidate carries analysis weight `1 / oversamplingFactor`."
    )
    lines.append(
        "- The stress row keeps `p = 1000 GeV` separate from ordinary fixtures."
    )
    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    args = parse_args()
    if args.seeds <= 0:
        print("--seeds must be positive", file=sys.stderr)
        return 1
    if args.oversampling_factor <= 0.0:
        print("--oversampling-factor must be positive", file=sys.stderr)
        return 1
    if args.min_shape_trials < 0 or args.shape_trials_per_photon < 0:
        print("shape trial counts must be non-negative", file=sys.stderr)
        return 1

    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.summary_csv.parent.mkdir(parents=True, exist_ok=True)
    args.aggregate_csv.parent.mkdir(parents=True, exist_ok=True)
    args.summary_md.parent.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, str]] = []
    for case in cases(args):
        for seed_index in range(args.seeds):
            row = run_case_seed(args, case, seed_index)
            rows.append(row)
            print(
                f"{case.label} seed={seed_index}: status={row['status']} "
                f"requested={row['requested_photons']} retained={row['retained_photons']} "
                f"envelope_ok={row['envelope_ok']} max_ratio={float(row['max_trial_bound_ratio']):.3g}"
            )

    fieldnames = [
        "case",
        "seed_index",
        "status",
        "validation_class",
        "validation_gate_ok",
        "lambda_event",
        "yield_method",
        "yield_log_kt_bins",
        "yield_grid_cells",
        "multiplicity_seed",
        "shape_seed",
        "oversampling_factor",
        "poisson_mean",
        "requested_photons",
        "retained_photons",
        "analysis_weight",
        "shape_status",
        "shape_trials",
        "shape_accepted_candidates",
        "enough_candidates",
        "envelope_ok",
        "failed_envelope_count",
        "max_trial_bound_ratio",
        "envelope_integral",
        "shape_grid_cells",
        "particle_count",
        "topology",
        "particle_file",
        "summary_csv",
        "candidates_csv",
    ]
    with args.summary_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    aggregate_rows = aggregate(rows)
    aggregate_fieldnames = [
        "case",
        "runs",
        "validationClass",
        "lambdaEvent",
        "oversamplingFactor",
        "poissonMean",
        "requestedTotal",
        "retainedTotal",
        "requestedMin",
        "requestedMax",
        "allEnoughCandidates",
        "allEnvelopeOk",
        "allStatusOk",
        "maxTrialBoundRatioMax",
    ]
    with args.aggregate_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=aggregate_fieldnames)
        writer.writeheader()
        writer.writerows(aggregate_rows)
    write_markdown(args.summary_md, aggregate_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
