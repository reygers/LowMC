#!/usr/bin/env python3
"""Compare the diagnostic leg-centred and tuned-grid angular proposals."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
import subprocess
import tempfile
import time


DEFAULT_FIXTURES = (
    "event4=fixtures/event4_particles.lowmc",
    "clustered8=fixtures/clustered8_particles.lowmc",
    "incoming_outgoing4=fixtures/incoming_outgoing4_particles.lowmc",
    "hepmc3_validation=fixtures/hepmc3_validation_event_particles.lowmc",
)
DEFAULT_SEEDS = (12345, 23456, 34567)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run a multi-seed weighted comparison of the analytic leg-centred "
            "proposal and the tuned mass-refined grid proposal."
        )
    )
    parser.add_argument("output", type=Path, help="aggregate output CSV")
    parser.add_argument(
        "--benchmark",
        type=Path,
        default=Path("build/lowmc_fixed_kt_weighted_samples"),
        help="path to lowmc_fixed_kt_weighted_samples",
    )
    parser.add_argument("--samples", type=int, default=100_000)
    parser.add_argument("--seed", type=int, action="append", dest="seeds")
    parser.add_argument(
        "--fixture",
        action="append",
        help="repeatable NAME=PATH particle fixture; defaults to four repository fixtures",
    )
    parser.add_argument("--leg-uniform-fraction", type=float, default=0.1)
    return parser.parse_args()


def parse_fixtures(values: list[str] | None) -> list[tuple[str, Path]]:
    fixtures = []
    for value in values or DEFAULT_FIXTURES:
        if "=" not in value:
            raise ValueError(f"fixture must have NAME=PATH form: {value}")
        name, path = value.split("=", 1)
        if not name or not path:
            raise ValueError(f"fixture must have NAME=PATH form: {value}")
        fixtures.append((name, Path(path)))
    return fixtures


def summarize_samples(path: Path) -> dict[str, float]:
    count = 0
    sum_weight = 0.0
    sum_weight2 = 0.0
    max_weight = 0.0
    inside = 0
    with path.open(newline="") as stream:
        for row in csv.DictReader(stream):
            weight = float(row["weight"])
            count += 1
            sum_weight += weight
            sum_weight2 += weight * weight
            max_weight = max(max_weight, weight)
            inside += int(row["insideAcceptance"])
    if count == 0:
        raise RuntimeError(f"benchmark wrote no samples: {path}")
    mean = sum_weight / count
    variance = max(0.0, sum_weight2 / count - mean * mean)
    return {
        "count": count,
        "sum_weight": sum_weight,
        "sum_weight2": sum_weight2,
        "mean_weight": mean,
        "standard_error": math.sqrt(variance / count),
        "effective_sample_fraction": (
            sum_weight * sum_weight / (count * sum_weight2) if sum_weight2 > 0.0 else 0.0
        ),
        "max_weight": max_weight,
        "inside_fraction": inside / count,
    }


def benchmark_command(
    executable: Path,
    sample_path: Path,
    fixture: Path,
    proposal: str,
    samples: int,
    seed: int,
    leg_uniform_fraction: float,
) -> list[str]:
    command = [
        str(executable),
        str(sample_path),
        "--particles",
        str(fixture),
        "--samples",
        str(samples),
        "--seed",
        str(seed),
        "--proposal",
        proposal,
    ]
    if proposal == "leg":
        command += ["--leg-uniform-fraction", str(leg_uniform_fraction)]
    else:
        command += [
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
    return command


def main() -> int:
    args = parse_args()
    if args.samples <= 0:
        raise ValueError("--samples must be positive")
    if not 0.0 < args.leg_uniform_fraction < 1.0:
        raise ValueError("--leg-uniform-fraction must satisfy 0 < f < 1")
    fixtures = parse_fixtures(args.fixture)
    seeds = args.seeds or list(DEFAULT_SEEDS)
    executable = args.benchmark.resolve()
    if not executable.is_file():
        raise FileNotFoundError(f"benchmark executable does not exist: {executable}")

    detailed_rows: list[dict[str, object]] = []
    totals: dict[tuple[str, str], dict[str, float]] = {}
    with tempfile.TemporaryDirectory(prefix="lowmc-leg-proposal-") as temporary:
        temporary_path = Path(temporary)
        for fixture_name, fixture_path in fixtures:
            resolved_fixture = fixture_path.resolve()
            if not resolved_fixture.is_file():
                raise FileNotFoundError(f"fixture does not exist: {resolved_fixture}")
            for proposal in ("leg", "grid"):
                for seed in seeds:
                    sample_path = temporary_path / f"{fixture_name}_{proposal}_{seed}.csv"
                    command = benchmark_command(
                        executable,
                        sample_path,
                        resolved_fixture,
                        proposal,
                        args.samples,
                        seed,
                        args.leg_uniform_fraction,
                    )
                    start = time.perf_counter()
                    subprocess.run(command, check=True, stdout=subprocess.DEVNULL,
                                   stderr=subprocess.DEVNULL)
                    elapsed = time.perf_counter() - start
                    summary = summarize_samples(sample_path)
                    detailed_rows.append(
                        {
                            "fixture": fixture_name,
                            "proposal": proposal,
                            "seed": seed,
                            "samples": summary["count"],
                            "mean_weight": summary["mean_weight"],
                            "standard_error": summary["standard_error"],
                            "effective_sample_fraction": summary[
                                "effective_sample_fraction"
                            ],
                            "max_weight": summary["max_weight"],
                            "inside_fraction": summary["inside_fraction"],
                            "wall_seconds": elapsed,
                        }
                    )
                    total = totals.setdefault(
                        (fixture_name, proposal),
                        {
                            "count": 0.0,
                            "sum_weight": 0.0,
                            "sum_weight2": 0.0,
                            "max_weight": 0.0,
                            "inside_count": 0.0,
                            "wall_seconds": 0.0,
                        },
                    )
                    total["count"] += summary["count"]
                    total["sum_weight"] += summary["sum_weight"]
                    total["sum_weight2"] += summary["sum_weight2"]
                    total["max_weight"] = max(total["max_weight"], summary["max_weight"])
                    total["inside_count"] += summary["inside_fraction"] * summary["count"]
                    total["wall_seconds"] += elapsed

    aggregate_rows = []
    for (fixture_name, proposal), total in totals.items():
        count = total["count"]
        mean = total["sum_weight"] / count
        variance = max(0.0, total["sum_weight2"] / count - mean * mean)
        aggregate_rows.append(
            {
                "fixture": fixture_name,
                "proposal": proposal,
                "seeds": len(seeds),
                "samples": int(count),
                "mean_weight": mean,
                "standard_error": math.sqrt(variance / count),
                "effective_sample_fraction": (
                    total["sum_weight"] ** 2 / (count * total["sum_weight2"])
                    if total["sum_weight2"] > 0.0
                    else 0.0
                ),
                "max_weight": total["max_weight"],
                "inside_fraction": total["inside_count"] / count,
                "wall_seconds": total["wall_seconds"],
            }
        )

    fields = (
        "fixture",
        "proposal",
        "seeds",
        "samples",
        "mean_weight",
        "standard_error",
        "effective_sample_fraction",
        "max_weight",
        "inside_fraction",
        "wall_seconds",
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(aggregate_rows)

    detailed_path = args.output.with_name(args.output.stem + "_detailed.csv")
    detailed_fields = (
        "fixture",
        "proposal",
        "seed",
        "samples",
        "mean_weight",
        "standard_error",
        "effective_sample_fraction",
        "max_weight",
        "inside_fraction",
        "wall_seconds",
    )
    with detailed_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=detailed_fields)
        writer.writeheader()
        writer.writerows(detailed_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
