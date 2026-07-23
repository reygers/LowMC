#!/usr/bin/env python3
"""Run a predeclared multi-seed PHOTOS/LowMC soft-energy scan."""

from __future__ import annotations

import argparse
import csv
import math
import os
import statistics
import subprocess
import sys
from collections import defaultdict
from pathlib import Path


DEFAULT_ENERGY_MAXIMA = (0.06, 0.075, 0.10, 0.15, 0.25, 0.50)
DEFAULT_SEEDS = ((1802, 9373), (24680, 13579), (12345, 23456))


def parse_seed_pair(text: str) -> tuple[int, int]:
    fields = text.split(",")
    if len(fields) != 2:
        raise argparse.ArgumentTypeError("seed pair must have the form seed1,seed2")
    try:
        seed1, seed2 = (int(field) for field in fields)
    except ValueError as error:
        raise argparse.ArgumentTypeError("seed pair values must be integers") from error
    if not (0 <= seed1 <= 31327 and 0 <= seed2 <= 30080):
        raise argparse.ArgumentTypeError("seed pair is outside PHOTOS ranges")
    return seed1, seed2


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run nested photon-energy windows for fixed-order PHOTOS and LowMC."
    )
    parser.add_argument("summary_csv", type=Path)
    parser.add_argument("output_png", type=Path)
    parser.add_argument(
        "--benchmark",
        type=Path,
        default=Path("build-photos/lowmc_photos_zmumu_benchmark"),
    )
    parser.add_argument("--output-dir", type=Path, default=Path("build-photos/soft_scan"))
    parser.add_argument("--events", type=int, default=100000)
    parser.add_argument(
        "--energy-max",
        type=float,
        action="append",
        dest="energy_maxima",
        help="repeat to replace the predeclared energy maxima",
    )
    parser.add_argument(
        "--seed-pair",
        type=parse_seed_pair,
        action="append",
        dest="seed_pairs",
        help="repeat to replace the three predeclared PHOTOS seed pairs",
    )
    return parser.parse_args()


def run_one(
    args: argparse.Namespace, energy_max: float, seed1: int, seed2: int
) -> dict[str, str]:
    scratch_photons = args.output_dir / "scratch_photos.csv"
    scratch_lowmc = args.output_dir / "scratch_lowmc.csv"
    scratch_summary = args.output_dir / "scratch_summary.csv"
    command = [
        str(args.benchmark),
        str(scratch_photons),
        str(scratch_summary),
        "--events",
        str(args.events),
        "--seed1",
        str(seed1),
        "--seed2",
        str(seed2),
        "--energy-max",
        f"{energy_max:.17g}",
        "--lowmc-output",
        str(scratch_lowmc),
        "--lowmc-samples",
        "32",
        "--lowmc-refinement-levels",
        "3",
    ]
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        if result.stdout:
            print(result.stdout, file=sys.stdout, end="")
        if result.stderr:
            print(result.stderr, file=sys.stderr, end="")
        raise subprocess.CalledProcessError(result.returncode, command)
    with scratch_summary.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 1:
        raise RuntimeError("benchmark summary did not contain exactly one row")
    row = rows[0]
    print(
        f"energy_max={energy_max:.6g} seeds={seed1},{seed2} "
        f"accepted={row['accepted_photons']} photos_rate={float(row['photos_accepted_rate']):.8g} "
        f"lowmc={float(row['lowmc_lambda_event']):.8g}"
    )
    return row


def aggregate(rows: list[dict[str, str]]) -> list[dict[str, float | int]]:
    grouped: dict[float, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[float(row["energy_max_gev"])].append(row)

    output: list[dict[str, float | int]] = []
    for energy_max in sorted(grouped):
        group = grouped[energy_max]
        total_events = sum(int(row["events"]) for row in group)
        total_accepted = sum(int(row["accepted_photons"]) for row in group)
        photos_rate = total_accepted / total_events
        photos_error = math.sqrt(total_accepted) / total_events
        lowmc_values = [float(row["lowmc_lambda_event"]) for row in group]
        lowmc_rate = statistics.fmean(lowmc_values)
        lowmc_over_photos = lowmc_rate / photos_rate
        ratio_error = lowmc_over_photos * photos_error / photos_rate
        per_seed_rates = [float(row["photos_accepted_rate"]) for row in group]
        output.append(
            {
                "energy_max_gev": energy_max,
                "seed_count": len(group),
                "events": total_events,
                "accepted_photons": total_accepted,
                "photos_rate": photos_rate,
                "photos_rate_error": photos_error,
                "photos_seed_rate_stddev": statistics.stdev(per_seed_rates)
                if len(per_seed_rates) > 1
                else 0.0,
                "lowmc_rate": lowmc_rate,
                "lowmc_relative_integration_difference": max(
                    float(row["lowmc_relative_integration_difference"]) for row in group
                ),
                "lowmc_over_photos": lowmc_over_photos,
                "ratio_error": ratio_error,
            }
        )
    return output


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def plot_scan(path: Path, rows: list[dict[str, float | int]]) -> None:
    cache_directory = path.parent / ".matplotlib"
    cache_directory.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache_directory))
    os.environ.setdefault("XDG_CACHE_HOME", str(cache_directory))
    import matplotlib.pyplot as plt

    energy = [float(row["energy_max_gev"]) for row in rows]
    ratio = [float(row["lowmc_over_photos"]) for row in rows]
    error = [float(row["ratio_error"]) for row in rows]
    figure, axis = plt.subplots(figsize=(7.2, 4.8), constrained_layout=True)
    axis.errorbar(energy, ratio, yerr=error, fmt="o-", color="#2563eb", capsize=3)
    axis.axhline(1.0, color="#dc2626", linestyle="--", linewidth=1.2)
    axis.set_xscale("log")
    axis.set_xlabel(r"maximum photon energy $E_\gamma^{\max}$ [GeV]")
    axis.set_ylabel("LowMC / PHOTOS")
    axis.set_title(
        r"Fixed $Z\to\mu^+\mu^-$ soft-limit scan"
        + f"\n{int(rows[0]['events'])} total decays per energy window"
    )
    axis.grid(True, which="both", alpha=0.2)
    figure.savefig(path, dpi=180)
    plt.close(figure)


def main() -> int:
    args = parse_args()
    if args.events <= 0:
        print("--events must be positive", file=sys.stderr)
        return 2
    energy_maxima = tuple(args.energy_maxima or DEFAULT_ENERGY_MAXIMA)
    seed_pairs = tuple(args.seed_pairs or DEFAULT_SEEDS)
    if any(value <= 0.05 for value in energy_maxima):
        print("every --energy-max must exceed the fixed kT minimum of 0.05 GeV", file=sys.stderr)
        return 2
    if not args.benchmark.exists():
        print(f"benchmark executable not found: {args.benchmark}", file=sys.stderr)
        return 2

    args.output_dir.mkdir(parents=True, exist_ok=True)
    detailed: list[dict[str, str]] = []
    try:
        for energy_max in sorted(set(energy_maxima)):
            for seed1, seed2 in seed_pairs:
                detailed.append(run_one(args, energy_max, seed1, seed2))
    except (OSError, RuntimeError, subprocess.CalledProcessError) as error:
        print(f"Soft scan failed: {error}", file=sys.stderr)
        return 1

    aggregated = aggregate(detailed)
    write_csv(args.summary_csv, aggregated)
    detailed_path = args.summary_csv.with_name(args.summary_csv.stem + "_detailed.csv")
    write_csv(detailed_path, detailed)
    try:
        plot_scan(args.output_png, aggregated)
    except ImportError:
        print("matplotlib is required for the convergence plot", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
