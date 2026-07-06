#!/usr/bin/env python3
"""Run an opt-in higher-statistics HepMC3 count/yield validation."""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
import time
from pathlib import Path


AGGREGATE_RE = re.compile(
    r"aggregate cases=(?P<cases>\d+) events=(?P<events>\d+) "
    r"expected_mean=(?P<expected_mean>[^ ]+) requested_photons=(?P<requested_photons>\d+) "
    r"difference=(?P<difference>[^ ]+) tolerance=(?P<tolerance>[^ \n]+)"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run a slower opt-in LowMC HepMC3 count/yield validation. "
            "This checks Poisson multiplicity and file bookkeeping, not the "
            "physical event-generator distribution."
        )
    )
    parser.add_argument(
        "--tool-dir",
        type=Path,
        default=Path("build-hepmc3/adapters/hepmc3"),
        help="directory containing the HepMC3 LowMC command-line tools",
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=None,
        help="optional pre-existing HepMC3 input file; otherwise generate a sample",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("build/hepmc3_count_yield_highstat"),
        help="directory where generated inputs and temporary outputs are written",
    )
    parser.add_argument(
        "--seeds",
        type=int,
        nargs="+",
        default=[1, 2, 3, 4, 5],
        help="multiplicity seeds scanned by the underlying count/yield checker",
    )
    parser.add_argument(
        "--enhancement-factors",
        type=float,
        nargs="+",
        default=[10.0, 25.0, 50.0],
        help="enhancement factors scanned by the underlying count/yield checker",
    )
    parser.add_argument(
        "--max-events",
        type=int,
        default=5,
        help="number of input events to process; use 0 for all events",
    )
    parser.add_argument(
        "--max-poisson-sigma",
        type=float,
        default=3.0,
        help="aggregate Poisson tolerance in sqrt(expected) units",
    )
    parser.add_argument("--shape-seed", type=int, default=67890)
    parser.add_argument("--min-shape-trials", type=int, default=800)
    parser.add_argument("--shape-trials-per-photon", type=int, default=120)
    return parser.parse_args()


def run_command(command: list[str]) -> str:
    result = subprocess.run(command, capture_output=True, text=True)
    if result.stdout:
        print(result.stdout, file=sys.stdout, end="")
    if result.stderr:
        print(result.stderr, file=sys.stderr, end="")
    if result.returncode != 0:
        raise subprocess.CalledProcessError(result.returncode, command)
    return result.stdout


def require_file(path: Path, description: str) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"{description} not found: {path}")


def generated_input_path(args: argparse.Namespace) -> Path:
    generator = args.tool_dir / "lowmc_hepmc3_generate_realistic_sample"
    require_file(generator, "generated-sample helper")
    input_path = args.output_dir / "generated_realistic_sample_highstat.hepmc3"
    run_command([str(generator), str(input_path)])
    return input_path


def main() -> int:
    args = parse_args()
    if args.max_events < 0:
        raise ValueError("--max-events must be non-negative")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    input_path = args.input if args.input is not None else generated_input_path(args)
    require_file(input_path, "input HepMC3 file")

    smoke_script = Path(__file__).resolve().with_name("hepmc3_count_yield_smoke.py")
    require_file(smoke_script, "count/yield checker")

    command = [
        sys.executable,
        str(smoke_script),
        "--tool-dir",
        str(args.tool_dir),
        "--input",
        str(input_path),
        "--output-dir",
        str(args.output_dir / "runs"),
        "--seeds",
        *[str(seed) for seed in args.seeds],
        "--enhancement-factors",
        *[f"{factor:.17g}" for factor in args.enhancement_factors],
        "--max-events",
        str(args.max_events),
        "--max-poisson-sigma",
        f"{args.max_poisson_sigma:.17g}",
        "--shape-seed",
        str(args.shape_seed),
        "--min-shape-trials",
        str(args.min_shape_trials),
        "--shape-trials-per-photon",
        str(args.shape_trials_per_photon),
    ]

    started = time.monotonic()
    checker_output = run_command(command)
    elapsed = time.monotonic() - started

    aggregate = None
    for match in AGGREGATE_RE.finditer(checker_output):
        aggregate = match.groupdict()
    if aggregate is None:
        raise AssertionError("count/yield checker did not print an aggregate summary")

    print(
        "highstat "
        f"input={input_path.name} cases={aggregate['cases']} events={aggregate['events']} "
        f"expected_mean={aggregate['expected_mean']} "
        f"requested_photons={aggregate['requested_photons']} "
        f"difference={aggregate['difference']} tolerance={aggregate['tolerance']} "
        f"runtime_seconds={elapsed:.3f}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
