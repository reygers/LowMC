#!/usr/bin/env python3
"""Run a minimal HepMC3 afterburner plus output-summary quickstart."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run the LowMC HepMC3 afterburner on the checked-in charge-conserving "
            "fixture and summarize the produced output file."
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
        default=Path("adapters/hepmc3/tests/fixtures/charge_conserving_event.hepmc3"),
        help="read-only HepMC3 fixture or input file",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("build/hepmc3_quickstart"),
        help="directory where the quickstart output file is written",
    )
    parser.add_argument(
        "--enhancement-factor",
        type=float,
        default=10.0,
        help="LowMC photon enhancement factor for the quickstart run",
    )
    parser.add_argument("--multiplicity-seed", type=int, default=2)
    parser.add_argument("--shape-seed", type=int, default=67890)
    parser.add_argument("--min-shape-trials", type=int, default=1000)
    parser.add_argument("--shape-trials-per-photon", type=int, default=200)
    parser.add_argument(
        "--no-atomic-output",
        action="store_true",
        help=(
            "use the diagnostic afterburner's streaming output mode instead of "
            "the quickstart default --atomic-output mode"
        ),
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
    if result.stderr:
        print(result.stderr, file=sys.stderr, end="")
    return result.stdout


def require_file(path: Path, description: str) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"{description} not found: {path}")


def main() -> int:
    args = parse_args()
    afterburner = args.tool_dir / "lowmc_hepmc3_afterburner_diagnostic"
    summary = args.tool_dir / "lowmc_hepmc3_output_summary_diagnostic"
    output_path = args.output_dir / "charge_conserving_with_lowmc.hepmc3"

    require_file(afterburner, "afterburner executable")
    require_file(summary, "summary executable")
    require_file(args.input, "input HepMC3 file")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    afterburner_command = [
        str(afterburner),
        "--max-events",
        "1",
        "--multiplicity-seed",
        str(args.multiplicity_seed),
        "--shape-seed",
        str(args.shape_seed),
        "--min-shape-trials",
        str(args.min_shape_trials),
        "--shape-trials-per-photon",
        str(args.shape_trials_per_photon),
        "--enhancement-factor",
        f"{args.enhancement_factor:.17g}",
        str(args.input),
        str(output_path),
    ]
    if not args.no_atomic_output:
        afterburner_command.insert(1, "--atomic-output")
    summary_command = [str(summary), str(output_path)]

    print(f"input={args.input}")
    print(f"output={output_path}")
    print(f"atomic_output={0 if args.no_atomic_output else 1}")
    print("running_afterburner=1")
    afterburner_output = run_command(afterburner_command)
    print(afterburner_output, end="")
    print("running_summary=1")
    summary_output = run_command(summary_command)
    print(summary_output, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
