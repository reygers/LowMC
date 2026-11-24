#!/usr/bin/env python3
"""Run a small HepMC3 afterburner count/yield consistency smoke test."""

from __future__ import annotations

import argparse
import math
import re
import subprocess
import sys
from pathlib import Path


KEY_VALUE_RE = re.compile(r"([A-Za-z_][A-Za-z0-9_]*)=([^ \n]+)")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run the LowMC HepMC3 afterburner for deterministic seeds and check "
            "that requested/added photon counts are consistent with reported "
            "lambda_event and enhancementFactor at smoke-test precision."
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
        nargs="+",
        default=[Path("adapters/hepmc3/tests/fixtures/charge_conserving_event.hepmc3")],
        help="one or more read-only HepMC3 fixtures or input files",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("build/hepmc3_count_yield_smoke"),
        help="directory where temporary afterburner outputs are written",
    )
    parser.add_argument(
        "--seeds",
        type=int,
        nargs="+",
        default=[1, 2, 3, 4],
        help="multiplicity seeds to scan for each enhancement factor",
    )
    parser.add_argument(
        "--enhancement-factors",
        type=float,
        nargs="+",
        default=[1.0, 10.0],
        help="enhancement factors to scan",
    )
    parser.add_argument("--shape-seed", type=int, default=67890)
    parser.add_argument("--min-shape-trials", type=int, default=1000)
    parser.add_argument("--shape-trials-per-photon", type=int, default=200)
    parser.add_argument(
        "--max-events",
        type=int,
        default=0,
        help="maximum events passed to the afterburner and summary diagnostics; 0 means all",
    )
    parser.add_argument(
        "--max-poisson-sigma",
        type=float,
        default=5.0,
        help="loose aggregate Poisson smoke-test tolerance in sqrt(expected) units",
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


def parse_key_values(line: str) -> dict[str, str]:
    return {match.group(1): match.group(2) for match in KEY_VALUE_RE.finditer(line)}


def diagnostic_lines(output: str, required_key: str) -> list[dict[str, str]]:
    lines = []
    for line in output.splitlines():
        values = parse_key_values(line)
        if required_key in values:
            lines.append(values)
    if not lines:
        raise ValueError(f"could not find diagnostic line containing {required_key!r}")
    return lines


def as_int(values: dict[str, str], key: str) -> int:
    return int(values[key])


def as_float(values: dict[str, str], key: str) -> float:
    return float(values[key])


def output_name(input_path: Path, enhancement_factor: float, seed: int) -> str:
    factor_text = f"{enhancement_factor:.12g}".replace(".", "p")
    return f"{input_path.stem}_f{factor_text}_seed{seed}.hepmc3"


def check_run(
    afterburner: Path,
    summary: Path,
    input_path: Path,
    output_dir: Path,
    enhancement_factor: float,
    seed: int,
    args: argparse.Namespace,
) -> tuple[int, float, int]:
    output_path = output_dir / output_name(input_path, enhancement_factor, seed)
    afterburner_command = [
        str(afterburner),
        "--failure-policy",
        "skip",
        "--multiplicity-seed",
        str(seed),
        "--shape-seed",
        str(args.shape_seed),
        "--min-shape-trials",
        str(args.min_shape_trials),
        "--shape-trials-per-photon",
        str(args.shape_trials_per_photon),
        "--enhancement-factor",
        f"{enhancement_factor:.17g}",
        str(input_path),
        str(output_path),
    ]
    summary_command = [str(summary), str(output_path)]
    if args.max_events > 0:
        afterburner_command[1:1] = ["--max-events", str(args.max_events)]
        summary_command[1:1] = ["--max-events", str(args.max_events)]

    afterburner_output = run_command(afterburner_command)
    summary_output = run_command(summary_command)
    afterburner_rows = diagnostic_lines(afterburner_output, "photon_status")
    summary_rows = diagnostic_lines(summary_output, "lowmc_photons")
    if len(afterburner_rows) != len(summary_rows):
        raise AssertionError(
            f"event-count mismatch for {input_path}: "
            f"afterburner_rows={len(afterburner_rows)} summary_rows={len(summary_rows)}"
        )

    run_requested = 0
    run_expected_mean = 0.0
    expected_mode = "unweighted" if enhancement_factor == 1.0 else "enhanced"
    for afterburner_values, summary_values in zip(afterburner_rows, summary_rows):
        event_index = as_int(afterburner_values, "event_index")
        summary_event_index = as_int(summary_values, "event_index")
        if summary_event_index != event_index:
            raise AssertionError(
                f"event-index mismatch for {input_path}: "
                f"afterburner={event_index} summary={summary_event_index}"
            )

        status = afterburner_values["photon_status"]
        requested = as_int(afterburner_values, "requested_photons")
        retained = as_int(afterburner_values, "retained_photons")
        added = as_int(afterburner_values, "photons_added")
        lowmc_photons = as_int(summary_values, "lowmc_photons")
        lambda_event = as_float(afterburner_values, "lambda_event")
        reported_factor = as_float(afterburner_values, "enhancement_factor")

        if abs(reported_factor - enhancement_factor) > 1.0e-12:
            raise AssertionError(
                f"reported enhancement factor {reported_factor} does not match "
                f"{enhancement_factor} for {input_path} event {event_index}"
            )

        if status in {"generated", "ok_zero_photons"}:
            if not math.isfinite(lambda_event) or lambda_event <= 0.0:
                raise AssertionError(
                    f"invalid lambda_event for {input_path} event {event_index}: {lambda_event}"
                )
            if retained != requested or added != requested or lowmc_photons != requested:
                raise AssertionError(
                    "count mismatch for "
                    f"input={input_path.name} factor={enhancement_factor} seed={seed} "
                    f"event={event_index}: requested={requested} retained={retained} "
                    f"added={added} summary={lowmc_photons}"
                )
            if summary_values["photon_mode"] != ("none" if requested == 0 else expected_mode):
                raise AssertionError(
                    f"unexpected summary mode for input={input_path.name} "
                    f"factor={enhancement_factor} seed={seed} event={event_index}: "
                    f"{summary_values['photon_mode']}"
                )
            run_requested += requested
            run_expected_mean += enhancement_factor * lambda_event
        elif status.startswith("skipped_"):
            if requested != 0 or retained != 0 or added != 0 or lowmc_photons != 0:
                raise AssertionError(
                    "skipped event should not contain photons for "
                    f"input={input_path.name} factor={enhancement_factor} seed={seed} "
                    f"event={event_index}: requested={requested} retained={retained} "
                    f"added={added} summary={lowmc_photons}"
                )
            if summary_values["photon_mode"] != "none":
                raise AssertionError(
                    f"skipped event should have photon_mode=none for {input_path.name} "
                    f"event={event_index}, got {summary_values['photon_mode']}"
                )
        else:
            raise AssertionError(
                f"unexpected afterburner status for {input_path.name} "
                f"seed={seed} event={event_index}: {status}"
            )

    print(
        "row "
        f"input={input_path.name} "
        f"enhancement_factor={enhancement_factor:.12g} seed={seed} "
        f"events={len(afterburner_rows)} expected_mean={run_expected_mean:.12g} "
        f"requested_photons={run_requested}"
    )
    return run_requested, run_expected_mean, len(afterburner_rows)


def main() -> int:
    args = parse_args()
    afterburner = args.tool_dir / "lowmc_hepmc3_afterburner_diagnostic"
    summary = args.tool_dir / "lowmc_hepmc3_output_summary_diagnostic"

    require_file(afterburner, "afterburner executable")
    require_file(summary, "summary executable")
    if args.max_events < 0:
        raise ValueError("--max-events must be non-negative")
    for input_path in args.input:
        require_file(input_path, "input HepMC3 file")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    total_requested = 0
    total_expected = 0.0
    total_events = 0
    total_cases = 0
    for input_path in args.input:
        for enhancement_factor in args.enhancement_factors:
            if enhancement_factor <= 0.0 or not math.isfinite(enhancement_factor):
                raise ValueError("enhancement factors must be positive and finite")
            for seed in args.seeds:
                requested, expected, events = check_run(
                    afterburner,
                    summary,
                    input_path,
                    args.output_dir,
                    enhancement_factor,
                    seed,
                    args,
                )
                total_requested += requested
                total_expected += expected
                total_events += events
                total_cases += 1

    sigma = math.sqrt(total_expected)
    tolerance = max(3.0, args.max_poisson_sigma * sigma)
    difference = abs(total_requested - total_expected)
    print(
        "aggregate "
        f"cases={total_cases} events={total_events} "
        f"expected_mean={total_expected:.12g} requested_photons={total_requested} "
        f"difference={difference:.12g} tolerance={tolerance:.12g}"
    )
    if difference > tolerance:
        raise AssertionError(
            f"aggregate requested photons {total_requested} differ from expected "
            f"{total_expected:.12g} by more than {tolerance:.12g}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
