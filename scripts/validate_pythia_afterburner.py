#!/usr/bin/env python3
"""Run bounded Pythia-produced HepMC3 afterburner validation."""

from __future__ import annotations

import argparse
import csv
import math
import re
import subprocess
import sys
import time
from pathlib import Path


KEY_VALUE_RE = re.compile(r"([A-Za-z_][A-Za-z0-9_]*)=([^ \n]+)")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate a bounded Pythia8 -> HepMC3 sample, run the LowMC "
            "HepMC3 afterburner on it, run the read-only output-summary "
            "diagnostic, and print a compact validation summary. This is an "
            "opt-in manual validation helper, not a default test-suite path."
        )
    )
    parser.add_argument(
        "--tool-dir",
        type=Path,
        default=Path("build-pythia8-hepmc3/adapters/hepmc3"),
        help="directory containing lowmc_hepmc3_afterburner_diagnostic and summary tool",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("build-pythia8-hepmc3/pythia_afterburner_validation"),
        help="directory for generated input, afterburner output, logs, and CSV summary",
    )
    parser.add_argument("--pythia-prefix", type=Path, default=Path.home() / "pythia8317")
    parser.add_argument("--hepmc3-prefix", type=Path, default=Path.home() / "hepmc3")
    parser.add_argument("--events", type=int, default=5)
    parser.add_argument("--pythia-seed", type=int, default=24680)
    parser.add_argument("--ecm", type=float, default=13000.0)
    parser.add_argument("--enhancement-factor", type=float, default=10.0)
    parser.add_argument("--multiplicity-seed", type=int, default=13579)
    parser.add_argument("--shape-seed", type=int, default=97531)
    parser.add_argument("--kt-min", type=float, default=0.002)
    parser.add_argument("--kt-max", type=float, default=0.2)
    parser.add_argument("--eta-min", type=float, default=-5.0)
    parser.add_argument("--eta-max", type=float, default=5.0)
    parser.add_argument("--phi-min", type=float, default=0.0)
    parser.add_argument("--phi-max", type=float, default=2.0 * math.pi)
    parser.add_argument(
        "--cxx",
        default="c++",
        help="C++ compiler used by generate_pythia_hepmc3_sample.py",
    )
    parser.add_argument(
        "--no-fail-on-mismatch",
        action="store_true",
        help="report summary mismatches without returning a failing exit code",
    )
    return parser.parse_args()


def require_file(path: Path, description: str) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"{description} not found: {path}")


def parse_key_values(line: str) -> dict[str, str]:
    return {match.group(1): match.group(2) for match in KEY_VALUE_RE.finditer(line)}


def rows_with_key(output: str, key: str) -> list[dict[str, str]]:
    return [
        values
        for line in output.splitlines()
        if key in (values := parse_key_values(line))
    ]


def run_logged(command: list[str], log_path: Path) -> tuple[str, float]:
    start = time.perf_counter()
    result = subprocess.run(command, capture_output=True, text=True)
    elapsed = time.perf_counter() - start
    log_path.write_text(result.stdout, encoding="utf-8")
    if result.stderr:
        stderr_path = log_path.with_suffix(log_path.suffix + ".stderr")
        stderr_path.write_text(result.stderr, encoding="utf-8")
    if result.returncode != 0:
        if result.stdout:
            print(result.stdout, file=sys.stdout, end="")
        if result.stderr:
            print(result.stderr, file=sys.stderr, end="")
        raise subprocess.CalledProcessError(result.returncode, command)
    return result.stdout, elapsed


def as_int(values: dict[str, str], key: str) -> int:
    return int(values[key])


def as_float(values: dict[str, str], key: str) -> float:
    return float(values[key])


def compare_rows(
    afterburner_rows: list[dict[str, str]],
    summary_rows: list[dict[str, str]],
) -> tuple[dict[str, object], list[dict[str, object]], list[str]]:
    mismatches: list[str] = []
    summary_by_event = {as_int(row, "event_index"): row for row in summary_rows}
    shape_status_counts: dict[str, int] = {}
    event_summaries: list[dict[str, object]] = []
    events = len(afterburner_rows)
    skipped = 0
    requested = 0
    retained = 0
    photons_added = 0
    summary_photons = 0
    summary_weight_sum = 0.0
    max_selected_outgoing_charged = 0
    max_requested_photons = 0
    max_retained_photons = 0
    max_trial_bound_ratio = 0.0
    event_seconds_sum = 0.0
    max_event_seconds = 0.0
    event_seconds_count = 0

    if len(summary_by_event) != len(summary_rows):
        mismatches.append("duplicate event_index in output-summary diagnostic")

    for row in afterburner_rows:
        event_index = as_int(row, "event_index")
        summary = summary_by_event.get(event_index)
        if summary is None:
            mismatches.append(f"missing output-summary row for event {event_index}")
            continue

        status = row["photon_status"]
        selected_outgoing = as_int(row, "selected_outgoing_charged")
        selected_incoming = as_int(row, "selected_incoming_charged")
        ignored_neutral = as_int(row, "ignored_neutral_candidates")
        rejected_charged = as_int(row, "rejected_charged_candidates")
        if status.startswith("skipped_"):
            skipped += 1
        requested_i = as_int(row, "requested_photons")
        retained_i = as_int(row, "retained_photons")
        added_i = as_int(row, "photons_added")
        summary_i = as_int(summary, "lowmc_photons")

        requested += requested_i
        retained += retained_i
        photons_added += added_i
        summary_photons += summary_i
        summary_weight_sum += as_float(summary, "photon_weight_sum")
        max_selected_outgoing_charged = max(max_selected_outgoing_charged, selected_outgoing)
        max_requested_photons = max(max_requested_photons, requested_i)
        max_retained_photons = max(max_retained_photons, retained_i)

        shape_status = row.get("shape_status", "none")
        shape_status_counts[shape_status] = shape_status_counts.get(shape_status, 0) + 1
        event_trial_bound_ratio = 0.0
        if "max_trial_bound_ratio" in row:
            event_trial_bound_ratio = as_float(row, "max_trial_bound_ratio")
            max_trial_bound_ratio = max(max_trial_bound_ratio, event_trial_bound_ratio)
        event_seconds: object = ""
        if "event_seconds" in row:
            event_seconds_value = as_float(row, "event_seconds")
            event_seconds = event_seconds_value
            event_seconds_sum += event_seconds_value
            max_event_seconds = max(max_event_seconds, event_seconds_value)
            event_seconds_count += 1

        event_summaries.append(
            {
                "event_index": event_index,
                "photon_status": status,
                "selected_incoming_charged": selected_incoming,
                "selected_outgoing_charged": selected_outgoing,
                "ignored_neutral_candidates": ignored_neutral,
                "rejected_charged_candidates": rejected_charged,
                "signed_charge_sum": row["signed_charge_sum"],
                "ordinary_charge_sum": row["ordinary_charge_sum"],
                "lambda_event": row["lambda_event"],
                "requested_photons": requested_i,
                "retained_photons": retained_i,
                "photons_added": added_i,
                "summary_lowmc_photons": summary_i,
                "summary_photon_weight_sum": as_float(summary, "photon_weight_sum"),
                "shape_status": shape_status,
                "max_trial_bound_ratio": event_trial_bound_ratio,
                "event_seconds": event_seconds,
                "event_metadata": summary.get("event_metadata", "0"),
                "photon_weight_attributes": summary.get("photon_weight_attributes", "0"),
            }
        )

        if retained_i != requested_i:
            mismatches.append(f"event {event_index}: retained != requested")
        if added_i != requested_i:
            mismatches.append(f"event {event_index}: photons_added != requested")
        if summary_i != added_i:
            mismatches.append(f"event {event_index}: summary photons != afterburner photons")
        if added_i > 0 and summary.get("event_metadata") != "1":
            mismatches.append(f"event {event_index}: missing LowMC event metadata")
        if added_i > 0 and as_int(summary, "photon_weight_attributes") != added_i:
            mismatches.append(f"event {event_index}: photon weight attribute count mismatch")

    aggregate = {
        "events": events,
        "skipped_events": skipped,
        "requested_photons": requested,
        "retained_photons": retained,
        "photons_added": photons_added,
        "summary_lowmc_photons": summary_photons,
        "summary_photon_weight_sum": summary_weight_sum,
        "max_selected_outgoing_charged": max_selected_outgoing_charged,
        "max_requested_photons": max_requested_photons,
        "max_retained_photons": max_retained_photons,
        "event_seconds_sum": event_seconds_sum,
        "max_event_seconds": max_event_seconds,
        "event_seconds_count": event_seconds_count,
        "shape_status_counts": ",".join(
            f"{key}:{shape_status_counts[key]}" for key in sorted(shape_status_counts)
        ),
        "max_trial_bound_ratio": max_trial_bound_ratio,
        "mismatch_count": len(mismatches),
    }
    return aggregate, event_summaries, mismatches


def write_summary_csv(path: Path, row: dict[str, object]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row))
        writer.writeheader()
        writer.writerow(row)


def write_event_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    if args.events <= 0:
        raise ValueError("--events must be positive")
    if args.enhancement_factor <= 0.0 or not math.isfinite(args.enhancement_factor):
        raise ValueError("--enhancement-factor must be positive and finite")
    if args.ecm <= 0.0:
        raise ValueError("--ecm must be positive")

    script_dir = Path(__file__).resolve().parent
    sample_helper = script_dir / "generate_pythia_hepmc3_sample.py"
    afterburner = args.tool_dir / "lowmc_hepmc3_afterburner_diagnostic"
    output_summary = args.tool_dir / "lowmc_hepmc3_output_summary_diagnostic"
    require_file(sample_helper, "Pythia-to-HepMC3 sample helper")
    require_file(afterburner, "HepMC3 afterburner diagnostic")
    require_file(output_summary, "HepMC3 output-summary diagnostic")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    generated_input = args.output_dir / f"pythia_softqcd_{args.events}.hepmc3"
    afterburner_output = args.output_dir / f"pythia_softqcd_{args.events}_lowmc.hepmc3"
    generator_work_dir = args.output_dir / "generator"

    generator_command = [
        sys.executable,
        str(sample_helper),
        "--pythia-prefix",
        str(args.pythia_prefix),
        "--hepmc3-prefix",
        str(args.hepmc3_prefix),
        "--output",
        str(generated_input),
        "--work-dir",
        str(generator_work_dir),
        "--events",
        str(args.events),
        "--seed",
        str(args.pythia_seed),
        "--ecm",
        f"{args.ecm:.17g}",
        "--cxx",
        args.cxx,
    ]
    _, generator_seconds = run_logged(generator_command, args.output_dir / "generate.log")

    afterburner_command = [
        str(afterburner),
        "--max-events",
        str(args.events),
        "--enhancement-factor",
        f"{args.enhancement_factor:.17g}",
        "--multiplicity-seed",
        str(args.multiplicity_seed),
        "--shape-seed",
        str(args.shape_seed),
        "--kt-min",
        f"{args.kt_min:.17g}",
        "--kt-max",
        f"{args.kt_max:.17g}",
        "--eta-min",
        f"{args.eta_min:.17g}",
        "--eta-max",
        f"{args.eta_max:.17g}",
        "--phi-min",
        f"{args.phi_min:.17g}",
        "--phi-max",
        f"{args.phi_max:.17g}",
        "--atomic-output",
        str(generated_input),
        str(afterburner_output),
    ]
    afterburner_log, afterburner_seconds = run_logged(
        afterburner_command,
        args.output_dir / "afterburner.log",
    )

    summary_command = [
        str(output_summary),
        "--max-events",
        str(args.events),
        str(afterburner_output),
    ]
    summary_log, summary_seconds = run_logged(summary_command, args.output_dir / "summary.log")

    afterburner_rows = rows_with_key(afterburner_log, "photon_status")
    summary_rows = rows_with_key(summary_log, "lowmc_photons")
    if not afterburner_rows:
        raise ValueError("afterburner log contains no event rows")
    if not summary_rows:
        raise ValueError("summary log contains no event rows")
    aggregate, event_summaries, mismatches = compare_rows(afterburner_rows, summary_rows)
    estimated_seconds_per_event = afterburner_seconds / float(aggregate["events"])
    aggregate.update(
        {
            "input": str(generated_input),
            "output": str(afterburner_output),
            "generator_seconds": generator_seconds,
            "afterburner_seconds": afterburner_seconds,
            "afterburner_seconds_per_event": estimated_seconds_per_event,
            "summary_seconds": summary_seconds,
            "enhancement_factor": args.enhancement_factor,
            "pythia_seed": args.pythia_seed,
            "multiplicity_seed": args.multiplicity_seed,
            "shape_seed": args.shape_seed,
        }
    )
    for row in event_summaries:
        row["estimated_afterburner_seconds_per_event"] = estimated_seconds_per_event
    write_summary_csv(args.output_dir / "validation_summary.csv", aggregate)
    write_event_csv(args.output_dir / "validation_events.csv", event_summaries)

    print(
        "validation "
        f"events={aggregate['events']} skipped_events={aggregate['skipped_events']} "
        f"photons_added={aggregate['photons_added']} "
        f"summary_lowmc_photons={aggregate['summary_lowmc_photons']} "
        f"summary_photon_weight_sum={aggregate['summary_photon_weight_sum']:.12g} "
        f"max_selected_outgoing_charged={aggregate['max_selected_outgoing_charged']} "
        f"max_requested_photons={aggregate['max_requested_photons']} "
        f"max_event_seconds={aggregate['max_event_seconds']:.6g} "
        f"shape_status_counts={aggregate['shape_status_counts']} "
        f"max_trial_bound_ratio={aggregate['max_trial_bound_ratio']:.12g} "
        f"mismatch_count={aggregate['mismatch_count']}"
    )
    print(
        "timing "
        f"generator_seconds={generator_seconds:.3f} "
        f"afterburner_seconds={afterburner_seconds:.3f} "
        f"summary_seconds={summary_seconds:.3f}"
    )
    print(
        "artifacts "
        f"input={generated_input} output={afterburner_output} "
        f"afterburner_log={args.output_dir / 'afterburner.log'} "
        f"summary_log={args.output_dir / 'summary.log'} "
        f"summary_csv={args.output_dir / 'validation_summary.csv'} "
        f"event_csv={args.output_dir / 'validation_events.csv'}"
    )
    for mismatch in mismatches:
        print(f"mismatch {mismatch}")
    if mismatches and not args.no_fail_on_mismatch:
        raise AssertionError(f"validation found {len(mismatches)} mismatch(es)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
