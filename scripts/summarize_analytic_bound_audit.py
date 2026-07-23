#!/usr/bin/env python3
"""Summarize held-out Pythia analytic Low-bound audit CSV files."""

import argparse
import csv
import math
import statistics
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Aggregate Pythia analytic-bound audit rows, apply the "
            "predeclared safety/tightness gates, and optionally plot them."
        )
    )
    parser.add_argument("inputs", nargs="+", type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--plot", type=Path)
    return parser.parse_args()


def percentile(values, fraction):
    ordered = sorted(values)
    if not ordered:
        return math.nan
    position = fraction * (len(ordered) - 1)
    lower = int(math.floor(position))
    upper = int(math.ceil(position))
    if lower == upper:
        return ordered[lower]
    weight = position - lower
    return ordered[lower] * (1.0 - weight) + ordered[upper] * weight


def read_rows(paths):
    rows = []
    for path in paths:
        with path.open(newline="", encoding="utf-8") as stream:
            rows.extend(csv.DictReader(stream))
    if not rows:
        raise ValueError("no analytic-bound audit rows found")
    return rows


def summarize_mode(mode, rows):
    selected = [row for row in rows if row["mode"] == mode]
    if not selected:
        raise ValueError(f"no rows found for mode {mode}")
    envelope_ratios = [
        float(row["analytic_cell_envelope_over_current"])
        for row in selected
    ]
    analytic_times = [
        float(row["analytic_seconds"]) for row in selected
    ]
    base_probe_times = [
        float(row["base_probe_seconds"]) for row in selected
    ]
    current_acceptance = [
        float(row["current_expected_acceptance"])
        for row in selected
    ]
    analytic_acceptance = [
        float(row["analytic_expected_acceptance"])
        for row in selected
    ]
    global_ratios = [
        float(row["analytic_global_weight_bound"])
        / float(row["current_envelope_integral"])
        for row in selected
    ]

    invalid_particles = sum(
        int(row["analytic_invalid_particles"]) for row in selected
    )
    fallback_cells = sum(
        int(row["analytic_fallback_cells"]) for row in selected
    )
    observed_violations = sum(
        int(row["analytic_observed_violations"]) for row in selected
    )
    trial_violations = sum(
        int(row["analytic_trial_violations"]) for row in selected
    )
    unavailable = sum(
        int(row["analytic_bound_available"]) == 0 for row in selected
    )
    max_observed_ratio = max(
        float(row["analytic_max_observed_over_bound"])
        for row in selected
    )
    max_trial_ratio = max(
        float(row["analytic_max_trial_over_bound"])
        for row in selected
    )
    safety_pass = (
        invalid_particles == 0
        and fallback_cells == 0
        and observed_violations == 0
        and trial_violations == 0
        and unavailable == 0
        and max_observed_ratio <= 1.0
        and max_trial_ratio <= 1.0
    )
    median_envelope_ratio = statistics.median(envelope_ratios)
    maximum_envelope_ratio = max(envelope_ratios)
    tightness_pass = (
        median_envelope_ratio <= 2.0
        and maximum_envelope_ratio <= 10.0
    )
    total_analytic_seconds = sum(analytic_times)
    total_base_probe_seconds = sum(base_probe_times)
    return {
        "mode": mode,
        "events": len(selected),
        "minimum_selected_charged": min(
            int(row["selected_charged"]) for row in selected
        ),
        "median_selected_charged": statistics.median(
            int(row["selected_charged"]) for row in selected
        ),
        "maximum_selected_charged": max(
            int(row["selected_charged"]) for row in selected
        ),
        "invalid_particles": invalid_particles,
        "fallback_cells": fallback_cells,
        "observed_bound_violations": observed_violations,
        "trial_bound_violations": trial_violations,
        "unavailable_events": unavailable,
        "maximum_observed_over_bound": max_observed_ratio,
        "maximum_trial_over_bound": max_trial_ratio,
        "safety_gate_pass": int(safety_pass),
        "median_cell_envelope_over_current": median_envelope_ratio,
        "p90_cell_envelope_over_current": percentile(
            envelope_ratios, 0.9
        ),
        "maximum_cell_envelope_over_current": maximum_envelope_ratio,
        "tightness_gate_pass": int(tightness_pass),
        "median_current_expected_acceptance": statistics.median(
            current_acceptance
        ),
        "median_analytic_expected_acceptance": statistics.median(
            analytic_acceptance
        ),
        "total_analytic_seconds": total_analytic_seconds,
        "total_base_probe_seconds": total_base_probe_seconds,
        "analytic_over_base_probe_time": (
            total_analytic_seconds / total_base_probe_seconds
        ),
        "median_global_bound_over_current_envelope": statistics.median(
            global_ratios
        ),
    }


def write_summary(path, summaries):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=list(summaries[0])
        )
        writer.writeheader()
        writer.writerows(summaries)


def make_plot(path, rows):
    import matplotlib.pyplot as plt

    modes = ["ee-hadrons", "softqcd"]
    labels = ["Hadronic Z", "SoftQCD"]
    colors = ["#0072B2", "#D55E00"]
    ratios = {
        mode: [
            float(row["analytic_cell_envelope_over_current"])
            for row in rows
            if row["mode"] == mode
        ]
        for mode in modes
    }
    time_ratios = {
        mode: [
            float(row["analytic_seconds"])
            / float(row["base_probe_seconds"])
            for row in rows
            if row["mode"] == mode
        ]
        for mode in modes
    }

    path.parent.mkdir(parents=True, exist_ok=True)
    figure, (bound_axis, time_axis) = plt.subplots(
        1, 2, figsize=(8.5, 4.2)
    )
    positions = [1, 2]
    bound_plot = bound_axis.boxplot(
        [ratios[mode] for mode in modes],
        positions=positions,
        patch_artist=True,
        tick_labels=labels,
    )
    for patch, color in zip(bound_plot["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    bound_axis.axhline(
        2.0,
        color="black",
        linestyle="--",
        linewidth=1,
        label="median gate",
    )
    bound_axis.axhline(
        10.0,
        color="black",
        linestyle=":",
        linewidth=1,
        label="event gate",
    )
    bound_axis.set(
        ylabel="analytic / current envelope integral",
        yscale="log",
    )
    bound_axis.grid(axis="y", alpha=0.22)
    bound_axis.legend(frameon=False, fontsize=8)

    time_plot = time_axis.boxplot(
        [time_ratios[mode] for mode in modes],
        positions=positions,
        patch_artist=True,
        tick_labels=labels,
    )
    for patch, color in zip(time_plot["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    time_axis.axhline(1.0, color="black", linewidth=1)
    time_axis.set(
        ylabel="analytic construction / base-probe time",
        yscale="log",
    )
    time_axis.grid(axis="y", alpha=0.22)
    figure.suptitle(
        "Held-out analytic Low-bound audit",
        fontsize=14,
    )
    figure.tight_layout()
    figure.savefig(path, dpi=180)
    plt.close(figure)


def main():
    args = parse_args()
    rows = read_rows(args.inputs)
    summaries = [
        summarize_mode(mode, rows)
        for mode in ("ee-hadrons", "softqcd")
    ]
    write_summary(args.summary, summaries)
    if args.plot:
        make_plot(args.plot, rows)
    for summary in summaries:
        print(
            f"mode={summary['mode']} "
            f"events={summary['events']} "
            f"safety_gate_pass={summary['safety_gate_pass']} "
            f"invalid_particles={summary['invalid_particles']} "
            f"fallback_cells={summary['fallback_cells']} "
            f"observed_violations="
            f"{summary['observed_bound_violations']} "
            f"trial_violations={summary['trial_bound_violations']} "
            f"max_observed_over_bound="
            f"{summary['maximum_observed_over_bound']:.6g} "
            f"max_trial_over_bound="
            f"{summary['maximum_trial_over_bound']:.6g} "
            f"median_envelope_ratio="
            f"{summary['median_cell_envelope_over_current']:.6g} "
            f"max_envelope_ratio="
            f"{summary['maximum_cell_envelope_over_current']:.6g} "
            f"tightness_gate_pass="
            f"{summary['tightness_gate_pass']} "
            f"analytic_over_base_probe_time="
            f"{summary['analytic_over_base_probe_time']:.6g}"
        )


if __name__ == "__main__":
    main()
