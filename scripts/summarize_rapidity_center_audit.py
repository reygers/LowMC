#!/usr/bin/env python3
"""Summarize paired pseudorapidity/rapidity shape-envelope audits."""

import argparse
import csv
import math
import statistics
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Pair shape-unweighting audit rows that differ only in the "
            "charged-particle eta-center convention."
        )
    )
    parser.add_argument("inputs", nargs="+", type=Path)
    parser.add_argument("--pairs", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--plot", type=Path)
    return parser.parse_args()


def read_rows(paths):
    rows = []
    for path in paths:
        with path.open(newline="", encoding="utf-8") as stream:
            for row in csv.DictReader(stream):
                if "eta_center" not in row:
                    raise ValueError(
                        f"{path} has no eta_center column; rerun the audit"
                    )
                rows.append(row)
    return rows


def finite_ratio(numerator, denominator):
    if denominator == 0.0:
        return math.nan
    return numerator / denominator


def pair_rows(rows):
    grouped = {}
    for row in rows:
        key = (
            row["mode"],
            row["seed"],
            row["event_index"],
            row["shape_seed"],
        )
        center = row["eta_center"]
        if center in grouped.setdefault(key, {}):
            raise ValueError(f"duplicate {center} row for {key}")
        grouped[key][center] = row

    pairs = []
    for key, centers in sorted(grouped.items()):
        if set(centers) != {"pseudorapidity", "rapidity"}:
            raise ValueError(f"incomplete eta-center pair for {key}")
        pseudo = centers["pseudorapidity"]
        rapid = centers["rapidity"]
        if pseudo["selected_charged"] != rapid["selected_charged"]:
            raise ValueError(f"particle selection differs for {key}")

        def number(row, field):
            return float(row[field])

        envelope_ratio = finite_ratio(
            number(rapid, "current_envelope_integral"),
            number(pseudo, "current_envelope_integral"),
        )
        pairs.append(
            {
                "mode": key[0],
                "seed": key[1],
                "event_index": key[2],
                "shape_seed": key[3],
                "selected_charged": int(pseudo["selected_charged"]),
                "median_abs_eta_center_difference": number(
                    pseudo, "median_abs_eta_center_difference"
                ),
                "max_abs_eta_center_difference": number(
                    pseudo, "max_abs_eta_center_difference"
                ),
                "pseudo_grid_cells": int(pseudo["grid_cells"]),
                "rapid_grid_cells": int(rapid["grid_cells"]),
                "grid_cell_difference": (
                    int(rapid["grid_cells"]) -
                    int(pseudo["grid_cells"])
                ),
                "pseudo_probe_evaluations": int(
                    pseudo["probe_evaluations"]
                ),
                "rapid_probe_evaluations": int(
                    rapid["probe_evaluations"]
                ),
                "probe_evaluation_difference": (
                    int(rapid["probe_evaluations"]) -
                    int(pseudo["probe_evaluations"])
                ),
                "promoted_cell_difference": (
                    int(rapid["promoted_cells"]) -
                    int(pseudo["promoted_cells"])
                ),
                "envelope_integral_ratio": envelope_ratio,
                "implied_acceptance_ratio": finite_ratio(
                    1.0, envelope_ratio
                ),
                "expected_acceptance_ratio": finite_ratio(
                    number(rapid, "current_expected_acceptance"),
                    number(pseudo, "current_expected_acceptance"),
                ),
                "observed_acceptance_difference": (
                    number(rapid, "current_observed_acceptance") -
                    number(pseudo, "current_observed_acceptance")
                ),
                "pseudo_envelope_failures": int(
                    pseudo["current_envelope_failures"]
                ),
                "rapid_envelope_failures": int(
                    rapid["current_envelope_failures"]
                ),
                "pseudo_max_trial_bound_ratio": number(
                    pseudo, "max_trial_bound_ratio"
                ),
                "rapid_max_trial_bound_ratio": number(
                    rapid, "max_trial_bound_ratio"
                ),
                "runtime_ratio": finite_ratio(
                    number(rapid, "total_measured_seconds"),
                    number(pseudo, "total_measured_seconds"),
                ),
            }
        )
    if not pairs:
        raise ValueError("no complete eta-center pairs found")
    return pairs


def median(rows, field):
    values = [
        float(row[field])
        for row in rows
        if math.isfinite(float(row[field]))
    ]
    return statistics.median(values) if values else math.nan


def summarize(pairs):
    summaries = []
    for mode in sorted({row["mode"] for row in pairs}):
        selected = [row for row in pairs if row["mode"] == mode]
        summaries.append(
            {
                "mode": mode,
                "events": len(selected),
                "median_selected_charged": statistics.median(
                    row["selected_charged"] for row in selected
                ),
                "median_abs_eta_center_difference": median(
                    selected, "median_abs_eta_center_difference"
                ),
                "maximum_abs_eta_center_difference": max(
                    row["max_abs_eta_center_difference"]
                    for row in selected
                ),
                "median_grid_cell_difference": median(
                    selected, "grid_cell_difference"
                ),
                "median_probe_evaluation_difference": median(
                    selected, "probe_evaluation_difference"
                ),
                "median_promoted_cell_difference": median(
                    selected, "promoted_cell_difference"
                ),
                "median_envelope_integral_ratio": median(
                    selected, "envelope_integral_ratio"
                ),
                "maximum_envelope_integral_ratio": max(
                    row["envelope_integral_ratio"] for row in selected
                ),
                "median_expected_acceptance_ratio": median(
                    selected, "expected_acceptance_ratio"
                ),
                "median_implied_acceptance_ratio": median(
                    selected, "implied_acceptance_ratio"
                ),
                "median_observed_acceptance_difference": median(
                    selected, "observed_acceptance_difference"
                ),
                "pseudo_envelope_failures": sum(
                    row["pseudo_envelope_failures"] for row in selected
                ),
                "rapid_envelope_failures": sum(
                    row["rapid_envelope_failures"] for row in selected
                ),
                "maximum_pseudo_trial_bound_ratio": max(
                    row["pseudo_max_trial_bound_ratio"]
                    for row in selected
                ),
                "maximum_rapid_trial_bound_ratio": max(
                    row["rapid_max_trial_bound_ratio"]
                    for row in selected
                ),
                "median_runtime_ratio": median(
                    selected, "runtime_ratio"
                ),
            }
        )
    return summaries


def write_csv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def make_plot(path, pairs):
    import matplotlib.pyplot as plt

    modes = sorted({row["mode"] for row in pairs})
    labels = [
        "Hadronic Z" if mode == "ee-hadrons" else "SoftQCD"
        for mode in modes
    ]
    colors = ["#0072B2", "#D55E00"]
    path.parent.mkdir(parents=True, exist_ok=True)
    figure, (envelope_axis, center_axis) = plt.subplots(
        1, 2, figsize=(8.5, 4.2)
    )
    plot = envelope_axis.boxplot(
        [
            [
                row["envelope_integral_ratio"]
                for row in pairs
                if row["mode"] == mode
            ]
            for mode in modes
        ],
        patch_artist=True,
        tick_labels=labels,
    )
    for patch, color in zip(plot["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    envelope_axis.axhline(1.0, color="black", linewidth=1)
    envelope_axis.set_ylabel(
        "true-rapidity / pseudorapidity envelope"
    )
    envelope_axis.grid(axis="y", alpha=0.22)

    for mode, label, color in zip(modes, labels, colors):
        selected = [row for row in pairs if row["mode"] == mode]
        center_axis.scatter(
            [
                row["max_abs_eta_center_difference"]
                for row in selected
            ],
            [row["envelope_integral_ratio"] for row in selected],
            label=label,
            color=color,
            alpha=0.75,
        )
    center_axis.axhline(1.0, color="black", linewidth=1)
    center_axis.set(
        xlabel="maximum |pseudorapidity - rapidity|",
        ylabel="envelope ratio",
    )
    center_axis.grid(alpha=0.22)
    center_axis.legend(frameon=False, fontsize=8)
    figure.suptitle("Charged-particle eta-center A/B audit")
    figure.tight_layout()
    figure.savefig(path, dpi=180)
    plt.close(figure)


def main():
    args = parse_args()
    pairs = pair_rows(read_rows(args.inputs))
    summaries = summarize(pairs)
    write_csv(args.pairs, pairs)
    write_csv(args.summary, summaries)
    if args.plot:
        make_plot(args.plot, pairs)
    for row in summaries:
        print(
            f"mode={row['mode']} events={row['events']} "
            f"median_center_difference="
            f"{row['median_abs_eta_center_difference']:.6g} "
            f"max_center_difference="
            f"{row['maximum_abs_eta_center_difference']:.6g} "
            f"median_envelope_ratio="
            f"{row['median_envelope_integral_ratio']:.6g} "
            f"median_acceptance_ratio="
            f"{row['median_expected_acceptance_ratio']:.6g} "
            f"pseudo_failures={row['pseudo_envelope_failures']} "
            f"rapid_failures={row['rapid_envelope_failures']} "
            f"median_runtime_ratio={row['median_runtime_ratio']:.6g}"
        )


if __name__ == "__main__":
    main()
