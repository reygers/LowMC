#!/usr/bin/env python3
"""Summarize paired Pythia leg-centred versus grid residual-tail runs."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
import statistics


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize Pythia analytic-leg versus tuned-grid tail CSV files."
    )
    parser.add_argument("inputs", type=Path, nargs="+")
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--paired", type=Path, required=True)
    parser.add_argument("--plot", type=Path)
    return parser.parse_args()


def quantile(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    if not ordered:
        return 0.0
    position = probability * (len(ordered) - 1)
    lower = int(position)
    upper = min(lower + 1, len(ordered) - 1)
    fraction = position - lower
    return ordered[lower] + fraction * (ordered[upper] - ordered[lower])


def load_rows(paths: list[Path]) -> list[dict[str, object]]:
    numeric_float = {
        "max_outgoing_beta",
        "max_outgoing_abs_eta",
        "mean_weight",
        "standard_error",
        "effective_sample_fraction",
        "max_weight",
        "q99_weight",
        "q999_weight",
        "q9999_weight",
        "inside_fraction",
        "build_milliseconds",
        "sample_milliseconds",
    }
    numeric_int = {
        "seed",
        "event_index",
        "selected_charged",
        "selected_incoming",
        "selected_outgoing",
        "samples",
        "proposal_elements",
    }
    rows: list[dict[str, object]] = []
    for path in paths:
        with path.open(newline="") as stream:
            for raw in csv.DictReader(stream):
                row: dict[str, object] = dict(raw)
                for field in numeric_float:
                    row[field] = float(raw[field])
                for field in numeric_int:
                    row[field] = int(raw[field])
                rows.append(row)
    return rows


def make_pairs(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, int, int], dict[str, dict[str, object]]] = {}
    for row in rows:
        key = (str(row["mode"]), int(row["seed"]), int(row["event_index"]))
        grouped.setdefault(key, {})[str(row["proposal"])] = row

    pairs = []
    for (mode, seed, event_index), proposals in sorted(grouped.items()):
        if set(proposals) != {"leg", "grid"}:
            raise ValueError(f"missing paired proposal for {mode}, seed={seed}, event={event_index}")
        leg = proposals["leg"]
        grid = proposals["grid"]
        leg_mean = float(leg["mean_weight"])
        grid_mean = float(grid["mean_weight"])
        average = 0.5 * (leg_mean + grid_mean)
        uncertainty = math.hypot(
            float(leg["standard_error"]), float(grid["standard_error"])
        )
        pairs.append(
            {
                "mode": mode,
                "seed": seed,
                "event_index": event_index,
                "selected_charged": int(leg["selected_charged"]),
                "selected_outgoing": int(leg["selected_outgoing"]),
                "max_outgoing_beta": float(leg["max_outgoing_beta"]),
                "max_outgoing_abs_eta": float(leg["max_outgoing_abs_eta"]),
                "leg_mean_weight": leg_mean,
                "grid_mean_weight": grid_mean,
                "relative_mean_difference": (
                    (leg_mean - grid_mean) / average if average != 0.0 else 0.0
                ),
                "mean_difference_pull": (
                    (leg_mean - grid_mean) / uncertainty if uncertainty > 0.0 else 0.0
                ),
                "leg_ess_fraction": float(leg["effective_sample_fraction"]),
                "grid_ess_fraction": float(grid["effective_sample_fraction"]),
                "leg_max_over_mean": float(leg["max_weight"]) / leg_mean,
                "grid_max_over_mean": float(grid["max_weight"]) / grid_mean,
                "leg_q9999_over_mean": float(leg["q9999_weight"]) / leg_mean,
                "grid_q9999_over_mean": float(grid["q9999_weight"]) / grid_mean,
                "leg_inside_fraction": float(leg["inside_fraction"]),
                "leg_build_milliseconds": float(leg["build_milliseconds"]),
                "grid_build_milliseconds": float(grid["build_milliseconds"]),
                "leg_sample_milliseconds": float(leg["sample_milliseconds"]),
                "grid_sample_milliseconds": float(grid["sample_milliseconds"]),
            }
        )
    return pairs


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def make_summary(rows: list[dict[str, object]], pairs: list[dict[str, object]]) -> list[dict[str, object]]:
    output = []
    modes = sorted({str(row["mode"]) for row in rows})
    for mode in modes:
        mode_pairs = [row for row in pairs if row["mode"] == mode]
        leg_better = sum(
            float(row["leg_ess_fraction"]) > float(row["grid_ess_fraction"])
            for row in mode_pairs
        ) / len(mode_pairs)
        pull_failures = sum(abs(float(row["mean_difference_pull"])) > 3.0 for row in mode_pairs)
        relative_differences = [float(row["relative_mean_difference"]) for row in mode_pairs]
        for proposal in ("leg", "grid"):
            selected = [
                row for row in rows
                if row["mode"] == mode and row["proposal"] == proposal
            ]
            ess = [float(row["effective_sample_fraction"]) for row in selected]
            max_over_mean = [
                float(row["max_weight"]) / float(row["mean_weight"]) for row in selected
            ]
            q9999_over_mean = [
                float(row["q9999_weight"]) / float(row["mean_weight"]) for row in selected
            ]
            effective_samples_per_second = [
                float(row["samples"])
                * float(row["effective_sample_fraction"])
                * 1000.0
                / (
                    float(row["build_milliseconds"])
                    + float(row["sample_milliseconds"])
                )
                for row in selected
            ]
            output.append(
                {
                    "mode": mode,
                    "proposal": proposal,
                    "events": len(selected),
                    "median_ess_fraction": statistics.median(ess),
                    "q10_ess_fraction": quantile(ess, 0.10),
                    "minimum_ess_fraction": min(ess),
                    "median_max_over_mean": statistics.median(max_over_mean),
                    "maximum_max_over_mean": max(max_over_mean),
                    "median_q9999_over_mean": statistics.median(q9999_over_mean),
                    "mean_inside_fraction": statistics.mean(
                        float(row["inside_fraction"]) for row in selected
                    ),
                    "mean_build_milliseconds": statistics.mean(
                        float(row["build_milliseconds"]) for row in selected
                    ),
                    "mean_sample_milliseconds": statistics.mean(
                        float(row["sample_milliseconds"]) for row in selected
                    ),
                    "median_effective_samples_per_second": statistics.median(
                        effective_samples_per_second
                    ),
                    "fraction_events_leg_has_higher_ess": leg_better,
                    "paired_abs_pull_gt_3": pull_failures,
                    "paired_median_relative_mean_difference": statistics.median(
                        relative_differences
                    ),
                }
            )
    return output


def make_plot(path: Path, rows: list[dict[str, object]], pairs: list[dict[str, object]]) -> None:
    import matplotlib.pyplot as plt

    colors = {"leg": "#d97706", "grid": "#2563eb"}
    markers = {"ee-hadrons": "o", "softqcd": "s"}
    figure, axes = plt.subplots(1, 3, figsize=(13.5, 4.2))
    for mode in markers:
        for proposal in colors:
            selected = [
                row for row in rows
                if row["mode"] == mode and row["proposal"] == proposal
            ]
            axes[0].scatter(
                [row["selected_charged"] for row in selected],
                [row["effective_sample_fraction"] for row in selected],
                s=20,
                alpha=0.7,
                color=colors[proposal],
                marker=markers[mode],
                label=f"{mode}, {proposal}",
            )
            axes[1].scatter(
                [row["selected_charged"] for row in selected],
                [float(row["max_weight"]) / float(row["mean_weight"]) for row in selected],
                s=20,
                alpha=0.7,
                color=colors[proposal],
                marker=markers[mode],
            )
    for mode, marker in markers.items():
        selected_pairs = [row for row in pairs if row["mode"] == mode]
        axes[2].scatter(
            [row["selected_charged"] for row in selected_pairs],
            [row["relative_mean_difference"] for row in selected_pairs],
            s=22,
            alpha=0.75,
            marker=marker,
            label=mode,
        )
    axes[0].set(xlabel="selected charged particles", ylabel="ESS / N", ylim=(0, 1.02))
    axes[1].set(xlabel="selected charged particles", ylabel="maximum weight / mean")
    axes[1].set_yscale("log")
    axes[2].axhline(0.0, color="black", linewidth=0.8)
    axes[2].set(xlabel="selected charged particles", ylabel="(leg - grid) / paired mean")
    axes[0].legend(fontsize=8)
    axes[2].legend(fontsize=8)
    figure.suptitle("Pythia residual-tail gate for analytic leg-centred proposal")
    figure.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=180)
    plt.close(figure)


def main() -> int:
    args = parse_args()
    rows = load_rows(args.inputs)
    if not rows:
        raise ValueError("no diagnostic rows found")
    pairs = make_pairs(rows)
    summary = make_summary(rows, pairs)
    write_csv(args.summary, summary)
    write_csv(args.paired, pairs)
    if args.plot:
        make_plot(args.plot, rows, pairs)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
