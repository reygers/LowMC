#!/usr/bin/env python3
"""Summarize independent calibration/trial tests of the analytic leg envelope."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
import statistics


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize Pythia analytic leg-envelope diagnostic CSV files."
    )
    parser.add_argument("inputs", type=Path, nargs="+")
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--plot", type=Path)
    return parser.parse_args()


def quantile(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    position = probability * (len(ordered) - 1)
    lower = int(position)
    upper = min(lower + 1, len(ordered) - 1)
    fraction = position - lower
    return ordered[lower] + fraction * (ordered[upper] - ordered[lower])


def load_rows(paths: list[Path]) -> list[dict[str, object]]:
    float_fields = {
        "safety_factor",
        "max_trial_bound_ratio",
        "expected_acceptance",
        "observed_acceptance",
        "proposal_build_milliseconds",
        "calibration_milliseconds",
        "trial_milliseconds",
    }
    int_fields = {
        "seed",
        "event_index",
        "calibration_samples",
        "trial_samples",
        "overshoots",
        "envelope_ok",
    }
    rows: list[dict[str, object]] = []
    for path in paths:
        with path.open(newline="") as stream:
            for raw in csv.DictReader(stream):
                row: dict[str, object] = dict(raw)
                for field in float_fields:
                    row[field] = float(raw[field])
                for field in int_fields:
                    row[field] = int(raw[field])
                rows.append(row)
    return rows


def make_summary(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    output = []
    keys = sorted(
        {(str(row["mode"]), float(row["safety_factor"])) for row in rows}
    )
    for mode, safety_factor in keys:
        selected = [
            row
            for row in rows
            if row["mode"] == mode and row["safety_factor"] == safety_factor
        ]
        ratios = [float(row["max_trial_bound_ratio"]) for row in selected]
        expected = [float(row["expected_acceptance"]) for row in selected]
        observed = [float(row["observed_acceptance"]) for row in selected]
        output.append(
            {
                "mode": mode,
                "safety_factor": safety_factor,
                "events": len(selected),
                "failed_events": sum(int(row["envelope_ok"]) == 0 for row in selected),
                "total_trial_samples": sum(int(row["trial_samples"]) for row in selected),
                "total_overshoots": sum(int(row["overshoots"]) for row in selected),
                "maximum_trial_bound_ratio": max(ratios),
                "q90_trial_bound_ratio": quantile(ratios, 0.90),
                "median_trial_bound_ratio": statistics.median(ratios),
                "mean_expected_acceptance": statistics.mean(expected),
                "minimum_expected_acceptance": min(expected),
                "mean_observed_acceptance": statistics.mean(observed),
                "mean_proposal_build_milliseconds": statistics.mean(
                    float(row["proposal_build_milliseconds"]) for row in selected
                ),
                "mean_calibration_milliseconds": statistics.mean(
                    float(row["calibration_milliseconds"]) for row in selected
                ),
                "mean_trial_milliseconds": statistics.mean(
                    float(row["trial_milliseconds"]) for row in selected
                ),
            }
        )
    return output


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def make_plot(path: Path, summary: list[dict[str, object]]) -> None:
    import matplotlib.pyplot as plt

    colors = {"ee-hadrons": "#2563eb", "softqcd": "#d97706"}
    figure, axes = plt.subplots(1, 2, figsize=(9.5, 4.0))
    for mode, color in colors.items():
        selected = [row for row in summary if row["mode"] == mode]
        factors = [float(row["safety_factor"]) for row in selected]
        axes[0].plot(
            factors,
            [row["maximum_trial_bound_ratio"] for row in selected],
            marker="o",
            color=color,
            label=mode,
        )
        axes[1].plot(
            factors,
            [row["mean_expected_acceptance"] for row in selected],
            marker="o",
            color=color,
            label=mode,
        )
    axes[0].axhline(1.0, color="black", linewidth=0.9, linestyle="--")
    axes[0].set(xlabel="safety factor", ylabel="largest trial weight / bound")
    axes[1].set(xlabel="safety factor", ylabel="mean expected acceptance")
    axes[1].set_ylim(bottom=0.0)
    axes[0].legend()
    figure.suptitle("Independent Pythia test of empirical analytic leg envelope")
    figure.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=180)
    plt.close(figure)


def main() -> int:
    args = parse_args()
    rows = load_rows(args.inputs)
    if not rows:
        raise ValueError("no diagnostic rows found")
    summary = make_summary(rows)
    write_csv(args.summary, summary)
    if args.plot:
        make_plot(args.plot, summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
