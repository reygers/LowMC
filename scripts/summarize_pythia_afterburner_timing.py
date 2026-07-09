#!/usr/bin/env python3
"""Summarize per-event timing from Pythia-produced afterburner validations."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path

STAGE_SECOND_KEYS = [
    "select_seconds",
    "yield_seconds",
    "poisson_seconds",
    "shape_seconds",
    "shape_grid_seconds",
    "shape_base_probe_seconds",
    "shape_calibration_seconds",
    "shape_probe_plan_seconds",
    "shape_boosted_probe_seconds",
    "shape_envelope_build_seconds",
    "shape_trial_seconds",
    "shape_finalization_seconds",
    "write_seconds",
]

TOP_LEVEL_STAGE_LABELS = [
    ("selection", "select_seconds"),
    ("yield", "yield_seconds"),
    ("poisson", "poisson_seconds"),
    ("shape", "shape_seconds"),
    ("write", "write_seconds"),
]

SHAPE_SUBSTAGE_LABELS = [
    ("grid/refinement/index", "shape_grid_seconds"),
    ("base probing", "shape_base_probe_seconds"),
    ("calibration trials", "shape_calibration_seconds"),
    ("probe planning", "shape_probe_plan_seconds"),
    ("boosted-cell probing", "shape_boosted_probe_seconds"),
    ("envelope build", "shape_envelope_build_seconds"),
    ("rejection trials", "shape_trial_seconds"),
    ("finalization", "shape_finalization_seconds"),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Aggregate validation_events.csv files written by "
            "scripts/validate_pythia_afterburner.py and summarize runtime "
            "versus event-complexity proxies."
        )
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        type=Path,
        help="validation_events.csv files or directories containing them",
    )
    parser.add_argument(
        "--output-csv",
        type=Path,
        default=None,
        help="optional combined per-event CSV to write",
    )
    parser.add_argument(
        "--summary-md",
        type=Path,
        default=None,
        help="optional Markdown summary to write",
    )
    parser.add_argument(
        "--top",
        type=int,
        default=5,
        help="number of slowest events to show",
    )
    return parser.parse_args()


def event_csv_path(path: Path) -> Path:
    if path.is_dir():
        return path / "validation_events.csv"
    return path


def read_rows(paths: list[Path]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for original in paths:
        path = event_csv_path(original)
        if not path.exists():
            raise FileNotFoundError(path)
        with path.open(newline="") as handle:
            for row in csv.DictReader(handle):
                row = dict(row)
                row["source"] = str(path)
                rows.append(row)
    return rows


def finite_float(row: dict[str, str], key: str) -> float | None:
    text = row.get(key, "")
    if text == "":
        return None
    value = float(text)
    if not math.isfinite(value):
        return None
    return value


def require_float(row: dict[str, str], key: str) -> float:
    value = finite_float(row, key)
    if value is None:
        raise ValueError(f"row is missing finite {key}: {row}")
    return value


def pearson(rows: list[dict[str, str]], x_key: str, y_key: str = "event_seconds") -> float:
    pairs: list[tuple[float, float]] = []
    for row in rows:
        x = finite_float(row, x_key)
        y = finite_float(row, y_key)
        if x is not None and y is not None:
            pairs.append((x, y))
    if len(pairs) < 2:
        return math.nan

    mean_x = sum(x for x, _ in pairs) / len(pairs)
    mean_y = sum(y for _, y in pairs) / len(pairs)
    cov = sum((x - mean_x) * (y - mean_y) for x, y in pairs)
    var_x = sum((x - mean_x) ** 2 for x, _ in pairs)
    var_y = sum((y - mean_y) ** 2 for _, y in pairs)
    if var_x <= 0.0 or var_y <= 0.0:
        return math.nan
    return cov / math.sqrt(var_x * var_y)


def numeric_summary(rows: list[dict[str, str]]) -> dict[str, float]:
    seconds = [require_float(row, "event_seconds") for row in rows]
    selected = [require_float(row, "selected_outgoing_charged") for row in rows]
    requested = [require_float(row, "requested_photons") for row in rows]
    retained = [require_float(row, "retained_photons") for row in rows]
    ratios = [require_float(row, "max_trial_bound_ratio") for row in rows]
    summary = {
        "events": float(len(rows)),
        "event_seconds_sum": sum(seconds),
        "event_seconds_mean": sum(seconds) / len(seconds),
        "event_seconds_max": max(seconds),
        "selected_outgoing_charged_min": min(selected),
        "selected_outgoing_charged_max": max(selected),
        "requested_photons_sum": sum(requested),
        "requested_photons_max": max(requested),
        "retained_photons_sum": sum(retained),
        "retained_photons_max": max(retained),
        "max_trial_bound_ratio_max": max(ratios),
        "corr_seconds_selected_outgoing_charged": pearson(rows, "selected_outgoing_charged"),
        "corr_seconds_requested_photons": pearson(rows, "requested_photons"),
        "corr_seconds_retained_photons": pearson(rows, "retained_photons"),
        "corr_seconds_max_trial_bound_ratio": pearson(rows, "max_trial_bound_ratio"),
    }
    for key in STAGE_SECOND_KEYS:
        values = [
            value
            for row in rows
            if (value := finite_float(row, key)) is not None
        ]
        stage_sum = sum(values)
        summary[f"{key}_sum"] = stage_sum
        summary[f"{key}_max"] = max(values) if values else math.nan
        summary[f"{key}_fraction"] = (
            stage_sum / summary["event_seconds_sum"]
            if values and summary["event_seconds_sum"] > 0.0 else math.nan
        )
        summary[f"{key}_count"] = float(len(values))
    return summary


def write_combined_csv(path: Path, rows: list[dict[str, str]]) -> None:
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_summary_md(
    path: Path,
    rows: list[dict[str, str]],
    summary: dict[str, float],
    top: int,
) -> None:
    slowest = sorted(rows, key=lambda row: require_float(row, "event_seconds"), reverse=True)[:top]
    lines = [
        "# Pythia Afterburner Timing Summary",
        "",
        f"Events: `{int(summary['events'])}`.",
        f"Total event time: `{summary['event_seconds_sum']:.3f} s`.",
        f"Mean event time: `{summary['event_seconds_mean']:.3f} s`.",
        f"Slowest event: `{summary['event_seconds_max']:.3f} s`.",
        "",
        "## Complexity Range",
        "",
        (
            f"- Selected outgoing charged multiplicity: "
            f"`{summary['selected_outgoing_charged_min']:.0f}.."
            f"{summary['selected_outgoing_charged_max']:.0f}`."
        ),
        (
            f"- Requested photons: total `{summary['requested_photons_sum']:.0f}`, "
            f"maximum `{summary['requested_photons_max']:.0f}`."
        ),
        (
            f"- Retained photons: total `{summary['retained_photons_sum']:.0f}`, "
            f"maximum `{summary['retained_photons_max']:.0f}`."
        ),
        f"- Maximum trial/bound ratio: `{summary['max_trial_bound_ratio_max']:.3g}`.",
        "",
    ]

    if summary["shape_seconds_count"] > 0:
        lines.extend(
            [
                "## Stage Timing",
                "",
                "| stage | total seconds | max event seconds | fraction of event time |",
                "| --- | ---: | ---: | ---: |",
            ]
        )
        for label, key in TOP_LEVEL_STAGE_LABELS:
            lines.append(
                "| "
                f"{label} | "
                f"{summary[f'{key}_sum']:.3f} | "
                f"{summary[f'{key}_max']:.3f} | "
                f"{summary[f'{key}_fraction']:.1%} |"
            )
        lines.append("")

    if summary["shape_trial_seconds_count"] > 0:
        lines.extend(
            [
                "## Shape Substage Timing",
                "",
                "| shape substage | total seconds | max event seconds | fraction of shape time |",
                "| --- | ---: | ---: | ---: |",
            ]
        )
        shape_total = summary["shape_seconds_sum"]
        for label, key in SHAPE_SUBSTAGE_LABELS:
            substage_sum = summary[f"{key}_sum"]
            fraction = substage_sum / shape_total if shape_total > 0.0 else math.nan
            fraction_text = "n/a" if math.isnan(fraction) else f"{fraction:.1%}"
            lines.append(
                "| "
                f"{label} | "
                f"{substage_sum:.3f} | "
                f"{summary[f'{key}_max']:.3f} | "
                f"{fraction_text} |"
            )
        lines.append("")

        counter_ranges = []
        for label, key in [
            ("grid cells", "shape_grid_cell_count"),
            ("probe evaluations", "shape_probe_evaluations"),
            ("calibration trials", "shape_calibration_count"),
            ("rejection trials", "shape_trial_count"),
            ("accepted candidates", "shape_accepted_count"),
            ("suspect cells", "shape_suspect_cell_count"),
            ("promoted cells", "shape_promoted_cell_count"),
        ]:
            values = [
                value
                for row in rows
                if (value := finite_float(row, key)) is not None
            ]
            if values:
                counter_ranges.append((label, min(values), max(values)))
        if counter_ranges:
            lines.extend(
                [
                    "## Shape Counter Ranges",
                    "",
                    "| counter | min | max |",
                    "| --- | ---: | ---: |",
                ]
            )
            for label, minimum, maximum in counter_ranges:
                lines.append(f"| {label} | {minimum:.0f} | {maximum:.0f} |")
            lines.append("")

    lines.extend(
        [
            "## Correlations With Event Time",
            "",
            "| proxy | Pearson r |",
            "| --- | ---: |",
        ]
    )
    for label, key in [
        ("selected outgoing charged", "corr_seconds_selected_outgoing_charged"),
        ("requested photons", "corr_seconds_requested_photons"),
        ("retained photons", "corr_seconds_retained_photons"),
        ("max trial/bound ratio", "corr_seconds_max_trial_bound_ratio"),
    ]:
        value = summary[key]
        formatted = "n/a" if math.isnan(value) else f"{value:.3f}"
        lines.append(f"| {label} | {formatted} |")

    lines.extend(
        [
            "",
            "## Slowest Events",
            "",
            (
                "| source | event | seconds | selected outgoing charged | "
                "requested | retained | max trial/bound ratio | status |"
            ),
            "| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
        ]
    )
    for row in slowest:
        lines.append(
            "| "
            f"{Path(row['source']).parent.name} | "
            f"{row['event_index']} | "
            f"{require_float(row, 'event_seconds'):.3f} | "
            f"{require_float(row, 'selected_outgoing_charged'):.0f} | "
            f"{require_float(row, 'requested_photons'):.0f} | "
            f"{require_float(row, 'retained_photons'):.0f} | "
            f"{require_float(row, 'max_trial_bound_ratio'):.3g} | "
            f"{row.get('shape_status', '')} |"
        )

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def print_summary(summary: dict[str, float]) -> None:
    correlations = [
        ("selected", summary["corr_seconds_selected_outgoing_charged"]),
        ("requested", summary["corr_seconds_requested_photons"]),
        ("ratio", summary["corr_seconds_max_trial_bound_ratio"]),
    ]
    correlation_text = " ".join(
        f"corr_{name}={'nan' if math.isnan(value) else f'{value:.3g}'}"
        for name, value in correlations
    )
    shape_fraction = summary["shape_seconds_fraction"]
    shape_fraction_text = "nan" if math.isnan(shape_fraction) else f"{shape_fraction:.3g}"
    dominant_shape_label = "none"
    dominant_shape_seconds = math.nan
    if summary["shape_trial_seconds_count"] > 0:
        dominant_shape_label, dominant_shape_key = max(
            SHAPE_SUBSTAGE_LABELS,
            key=lambda item: summary[f"{item[1]}_sum"],
        )
        dominant_shape_seconds = summary[f"{dominant_shape_key}_sum"]
    dominant_shape_text = (
        "none"
        if math.isnan(dominant_shape_seconds)
        else f"{dominant_shape_label}:{dominant_shape_seconds:.6g}"
    )
    print(
        f"events={int(summary['events'])} "
        f"event_seconds_sum={summary['event_seconds_sum']:.6g} "
        f"event_seconds_mean={summary['event_seconds_mean']:.6g} "
        f"event_seconds_max={summary['event_seconds_max']:.6g} "
        f"selected_outgoing_range="
        f"{summary['selected_outgoing_charged_min']:.0f}.."
        f"{summary['selected_outgoing_charged_max']:.0f} "
        f"requested_photons_sum={summary['requested_photons_sum']:.0f} "
        f"max_trial_bound_ratio={summary['max_trial_bound_ratio_max']:.6g} "
        f"shape_seconds_fraction={shape_fraction_text} "
        f"dominant_shape_substage={dominant_shape_text} "
        f"{correlation_text}"
    )


def main() -> int:
    args = parse_args()
    if args.top <= 0:
        print("--top must be positive", file=sys.stderr)
        return 1
    try:
        rows = read_rows(args.inputs)
        if not rows:
            raise ValueError("no event rows found")
        for row in rows:
            require_float(row, "event_seconds")
        summary = numeric_summary(rows)
    except (OSError, ValueError, KeyError) as error:
        print(f"Could not summarize Pythia afterburner timing: {error}", file=sys.stderr)
        return 1

    if args.output_csv is not None:
        write_combined_csv(args.output_csv, rows)
    if args.summary_md is not None:
        write_summary_md(args.summary_md, rows, summary, args.top)
    print_summary(summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
