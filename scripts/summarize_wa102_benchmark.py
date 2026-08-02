#!/usr/bin/env python3
"""Summarize the paired WA102 direct-integral/LowMC benchmark."""

import argparse
import csv
import math
from pathlib import Path


FIELDS = {
    "direct_truth": "direct_truth_yield_event",
    "direct_isolated": "direct_isolated_yield_event",
    "direct_truth_smeared": "direct_truth_smeared_yield_event",
    "direct_smeared": "direct_isolated_smeared_yield_event",
    "direct_coarse_truth_smeared": (
        "direct_coarse_truth_smeared_yield_event"
    ),
    "direct_coarse_smeared": (
        "direct_coarse_isolated_smeared_yield_event"
    ),
    "mc_truth": "lowmc_truth_yield_event",
    "mc_isolated": "lowmc_isolated_yield_event",
    "mc_truth_smeared": "lowmc_truth_smeared_yield_event",
    "mc_truth_smeared_se": "lowmc_truth_smeared_standard_error",
    "mc_smeared": "lowmc_isolated_smeared_yield_event",
    "mc_smeared_se": "lowmc_isolated_smeared_standard_error",
}

COMPARISON_BIN_INDICES = tuple(range(1, 9))


def parse_args():
    project_root = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser(
        description=(
            "Aggregate event-level WA102 benchmark CSV files, compare with "
            "digitized Figure 2d, and optionally make a final plot."
        )
    )
    parser.add_argument("inputs", nargs="+", type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--plot", type=Path)
    parser.add_argument(
        "--isolation-plot",
        type=Path,
        help="Write the approximate isolation study as a separate diagnostic.",
    )
    parser.add_argument(
        "--reference",
        type=Path,
        default=(
            project_root
            / "fixtures/reference/wa102_figure2d_pt_digitized.csv"
        ),
    )
    return parser.parse_args()


def read_reference(path):
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != 10:
        raise ValueError(
            f"expected 10 WA102 reference bins in {path}, found {len(rows)}"
        )
    return rows


def read_events(paths):
    events = []
    metadata_fields = (
        "beam_momentum_GeV",
        "lead_distance_m",
        "isolation_distance_m",
        "magnetic_field_T",
        "sheet_half_size_m",
        "pt_resolution_GeV",
    )
    metadata = None
    for input_index, path in enumerate(paths):
        grouped = {}
        row_counts = {}
        with path.open(newline="", encoding="utf-8") as stream:
            for row in csv.DictReader(stream):
                key = (
                    input_index,
                    row["seed"],
                    row["event_index"],
                    row["sampling_seed"],
                )
                row_metadata = tuple(row[name] for name in metadata_fields)
                if metadata is None:
                    metadata = row_metadata
                elif row_metadata != metadata:
                    raise ValueError(
                        "input files use different WA102 geometry settings"
                    )
                if key not in grouped:
                    grouped[key] = {
                        name: [0.0] * 10 for name in FIELDS
                    }
                    grouped[key]["mc_smeared_total_se"] = float(
                        row[
                            "lowmc_isolated_smeared_total_standard_error_event"
                        ]
                    )
                    grouped[key]["mc_truth_smeared_total_se"] = float(
                        row[
                            "lowmc_truth_smeared_total_standard_error_event"
                        ]
                    )
                    grouped[key]["mc_truth_smeared_comparison_se"] = float(
                        row[
                            "lowmc_truth_smeared_comparison_standard_error_event"
                        ]
                    )
                    grouped[key]["track_hits"] = int(row["track_hits"])
                    grouped[key]["charged_multiplicity"] = int(
                        row["charged_multiplicity"]
                    )
                    row_counts[key] = 0
                event = grouped[key]
                bin_index = int(row["bin_index"])
                for name, column in FIELDS.items():
                    event[name][bin_index] = float(row[column])
                row_counts[key] += 1
        for key, event in grouped.items():
            if row_counts[key] != 10:
                raise ValueError(
                    f"expected 10 rows for event {key[2]} in {path}, "
                    f"found {row_counts[key]}"
                )
            events.append(event)
    if not events or metadata is None:
        raise ValueError("no WA102 event rows found")
    return events, dict(zip(metadata_fields, map(float, metadata)))


def mean_value(events, field, bin_index=None):
    values = [
        (
            sum(event[field])
            if bin_index is None
            else event[field][bin_index]
        )
        for event in events
    ]
    return sum(values) / len(values)


def jackknife_mean_se(events, field, bin_index=None):
    if len(events) < 2:
        return math.nan
    values = [
        (
            sum(event[field])
            if bin_index is None
            else event[field][bin_index]
        )
        for event in events
    ]
    mean = sum(values) / len(values)
    return math.sqrt(
        sum((value - mean) ** 2 for value in values)
        / (len(values) * (len(values) - 1))
    )


def selected_mean(events, field, indices):
    return sum(
        sum(event[field][index] for index in indices) for event in events
    ) / len(events)


def selected_jackknife_mean_se(events, field, indices):
    if len(events) < 2:
        return math.nan
    values = [
        sum(event[field][index] for index in indices) for event in events
    ]
    mean = sum(values) / len(values)
    return math.sqrt(
        sum((value - mean) ** 2 for value in values)
        / (len(values) * (len(values) - 1))
    )


def mc_integration_se(
    events,
    bin_index=None,
    field="mc_truth_smeared_se",
    total_field="mc_truth_smeared_total_se",
):
    if bin_index is None:
        variance = sum(event[total_field] ** 2 for event in events)
    else:
        variance = sum(
            event[field][bin_index] ** 2 for event in events
        )
    return math.sqrt(variance) / len(events)


def safe_ratio(numerator, denominator):
    return numerator / denominator if denominator else math.nan


def write_summary(path, events, reference):
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    for bin_index, published in enumerate(reference):
        direct = mean_value(events, "direct_truth_smeared", bin_index)
        coarse = mean_value(
            events, "direct_coarse_truth_smeared", bin_index
        )
        mc = mean_value(events, "mc_truth_smeared", bin_index)
        brems = float(published["bremsstrahlung_per_event_digitized"])
        rows.append(
            {
                "bin_index": bin_index,
                "pt_min_GeV": published["pt_min_GeV"],
                "pt_max_GeV": published["pt_max_GeV"],
                "wa102_data_minus_background_per_event": float(
                    published["data_minus_background_per_event"]
                ),
                "wa102_data_statistical_uncertainty_per_event": float(
                    published["data_statistical_uncertainty_per_event"]
                ),
                "wa102_bremsstrahlung_per_event_digitized": brems,
                "direct_truth_per_event": mean_value(
                    events, "direct_truth", bin_index
                ),
                "direct_truth_event_jackknife_se": jackknife_mean_se(
                    events, "direct_truth", bin_index
                ),
                "direct_truth_smeared_per_event": mean_value(
                    events, "direct_truth_smeared", bin_index
                ),
                "direct_truth_smeared_event_jackknife_se": (
                    jackknife_mean_se(
                        events, "direct_truth_smeared", bin_index
                    )
                ),
                "direct_isolated_per_event": mean_value(
                    events, "direct_isolated", bin_index
                ),
                "direct_isolated_smeared_per_event": mean_value(
                    events, "direct_smeared", bin_index
                ),
                "direct_isolated_smeared_event_jackknife_se": (
                    jackknife_mean_se(events, "direct_smeared", bin_index)
                ),
                "direct_coarse_truth_smeared_per_event": coarse,
                "direct_coarse_relative_difference": (
                    abs(direct - coarse) / abs(direct) if direct else 0.0
                ),
                "lowmc_truth_per_event": mean_value(
                    events, "mc_truth", bin_index
                ),
                "lowmc_isolated_per_event": mean_value(
                    events, "mc_isolated", bin_index
                ),
                "lowmc_isolated_smeared_per_event": mean_value(
                    events, "mc_smeared", bin_index
                ),
                "lowmc_truth_smeared_per_event": mc,
                "lowmc_truth_smeared_event_jackknife_se": jackknife_mean_se(
                    events, "mc_truth_smeared", bin_index
                ),
                "lowmc_truth_smeared_integration_se": mc_integration_se(
                    events, bin_index
                ),
                "lowmc_over_direct": safe_ratio(mc, direct),
                "direct_over_wa102_brems_digitized": safe_ratio(
                    direct, brems
                ),
                "lowmc_over_wa102_brems_digitized": safe_ratio(
                    mc, brems
                ),
            }
        )
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return rows


def make_plot(path, rows, events, metadata):
    import matplotlib.pyplot as plt

    centers = [
        500.0 * (float(row["pt_min_GeV"]) + float(row["pt_max_GeV"]))
        for row in rows
    ]
    data = [
        row["wa102_data_minus_background_per_event"] for row in rows
    ]
    data_error = [
        row["wa102_data_statistical_uncertainty_per_event"]
        for row in rows
    ]
    brems = [
        row["wa102_bremsstrahlung_per_event_digitized"] for row in rows
    ]
    direct = [
        row["direct_truth_smeared_per_event"] for row in rows
    ]
    direct_error = [
        row["direct_truth_smeared_event_jackknife_se"] for row in rows
    ]
    mc = [
        row["lowmc_truth_smeared_per_event"] for row in rows
    ]
    mc_error = [
        math.hypot(
            row["lowmc_truth_smeared_event_jackknife_se"],
            row["lowmc_truth_smeared_integration_se"],
        )
        for row in rows
    ]

    path.parent.mkdir(parents=True, exist_ok=True)
    figure, (axis, ratio_axis) = plt.subplots(
        2,
        1,
        figsize=(7.4, 6.4),
        sharex=True,
        gridspec_kw={"height_ratios": [3.0, 1.15], "hspace": 0.06},
    )
    axis.errorbar(
        centers,
        data,
        xerr=1.0,
        yerr=data_error,
        fmt="o",
        color="black",
        capsize=2,
        label=r"WA102 Data $-$ BG (digitized)",
    )
    axis.plot(
        centers,
        brems,
        "^",
        markerfacecolor="none",
        color="#666666",
        label="WA102 Brems (digitized)",
    )
    axis.errorbar(
        [value - 0.22 for value in centers],
        direct,
        yerr=direct_error,
        fmt="s",
        color="#0072B2",
        capsize=2,
        label="Direct Low integral (no isolation)",
    )
    axis.errorbar(
        [value + 0.22 for value in centers],
        mc,
        yerr=mc_error,
        fmt="D",
        color="#D55E00",
        capsize=2,
        markersize=4,
        label="LowMC sampling (no isolation)",
    )
    axis.set(
        ylabel="photons / event / 2 MeV",
        xlim=(0.0, 20.0),
    )
    axis.grid(alpha=0.22)
    axis.legend(frameon=False, fontsize=8.5)

    direct_ratio = [
        safe_ratio(value, reference)
        for value, reference in zip(direct, brems)
    ]
    mc_ratio = [
        safe_ratio(value, reference)
        for value, reference in zip(mc, brems)
    ]
    ratio_axis.axhspan(
        0.95,
        1.05,
        color="0.85",
        alpha=0.8,
        label="WA102 5% Brems systematic",
    )
    ratio_axis.errorbar(
        [value - 0.18 for value in centers],
        direct_ratio,
        yerr=[
            safe_ratio(error, reference)
            for error, reference in zip(direct_error, brems)
        ],
        fmt="s",
        color="#0072B2",
        capsize=2,
    )
    ratio_axis.errorbar(
        [value + 0.18 for value in centers],
        mc_ratio,
        yerr=[
            safe_ratio(error, reference)
            for error, reference in zip(mc_error, brems)
        ],
        fmt="D",
        color="#D55E00",
        capsize=2,
        markersize=4,
    )
    ratio_axis.axhline(1.0, color="black", linewidth=1)
    for plot_axis in (axis, ratio_axis):
        plot_axis.axvspan(0.0, 2.0, color="0.92", zorder=-5)
        plot_axis.axvspan(18.0, 20.0, color="0.92", zorder=-5)
    ratio_axis.set(
        xlabel=r"$p_{\mathrm{T},\gamma}$ [MeV/c]",
        ylabel="ratio to\nWA102 Brems",
        xlim=(0.0, 20.0),
    )
    ratio_axis.grid(alpha=0.22)

    direct_truth_smeared = mean_value(events, "direct_truth_smeared")
    mc_truth_smeared = mean_value(events, "mc_truth_smeared")
    direct_window = selected_mean(
        events, "direct_truth_smeared", COMPARISON_BIN_INDICES
    )
    mc_window = selected_mean(
        events, "mc_truth_smeared", COMPARISON_BIN_INDICES
    )
    wa102_window = sum(brems[index] for index in COMPARISON_BIN_INDICES)
    axis.text(
        0.98,
        0.97,
        (
            f"{len(events)} selected Pythia events\n"
            f"0-20 MeV: direct={direct_truth_smeared:.4g}, "
            f"LowMC={mc_truth_smeared:.4g}\n"
            f"2-18 MeV: direct/WA102={direct_window/wa102_window:.3f}\n"
            f"2-18 MeV: LowMC/direct={mc_window/direct_window:.4f}"
        ),
        transform=axis.transAxes,
        ha="right",
        va="top",
        fontsize=8.2,
        bbox={
            "boxstyle": "round,pad=0.3",
            "facecolor": "white",
            "edgecolor": "0.75",
            "alpha": 0.92,
        },
    )
    figure.suptitle(
        (
            "WA102 fixed-target inner-bremsstrahlung benchmark\n"
            "no isolation; 2.9 MeV/c average resolution model"
        ),
        fontsize=13,
    )
    figure.subplots_adjust(
        left=0.14,
        right=0.98,
        bottom=0.11,
        top=0.88,
    )
    figure.savefig(path, dpi=180)
    plt.close(figure)


def make_isolation_plot(path, rows, metadata):
    import matplotlib.pyplot as plt

    centers = [
        500.0 * (float(row["pt_min_GeV"]) + float(row["pt_max_GeV"]))
        for row in rows
    ]
    truth = [row["direct_truth_smeared_per_event"] for row in rows]
    isolated = [
        row["direct_isolated_smeared_per_event"] for row in rows
    ]
    efficiency = [safe_ratio(value, base) for value, base in zip(isolated, truth)]

    path.parent.mkdir(parents=True, exist_ok=True)
    figure, (axis, efficiency_axis) = plt.subplots(
        2,
        1,
        figsize=(7.2, 6.0),
        sharex=True,
        gridspec_kw={"height_ratios": [3.0, 1.15], "hspace": 0.06},
    )
    axis.plot(centers, truth, "o-", color="#0072B2", label="No isolation")
    axis.plot(
        centers,
        isolated,
        "s-",
        color="#D55E00",
        label="Approximate isolation",
    )
    axis.set(ylabel="photons / event / 2 MeV", xlim=(0.0, 20.0))
    axis.grid(alpha=0.22)
    axis.legend(frameon=False)
    efficiency_axis.plot(centers, efficiency, "o-", color="#6A3D9A")
    efficiency_axis.axhline(1.0, color="black", linewidth=1)
    efficiency_axis.set(
        xlabel=r"$p_{\mathrm{T},\gamma}$ [MeV/c]",
        ylabel="isolation\nefficiency",
        ylim=(0.0, 1.08),
        xlim=(0.0, 20.0),
    )
    efficiency_axis.grid(alpha=0.22)
    figure.suptitle(
        (
            "WA102 approximate isolation sensitivity\n"
            f"L={metadata['lead_distance_m']:.2f} m, "
            f"d={1000.0 * metadata['isolation_distance_m']:.1f} mm, "
            f"B={metadata['magnetic_field_T']:.1f} T"
        ),
        fontsize=13,
    )
    figure.subplots_adjust(
        left=0.14,
        right=0.98,
        bottom=0.11,
        top=0.88,
    )
    figure.savefig(path, dpi=180)
    plt.close(figure)


def main():
    args = parse_args()
    events, metadata = read_events(args.inputs)
    reference = read_reference(args.reference)
    rows = write_summary(args.summary, events, reference)
    if args.plot:
        make_plot(args.plot, rows, events, metadata)
    if args.isolation_plot:
        make_isolation_plot(args.isolation_plot, rows, metadata)

    direct_truth = mean_value(events, "direct_truth")
    direct = mean_value(events, "direct_truth_smeared")
    direct_event_se = jackknife_mean_se(events, "direct_truth_smeared")
    direct_coarse = mean_value(events, "direct_coarse_truth_smeared")
    mc = mean_value(events, "mc_truth_smeared")
    mc_event_se = jackknife_mean_se(events, "mc_truth_smeared")
    mc_numeric_se = mc_integration_se(events)
    direct_isolated = mean_value(events, "direct_isolated")
    direct_isolated_smeared = mean_value(events, "direct_smeared")

    direct_window = selected_mean(
        events, "direct_truth_smeared", COMPARISON_BIN_INDICES
    )
    direct_window_event_se = selected_jackknife_mean_se(
        events, "direct_truth_smeared", COMPARISON_BIN_INDICES
    )
    direct_coarse_window = selected_mean(
        events, "direct_coarse_truth_smeared", COMPARISON_BIN_INDICES
    )
    mc_window = selected_mean(
        events, "mc_truth_smeared", COMPARISON_BIN_INDICES
    )
    mc_window_event_se = selected_jackknife_mean_se(
        events, "mc_truth_smeared", COMPARISON_BIN_INDICES
    )
    mc_window_numeric_se = math.sqrt(
        sum(
            event["mc_truth_smeared_comparison_se"] ** 2
            for event in events
        )
    ) / len(events)
    wa102_brems_displayed = sum(
        float(row["bremsstrahlung_per_event_digitized"])
        for row in reference
    )
    wa102_data = sum(
        float(row["data_minus_background_per_event"])
        for row in reference
    )
    wa102_brems_window = sum(
        float(reference[index]["bremsstrahlung_per_event_digitized"])
        for index in COMPARISON_BIN_INDICES
    )
    direct_window_shape = [
        rows[index]["direct_truth_smeared_per_event"] / direct_window
        for index in COMPARISON_BIN_INDICES
    ]
    reference_window_shape = [
        rows[index]["wa102_bremsstrahlung_per_event_digitized"]
        / wa102_brems_window
        for index in COMPARISON_BIN_INDICES
    ]
    window_shape_l1 = 0.5 * sum(
        abs(value - reference)
        for value, reference in zip(
            direct_window_shape, reference_window_shape
        )
    )
    print(
        f"events={len(events)} "
        f"mean_charged_multiplicity="
        f"{sum(event['charged_multiplicity'] for event in events)/len(events):.4g} "
        f"mean_track_hits="
        f"{sum(event['track_hits'] for event in events)/len(events):.4g} "
        f"lead_distance_m={metadata['lead_distance_m']:.6g} "
        f"magnetic_field_T={metadata['magnetic_field_T']:.6g} "
        f"direct_truth_per_event={direct_truth:.8g} "
        f"direct_truth_smeared_per_event={direct:.8g} "
        f"direct_event_jackknife_se={direct_event_se:.3g} "
        f"direct_coarse_relative_difference="
        f"{safe_ratio(abs(direct-direct_coarse), abs(direct)):.3%} "
        f"lowmc_truth_smeared_per_event={mc:.8g} "
        f"lowmc_event_jackknife_se={mc_event_se:.3g} "
        f"lowmc_integration_se={mc_numeric_se:.3g} "
        f"lowmc_over_direct={safe_ratio(mc, direct):.6f} "
        f"direct_isolated_per_event={direct_isolated:.8g} "
        f"isolation_efficiency={safe_ratio(direct_isolated, direct_truth):.6f} "
        f"direct_isolated_smeared_per_event={direct_isolated_smeared:.8g} "
        f"wa102_table1_brems_per_event=0.0116 "
        f"wa102_digitized_brems_displayed_per_event="
        f"{wa102_brems_displayed:.8g} "
        f"direct_truth_smeared_over_wa102_digitized="
        f"{safe_ratio(direct, wa102_brems_displayed):.6f} "
        f"comparison_pt_min_GeV=0.002 comparison_pt_max_GeV=0.018 "
        f"comparison_direct_per_event={direct_window:.8g} "
        f"comparison_direct_event_jackknife_se={direct_window_event_se:.3g} "
        f"comparison_direct_coarse_relative_difference="
        f"{safe_ratio(abs(direct_window-direct_coarse_window), abs(direct_window)):.3%} "
        f"comparison_lowmc_per_event={mc_window:.8g} "
        f"comparison_lowmc_event_jackknife_se={mc_window_event_se:.3g} "
        f"comparison_lowmc_integration_se={mc_window_numeric_se:.3g} "
        f"comparison_lowmc_over_direct={safe_ratio(mc_window, direct_window):.6f} "
        f"comparison_wa102_brems_per_event={wa102_brems_window:.8g} "
        f"comparison_direct_over_wa102="
        f"{safe_ratio(direct_window, wa102_brems_window):.6f} "
        f"comparison_shape_l1={window_shape_l1:.6f} "
        f"wa102_data_per_event={wa102_data:.8g}"
    )


if __name__ == "__main__":
    main()
