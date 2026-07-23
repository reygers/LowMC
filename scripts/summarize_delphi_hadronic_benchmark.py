#!/usr/bin/env python3
"""Summarize the paired DELPHI hadronic direct-integral/LowMC benchmark."""

import argparse
import csv
import math
from pathlib import Path


def parse_args():
    project_root = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser(
        description=(
            "Aggregate event-level DELPHI hadronic benchmark CSV files, "
            "compare with Table 4, and optionally make an unscaled plot."
        )
    )
    parser.add_argument("inputs", nargs="+", type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--plot", type=Path)
    parser.add_argument(
        "--final-plot",
        type=Path,
        help=(
            "Write a final pT comparison with a LowMC/DELPHI ratio panel."
        ),
    )
    parser.add_argument(
        "--joint-inputs",
        nargs="+",
        type=Path,
        help="Event-level joint pT-theta CSV files from the benchmark.",
    )
    parser.add_argument("--theta-summary", type=Path)
    parser.add_argument("--theta-plot", type=Path)
    parser.add_argument(
        "--reference",
        type=Path,
        default=project_root
        / "fixtures/reference/delphi_hadronic_bremsstrahlung_table4.csv",
    )
    parser.add_argument(
        "--theta-reference",
        type=Path,
        default=project_root
        / "fixtures/reference/"
        "delphi_hadronic_bremsstrahlung_figure4b_theta_digitized.csv",
    )
    return parser.parse_args()


def read_events(paths):
    events = []
    for input_index, path in enumerate(paths):
        grouped = {}
        with path.open(newline="", encoding="utf-8") as stream:
            for row in csv.DictReader(stream):
                key = (
                    input_index,
                    row["seed"],
                    row["event_index"],
                    row["sampling_seed"],
                )
                if key not in grouped:
                    grouped[key] = {
                        "jets": int(row["accepted_jets"]),
                        "direct": [0.0] * 10,
                        "coarse": [0.0] * 10,
                        "mc": [0.0] * 10,
                        "mc_se": [0.0] * 10,
                        "mc_total_se": float(
                            row.get("lowmc_mc_total_standard_error_event", "nan")
                        ),
                    }
                event = grouped[key]
                bin_index = int(row["bin_index"])
                event["direct"][bin_index] = float(row["direct_yield_event"])
                event["coarse"][bin_index] = float(
                    row["direct_coarse_yield_event"]
                )
                event["mc"][bin_index] = float(row["lowmc_mc_yield_event"])
                event["mc_se"][bin_index] = float(
                    row["lowmc_mc_standard_error"]
                )
        events.extend(grouped.values())
    if not events:
        raise ValueError("no event rows found")
    return events


def ratio_of_sums(events, field, bin_index=None):
    numerator = 0.0
    denominator = 0
    for event in events:
        values = event[field]
        numerator += (
            sum(values) if bin_index is None else values[bin_index]
        )
        denominator += event["jets"]
    return numerator / denominator


def safe_ratio(numerator, denominator):
    return numerator / denominator if denominator else math.nan


def jackknife_ratio_se(events, field, bin_index=None):
    if len(events) < 2:
        return math.nan
    total_numerator = sum(
        sum(event[field])
        if bin_index is None
        else event[field][bin_index]
        for event in events
    )
    total_jets = sum(event["jets"] for event in events)
    leave_one_out = []
    for event in events:
        event_numerator = (
            sum(event[field])
            if bin_index is None
            else event[field][bin_index]
        )
        leave_one_out.append(
            (total_numerator - event_numerator)
            / (total_jets - event["jets"])
        )
    mean = sum(leave_one_out) / len(leave_one_out)
    return math.sqrt(
        (len(events) - 1)
        / len(events)
        * sum((value - mean) ** 2 for value in leave_one_out)
    )


def mc_integration_se(events, bin_index=None):
    total_jets = sum(event["jets"] for event in events)
    if bin_index is not None:
        variance = sum(
            event["mc_se"][bin_index] ** 2 for event in events
        )
        return math.sqrt(variance) / total_jets
    total_errors = [event["mc_total_se"] for event in events]
    if any(not math.isfinite(value) for value in total_errors):
        return math.nan
    return math.sqrt(sum(value * value for value in total_errors)) / total_jets


def read_reference(path):
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != 10:
        raise ValueError(f"expected 10 DELPHI bins in {path}, found {len(rows)}")
    return rows


def read_joint_events(paths):
    events = []
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
                if key not in grouped:
                    grouped[key] = {
                        "jets": int(row["accepted_jets"]),
                        "direct": [0.0] * 25,
                        "coarse": [0.0] * 25,
                        "mc": [0.0] * 25,
                        "mc_second_moment": [0.0] * 25,
                        "direct_full": [0.0] * 25,
                        "coarse_full": [0.0] * 25,
                        "mc_full": [0.0] * 25,
                        "mc_full_second_moment": [0.0] * 25,
                        "mc_full_se_reference": [math.nan] * 25,
                        "mc_samples": int(row["mc_samples"]),
                    }
                    row_counts[key] = 0
                event = grouped[key]
                if event["mc_samples"] != int(row["mc_samples"]):
                    raise ValueError(
                        f"inconsistent sample count for event {key[2]} "
                        f"in {path}"
                    )
                pt_bin = int(row["pt_bin_index"])
                theta_bin = int(row["theta_bin_index"])
                direct = float(row["direct_yield_event"])
                coarse = float(row["direct_coarse_yield_event"])
                mc = float(row["lowmc_mc_yield_event"])
                mc_se = float(row["lowmc_mc_joint_standard_error"])
                mc_full_se_reference = float(
                    row["lowmc_mc_theta_standard_error_event"]
                )
                previous_reference = event["mc_full_se_reference"][theta_bin]
                if math.isfinite(previous_reference) and not math.isclose(
                    previous_reference,
                    mc_full_se_reference,
                    rel_tol=1.0e-12,
                    abs_tol=1.0e-15,
                ):
                    raise ValueError(
                        f"inconsistent full-theta MC error for event "
                        f"{key[2]} in {path}"
                    )
                event["mc_full_se_reference"][
                    theta_bin
                ] = mc_full_se_reference
                mc_second_moment = (
                    event["mc_samples"] * mc_se * mc_se + mc * mc
                )
                event["direct_full"][theta_bin] += direct
                event["coarse_full"][theta_bin] += coarse
                event["mc_full"][theta_bin] += mc
                event["mc_full_second_moment"][
                    theta_bin
                ] += mc_second_moment
                if pt_bin < 10:
                    event["direct"][theta_bin] += direct
                    event["coarse"][theta_bin] += coarse
                    event["mc"][theta_bin] += mc
                    event["mc_second_moment"][
                        theta_bin
                    ] += mc_second_moment
                row_counts[key] += 1
        for key, event in grouped.items():
            if row_counts[key] != 50 * 25:
                raise ValueError(
                    f"expected 1250 joint rows for event {key[2]} in "
                    f"{path}, found {row_counts[key]}"
                )
            samples = event["mc_samples"]
            event["mc_se"] = [
                math.sqrt(
                    max(
                        0.0,
                        second_moment - mean * mean,
                    )
                    / samples
                )
                for second_moment, mean in zip(
                    event.pop("mc_second_moment"),
                    event["mc"],
                )
            ]
            event["mc_full_se"] = [
                math.sqrt(
                    max(
                        0.0,
                        second_moment - mean * mean,
                    )
                    / samples
                )
                for second_moment, mean in zip(
                    event.pop("mc_full_second_moment"),
                    event["mc_full"],
                )
            ]
            for theta_bin, (projected, reference) in enumerate(
                zip(
                    event["mc_full_se"],
                    event.pop("mc_full_se_reference"),
                )
            ):
                if not math.isclose(
                    projected,
                    reference,
                    rel_tol=1.0e-10,
                    abs_tol=1.0e-14,
                ):
                    raise ValueError(
                        f"joint MC moments do not reproduce the full-theta "
                        f"error in bin {theta_bin} for event {key[2]} "
                        f"in {path}"
                    )
            events.append(event)
    if not events:
        raise ValueError("no joint event rows found")
    return events


def read_theta_reference(path):
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != 25:
        raise ValueError(
            f"expected 25 digitized DELPHI theta bins in {path}, "
            f"found {len(rows)}"
        )
    return rows


def theta_mc_integration_se(events, bin_index, field="mc_se"):
    total_jets = sum(event["jets"] for event in events)
    variance = sum(event[field][bin_index] ** 2 for event in events)
    return math.sqrt(variance) / total_jets


def validate_pt80_projection(events, joint_events):
    if len(events) != len(joint_events):
        raise ValueError(
            "primary and joint inputs contain different event counts"
        )
    if sum(event["jets"] for event in events) != sum(
        event["jets"] for event in joint_events
    ):
        raise ValueError(
            "primary and joint inputs contain different accepted-jet counts"
        )
    for field in ("direct", "coarse", "mc"):
        primary_total = sum(sum(event[field]) for event in events)
        joint_total = sum(sum(event[field]) for event in joint_events)
        tolerance = 1.0e-12 * max(1.0, abs(primary_total))
        if abs(primary_total - joint_total) > tolerance:
            raise ValueError(
                f"joint pT < 80 MeV theta projection disagrees with "
                f"the primary {field} result"
            )


def write_summary(path, events, reference):
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    for bin_index, delphi in enumerate(reference):
        direct = ratio_of_sums(events, "direct", bin_index)
        direct_event_se = jackknife_ratio_se(events, "direct", bin_index)
        coarse = ratio_of_sums(events, "coarse", bin_index)
        mc = ratio_of_sums(events, "mc", bin_index)
        mc_event_se = jackknife_ratio_se(events, "mc", bin_index)
        mc_numeric_se = mc_integration_se(events, bin_index)
        published = float(delphi["bremsstrahlung_per_jet"])
        published_stat = float(
            delphi["bremsstrahlung_statistical_uncertainty"]
        )
        published_syst = float(
            delphi["bremsstrahlung_systematic_uncertainty"]
        )
        comparison_sigma = math.sqrt(
            direct_event_se**2 + published_stat**2 + published_syst**2
        )
        rows.append(
            {
                "bin_index": bin_index,
                "pt_min_GeV": delphi["pt_min_GeV"],
                "pt_max_GeV": delphi["pt_max_GeV"],
                "delphi_bremsstrahlung_per_jet": published,
                "delphi_statistical_uncertainty": published_stat,
                "delphi_systematic_uncertainty": published_syst,
                "direct_per_jet": direct,
                "direct_event_jackknife_se": direct_event_se,
                "direct_coarse_per_jet": coarse,
                "direct_coarse_relative_difference": (
                    abs(direct - coarse) / abs(direct) if direct else 0.0
                ),
                "lowmc_mc_per_jet": mc,
                "lowmc_mc_event_jackknife_se": mc_event_se,
                "lowmc_mc_integration_se": mc_numeric_se,
                "lowmc_mc_over_direct": mc / direct if direct else math.nan,
                "direct_minus_delphi_pull_diagonal": (
                    (direct - published) / comparison_sigma
                    if comparison_sigma
                    else math.nan
                ),
            }
        )
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return rows


def make_plot(path, rows):
    import matplotlib.pyplot as plt

    centers = [
        500.0 * (float(row["pt_min_GeV"]) + float(row["pt_max_GeV"]))
        for row in rows
    ]
    delphi = [row["delphi_bremsstrahlung_per_jet"] for row in rows]
    delphi_error = [
        math.hypot(
            row["delphi_statistical_uncertainty"],
            row["delphi_systematic_uncertainty"],
        )
        for row in rows
    ]
    direct = [row["direct_per_jet"] for row in rows]
    direct_error = [row["direct_event_jackknife_se"] for row in rows]
    mc = [row["lowmc_mc_per_jet"] for row in rows]
    mc_error = [
        math.hypot(
            row["lowmc_mc_event_jackknife_se"],
            row["lowmc_mc_integration_se"],
        )
        for row in rows
    ]

    path.parent.mkdir(parents=True, exist_ok=True)
    figure, axis = plt.subplots(figsize=(7.2, 4.8))
    axis.errorbar(
        [value - 0.7 for value in centers],
        delphi,
        yerr=delphi_error,
        fmt="o",
        color="black",
        capsize=2,
        label="DELPHI Table 4 Brems",
    )
    axis.errorbar(
        centers,
        direct,
        yerr=direct_error,
        fmt="s",
        color="#0072B2",
        capsize=2,
        label="Direct Low integral (event jackknife)",
    )
    axis.errorbar(
        [value + 0.7 for value in centers],
        mc,
        yerr=mc_error,
        fmt="^",
        color="#D55E00",
        capsize=2,
        label="LowMC sampling (event + integration error)",
    )
    axis.set(
        xlabel=r"$p_{\mathrm{T},\gamma}$ relative to closest jet [MeV]",
        ylabel="photons / accepted jet / 8 MeV",
        xlim=(0.0, 80.0),
        title="Generator-level DELPHI hadronic benchmark (unscaled)",
    )
    axis.grid(alpha=0.22)
    axis.legend(frameon=False, fontsize=9)
    figure.tight_layout()
    figure.savefig(path, dpi=180)
    plt.close(figure)


def make_final_plot(path, rows, events):
    import matplotlib.pyplot as plt

    centers = [
        500.0 * (float(row["pt_min_GeV"]) + float(row["pt_max_GeV"]))
        for row in rows
    ]
    delphi = [row["delphi_bremsstrahlung_per_jet"] for row in rows]
    delphi_error = [
        math.hypot(
            row["delphi_statistical_uncertainty"],
            row["delphi_systematic_uncertainty"],
        )
        for row in rows
    ]
    direct = [row["direct_per_jet"] for row in rows]
    direct_error = [row["direct_event_jackknife_se"] for row in rows]
    mc = [row["lowmc_mc_per_jet"] for row in rows]
    mc_error = [
        math.hypot(
            row["lowmc_mc_event_jackknife_se"],
            row["lowmc_mc_integration_se"],
        )
        for row in rows
    ]
    direct_ratio = [
        value / reference for value, reference in zip(direct, delphi)
    ]
    mc_ratio = [
        value / reference for value, reference in zip(mc, delphi)
    ]
    direct_ratio_error = [
        ratio
        * math.hypot(
            prediction_error / prediction,
            reference_error / reference,
        )
        for ratio, prediction, prediction_error, reference, reference_error
        in zip(
            direct_ratio,
            direct,
            direct_error,
            delphi,
            delphi_error,
        )
    ]
    mc_ratio_error = [
        ratio
        * math.hypot(
            prediction_error / prediction,
            reference_error / reference,
        )
        for ratio, prediction, prediction_error, reference, reference_error
        in zip(mc_ratio, mc, mc_error, delphi, delphi_error)
    ]

    path.parent.mkdir(parents=True, exist_ok=True)
    figure, (axis, ratio_axis) = plt.subplots(
        2,
        1,
        figsize=(7.2, 6.2),
        sharex=True,
        gridspec_kw={"height_ratios": [3.0, 1.15], "hspace": 0.06},
    )
    axis.errorbar(
        [value - 0.7 for value in centers],
        delphi,
        yerr=delphi_error,
        fmt="o",
        color="black",
        capsize=2,
        label="DELPHI Table 4 Brems",
    )
    axis.errorbar(
        centers,
        direct,
        yerr=direct_error,
        fmt="s",
        color="#0072B2",
        capsize=2,
        label="Direct Low integral",
    )
    axis.errorbar(
        [value + 0.7 for value in centers],
        mc,
        yerr=mc_error,
        fmt="^",
        color="#D55E00",
        capsize=2,
        label="LowMC sampling",
    )
    axis.set(
        ylabel="photons / accepted jet / 8 MeV",
        xlim=(0.0, 80.0),
    )
    axis.grid(alpha=0.22)
    axis.legend(frameon=False, fontsize=9)

    ratio_axis.errorbar(
        [value - 0.35 for value in centers],
        direct_ratio,
        yerr=direct_ratio_error,
        fmt="s",
        color="#0072B2",
        capsize=2,
        label="direct / DELPHI",
    )
    ratio_axis.errorbar(
        [value + 0.35 for value in centers],
        mc_ratio,
        yerr=mc_ratio_error,
        fmt="^",
        color="#D55E00",
        capsize=2,
        label="LowMC / DELPHI",
    )
    ratio_axis.axhline(1.0, color="black", linewidth=1)
    ratio_axis.set(
        xlabel=r"$p_{\mathrm{T},\gamma}$ relative to closest jet [MeV]",
        ylabel="ratio to\nDELPHI",
        xlim=(0.0, 80.0),
    )
    ratio_axis.grid(alpha=0.22)

    delphi_total = sum(delphi)
    direct_total = ratio_of_sums(events, "direct")
    mc_total = ratio_of_sums(events, "mc")
    axis.text(
        0.025,
        0.045,
        (
            f"{len(events)} events, "
            f"{sum(event['jets'] for event in events)} jets\n"
            f"direct / DELPHI = {direct_total/delphi_total:.4f}\n"
            f"LowMC / DELPHI = {mc_total/delphi_total:.4f}"
        ),
        transform=axis.transAxes,
        fontsize=8.5,
        va="bottom",
        bbox={
            "boxstyle": "round,pad=0.3",
            "facecolor": "white",
            "edgecolor": "0.75",
            "alpha": 0.9,
        },
    )
    figure.suptitle(
        "LowMC reproduction of DELPHI inner bremsstrahlung (unscaled)",
        fontsize=15,
    )
    figure.subplots_adjust(
        left=0.14,
        right=0.98,
        bottom=0.11,
        top=0.92,
    )
    figure.savefig(path, dpi=180)
    plt.close(figure)


def write_theta_summary(path, events, reference):
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    reference_total = sum(
        float(row["bremsstrahlung_per_jet"]) for row in reference
    )
    direct_total = ratio_of_sums(events, "direct_full")
    direct_pt80_total = ratio_of_sums(events, "direct")
    for bin_index, delphi in enumerate(reference):
        published = float(delphi["bremsstrahlung_per_jet"])
        direct = ratio_of_sums(events, "direct_full", bin_index)
        direct_event_se = jackknife_ratio_se(
            events, "direct_full", bin_index
        )
        coarse = ratio_of_sums(events, "coarse_full", bin_index)
        mc = ratio_of_sums(events, "mc_full", bin_index)
        rows.append(
            {
                "theta_bin_index": bin_index,
                "theta_min_rad": delphi["theta_min_rad"],
                "theta_max_rad": delphi["theta_max_rad"],
                "delphi_bremsstrahlung_per_jet_digitized": published,
                "direct_per_jet": direct,
                "direct_event_jackknife_se": direct_event_se,
                "direct_coarse_per_jet": coarse,
                "direct_coarse_relative_difference": (
                    abs(direct - coarse) / abs(direct) if direct else 0.0
                ),
                "lowmc_mc_per_jet": mc,
                "lowmc_mc_event_jackknife_se": jackknife_ratio_se(
                    events, "mc_full", bin_index
                ),
                "lowmc_mc_integration_se": theta_mc_integration_se(
                    events, bin_index, "mc_full_se"
                ),
                "lowmc_mc_over_direct": mc / direct if direct else math.nan,
                "direct_over_delphi_digitized": (
                    direct / published if published else math.nan
                ),
                "direct_shape_fraction": safe_ratio(
                    direct, direct_total
                ),
                "delphi_digitized_shape_fraction": (
                    published / reference_total
                ),
                "direct_pt80_diagnostic_per_jet": ratio_of_sums(
                    events, "direct", bin_index
                ),
                "direct_pt80_diagnostic_event_jackknife_se": (
                    jackknife_ratio_se(events, "direct", bin_index)
                ),
                "direct_pt80_diagnostic_shape_fraction": safe_ratio(
                    ratio_of_sums(events, "direct", bin_index),
                    direct_pt80_total,
                ),
                "lowmc_mc_pt80_diagnostic_per_jet": ratio_of_sums(
                    events, "mc", bin_index
                ),
                "lowmc_mc_pt80_diagnostic_event_jackknife_se": (
                    jackknife_ratio_se(events, "mc", bin_index)
                ),
                "lowmc_mc_pt80_diagnostic_integration_se": (
                    theta_mc_integration_se(
                        events, bin_index
                    )
                ),
            }
        )
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return rows


def make_theta_plot(path, rows):
    import matplotlib.pyplot as plt

    centers = [
        500.0
        * (float(row["theta_min_rad"]) + float(row["theta_max_rad"]))
        for row in rows
    ]
    delphi = [
        row["delphi_bremsstrahlung_per_jet_digitized"] for row in rows
    ]
    direct = [row["direct_per_jet"] for row in rows]
    direct_error = [row["direct_event_jackknife_se"] for row in rows]
    mc = [row["lowmc_mc_per_jet"] for row in rows]
    mc_error = [
        math.hypot(
            row["lowmc_mc_event_jackknife_se"],
            row["lowmc_mc_integration_se"],
        )
        for row in rows
    ]

    path.parent.mkdir(parents=True, exist_ok=True)
    figure, (axis, ratio_axis) = plt.subplots(
        2,
        1,
        figsize=(7.2, 6.2),
        sharex=True,
        constrained_layout=True,
        gridspec_kw={"height_ratios": (3.0, 1.0), "hspace": 0.08},
    )
    axis.plot(
        centers,
        delphi,
        "^",
        markerfacecolor="none",
        color="black",
        label="DELPHI Figure 4b Brems (vector-digitized)",
    )
    axis.errorbar(
        [value - 1.0 for value in centers],
        direct,
        yerr=direct_error,
        fmt="s",
        color="#0072B2",
        capsize=2,
        label=r"Direct Low integral, $p_T<400$ MeV",
    )
    axis.errorbar(
        [value + 1.0 for value in centers],
        mc,
        yerr=mc_error,
        fmt="^",
        color="#D55E00",
        capsize=2,
        label=r"LowMC sampling, $p_T<400$ MeV",
    )
    upper_limit = 1.1 * max(
        max(delphi),
        max(value + error for value, error in zip(direct, direct_error)),
        max(value + error for value, error in zip(mc, mc_error)),
    )
    axis.set(
        ylabel="photons / accepted jet / 16 mrad",
        ylim=(0.0, upper_limit),
    )
    axis.grid(alpha=0.22)
    axis.legend(frameon=False, fontsize=8)

    ratio = [
        row["direct_over_delphi_digitized"] for row in rows
    ]
    ratio_error = [
        row["direct_event_jackknife_se"]
        / row["delphi_bremsstrahlung_per_jet_digitized"]
        for row in rows
    ]
    ratio_axis.errorbar(
        centers,
        ratio,
        yerr=ratio_error,
        fmt="s",
        color="#0072B2",
        capsize=2,
    )
    ratio_axis.axhline(1.0, color="black", linewidth=1)
    ratio_axis.set(
        xlabel=r"$\theta_\gamma$ relative to closest jet [mrad]",
        ylabel="direct /\nDELPHI",
        xlim=(0.0, 400.0),
    )
    ratio_axis.grid(alpha=0.22)
    figure.suptitle(
        r"DELPHI Figure 4b angular benchmark, $p_T<400$ MeV (unscaled)",
        fontsize=16,
    )
    figure.savefig(path, dpi=180)
    plt.close(figure)


def main():
    args = parse_args()
    events = read_events(args.inputs)
    reference = read_reference(args.reference)
    rows = write_summary(args.summary, events, reference)
    if args.plot:
        make_plot(args.plot, rows)
    if args.final_plot:
        make_final_plot(args.final_plot, rows, events)

    if args.joint_inputs:
        joint_events = read_joint_events(args.joint_inputs)
        validate_pt80_projection(events, joint_events)
        theta_reference = read_theta_reference(args.theta_reference)
        theta_summary_path = args.theta_summary
        if theta_summary_path is None:
            theta_summary_path = args.summary.with_name(
                args.summary.stem + "_theta" + args.summary.suffix
            )
        theta_rows = write_theta_summary(
            theta_summary_path, joint_events, theta_reference
        )
        theta_plot_path = args.theta_plot
        if theta_plot_path is None and args.plot:
            theta_plot_path = args.plot.with_name(
                args.plot.stem + "_theta" + args.plot.suffix
            )
        if theta_plot_path:
            make_theta_plot(theta_plot_path, theta_rows)

        theta_direct = ratio_of_sums(joint_events, "direct_full")
        theta_direct_se = jackknife_ratio_se(
            joint_events, "direct_full"
        )
        theta_mc = ratio_of_sums(joint_events, "mc_full")
        theta_direct_pt80 = ratio_of_sums(joint_events, "direct")
        theta_mc_pt80 = ratio_of_sums(joint_events, "mc")
        theta_published = sum(
            float(row["bremsstrahlung_per_jet"])
            for row in theta_reference
        )
        theta_shape_l1 = 0.5 * sum(
            abs(
                row["direct_shape_fraction"]
                - row["delphi_digitized_shape_fraction"]
            )
            for row in theta_rows
        )
        print(
            f"theta_events={len(joint_events)} "
            f"theta_jets={sum(event['jets'] for event in joint_events)} "
            f"theta_pt400_direct_per_jet={theta_direct:.8g} "
            f"theta_pt400_direct_event_jackknife_se={theta_direct_se:.3g} "
            f"theta_pt400_lowmc_mc_per_jet={theta_mc:.8g} "
            f"theta_pt400_lowmc_mc_over_direct="
            f"{safe_ratio(theta_mc, theta_direct):.6f} "
            f"theta_delphi_digitized_per_jet={theta_published:.8g} "
            f"theta_pt400_direct_over_delphi_digitized="
            f"{theta_direct/theta_published:.6f} "
            f"theta_pt400_shape_l1={theta_shape_l1:.6f} "
            f"theta_pt80_diagnostic_direct_per_jet="
            f"{theta_direct_pt80:.8g} "
            f"theta_pt80_diagnostic_lowmc_mc_per_jet="
            f"{theta_mc_pt80:.8g} "
            f"theta_pt80_diagnostic_lowmc_mc_over_direct="
            f"{safe_ratio(theta_mc_pt80, theta_direct_pt80):.6f}"
        )

    event_count = len(events)
    jet_count = sum(event["jets"] for event in events)
    direct = ratio_of_sums(events, "direct")
    direct_event_se = jackknife_ratio_se(events, "direct")
    coarse = ratio_of_sums(events, "coarse")
    mc = ratio_of_sums(events, "mc")
    mc_event_se = jackknife_ratio_se(events, "mc")
    mc_numeric_se = mc_integration_se(events)
    published = sum(
        float(row["bremsstrahlung_per_jet"]) for row in reference
    )
    chi2 = sum(
        row["direct_minus_delphi_pull_diagonal"] ** 2 for row in rows
    )
    print(
        f"events={event_count} jets={jet_count} "
        f"direct_per_jet={direct:.8g} "
        f"direct_event_jackknife_se={direct_event_se:.3g} "
        f"direct_coarse_relative_difference="
        f"{safe_ratio(abs(direct-coarse), abs(direct)):.3%} "
        f"lowmc_mc_per_jet={mc:.8g} "
        f"lowmc_mc_event_jackknife_se={mc_event_se:.3g} "
        f"lowmc_mc_integration_se={mc_numeric_se:.3g} "
        f"lowmc_mc_over_direct={safe_ratio(mc, direct):.6f} "
        f"delphi_per_jet={published:.8g} "
        f"direct_over_delphi={direct/published:.6f} "
        f"diagonal_chi2_ndf={chi2:.3f}/10"
    )


if __name__ == "__main__":
    main()
