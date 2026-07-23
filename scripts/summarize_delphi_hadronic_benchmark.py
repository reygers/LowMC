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
                        "mc_se": [math.nan] * 25,
                    }
                    row_counts[key] = 0
                event = grouped[key]
                theta_bin = int(row["theta_bin_index"])
                event["direct"][theta_bin] += float(
                    row["direct_yield_event"]
                )
                event["coarse"][theta_bin] += float(
                    row["direct_coarse_yield_event"]
                )
                event["mc"][theta_bin] += float(
                    row["lowmc_mc_yield_event"]
                )
                event["mc_se"][theta_bin] = float(
                    row["lowmc_mc_theta_standard_error_event"]
                )
                row_counts[key] += 1
        for key, event in grouped.items():
            if row_counts[key] != 50 * 25:
                raise ValueError(
                    f"expected 1250 joint rows for event {key[2]} in "
                    f"{path}, found {row_counts[key]}"
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


def theta_mc_integration_se(events, bin_index):
    total_jets = sum(event["jets"] for event in events)
    variance = sum(event["mc_se"][bin_index] ** 2 for event in events)
    return math.sqrt(variance) / total_jets


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


def write_theta_summary(path, events, reference):
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    reference_total = sum(
        float(row["bremsstrahlung_per_jet"]) for row in reference
    )
    direct_total = ratio_of_sums(events, "direct")
    for bin_index, delphi in enumerate(reference):
        published = float(delphi["bremsstrahlung_per_jet"])
        direct = ratio_of_sums(events, "direct", bin_index)
        direct_event_se = jackknife_ratio_se(events, "direct", bin_index)
        coarse = ratio_of_sums(events, "coarse", bin_index)
        mc = ratio_of_sums(events, "mc", bin_index)
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
                    events, "mc", bin_index
                ),
                "lowmc_mc_integration_se": theta_mc_integration_se(
                    events, bin_index
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
        label="Direct Low integral (event jackknife)",
    )
    axis.errorbar(
        [value + 1.0 for value in centers],
        mc,
        yerr=mc_error,
        fmt="^",
        color="#D55E00",
        capsize=2,
        label="LowMC sampling (event + integration error)",
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
        "DELPHI hadronic angular benchmark (unscaled)",
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

    if args.joint_inputs:
        joint_events = read_joint_events(args.joint_inputs)
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

        theta_direct = ratio_of_sums(joint_events, "direct")
        theta_direct_se = jackknife_ratio_se(joint_events, "direct")
        theta_mc = ratio_of_sums(joint_events, "mc")
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
            f"theta_direct_per_jet={theta_direct:.8g} "
            f"theta_direct_event_jackknife_se={theta_direct_se:.3g} "
            f"theta_lowmc_mc_per_jet={theta_mc:.8g} "
            f"theta_lowmc_mc_over_direct="
            f"{safe_ratio(theta_mc, theta_direct):.6f} "
            f"theta_delphi_digitized_per_jet={theta_published:.8g} "
            f"theta_direct_over_delphi_digitized="
            f"{theta_direct/theta_published:.6f} "
            f"theta_shape_l1={theta_shape_l1:.6f}"
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
