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
        "--reference",
        type=Path,
        default=project_root
        / "fixtures/reference/delphi_hadronic_bremsstrahlung_table4.csv",
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
    mc_error = [row["lowmc_mc_integration_se"] for row in rows]

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
        label="LowMC sampling (integration error)",
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


def main():
    args = parse_args()
    events = read_events(args.inputs)
    reference = read_reference(args.reference)
    rows = write_summary(args.summary, events, reference)
    if args.plot:
        make_plot(args.plot, rows)

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
        f"direct_coarse_relative_difference={abs(direct-coarse)/abs(direct):.3%} "
        f"lowmc_mc_per_jet={mc:.8g} "
        f"lowmc_mc_event_jackknife_se={mc_event_se:.3g} "
        f"lowmc_mc_integration_se={mc_numeric_se:.3g} "
        f"lowmc_mc_over_direct={mc/direct:.6f} "
        f"delphi_per_jet={published:.8g} "
        f"direct_over_delphi={direct/published:.6f} "
        f"diagonal_chi2_ndf={chi2:.3f}/10"
    )


if __name__ == "__main__":
    main()
