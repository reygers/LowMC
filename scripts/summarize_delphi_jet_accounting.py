#!/usr/bin/env python3
"""Summarize same-event DELPHI Lund/Durham/JADE jet accounting."""

import argparse
import csv
import math
from pathlib import Path


ALGORITHMS = ("Lund", "Durham", "Jade")


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Aggregate accounting-only DELPHI benchmark CSV files and "
            "compare jet multiplicities on identical Lund-selected events."
        )
    )
    parser.add_argument("inputs", nargs="+", type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--plot", type=Path)
    return parser.parse_args()


def read_events(paths):
    events = {}
    for path in paths:
        with path.open(newline="", encoding="utf-8") as stream:
            for row in csv.DictReader(stream):
                key = (int(row["seed"]), int(row["event_index"]))
                algorithm = row["jet_algorithm"]
                if algorithm not in ALGORITHMS:
                    raise ValueError(
                        f"unexpected jet algorithm {algorithm} in {path}"
                    )
                if algorithm in events.setdefault(key, {}):
                    raise ValueError(
                        f"duplicate {algorithm} row for event {key}"
                    )
                events[key][algorithm] = row
    if not events:
        raise ValueError("no jet-accounting rows found")
    for key, algorithms in events.items():
        missing = set(ALGORITHMS) - set(algorithms)
        if missing:
            raise ValueError(
                f"event {key} is missing algorithms: {sorted(missing)}"
            )
    return events


def percentile(values, fraction):
    if not values:
        return math.nan
    ordered = sorted(values)
    position = fraction * (len(ordered) - 1)
    lower = int(math.floor(position))
    upper = int(math.ceil(position))
    if lower == upper:
        return ordered[lower]
    weight = position - lower
    return ordered[lower] * (1.0 - weight) + ordered[upper] * weight


def summarize(events):
    lund_jets = sum(
        int(rows["Lund"]["accepted_jets"])
        for rows in events.values()
    )
    result = []
    for algorithm in ALGORITHMS:
        rows = [event[algorithm] for event in events.values()]
        accepted_jets = sum(int(row["accepted_jets"]) for row in rows)
        total_jets = sum(int(row["total_jets"]) for row in rows)
        accepted_events = sum(int(row["event_accepted"]) for row in rows)
        fewer = 0
        same = 0
        more = 0
        for event in events.values():
            baseline = int(event["Lund"]["accepted_jets"])
            current = int(event[algorithm]["accepted_jets"])
            if current < baseline:
                fewer += 1
            elif current > baseline:
                more += 1
            else:
                same += 1
        nearest_angles = []
        max_angles = []
        if algorithm != "Lund":
            for row in rows:
                mean_angle = float(
                    row["mean_algorithm_to_lund_angle_rad"]
                )
                max_angle = float(
                    row["max_algorithm_to_lund_angle_rad"]
                )
                if math.isfinite(mean_angle):
                    nearest_angles.append(mean_angle)
                if math.isfinite(max_angle):
                    max_angles.append(max_angle)
        result.append(
            {
                "jet_algorithm": algorithm,
                "lund_selected_events": len(events),
                "accepted_events": accepted_events,
                "total_jets": total_jets,
                "accepted_jets": accepted_jets,
                "mean_total_jets_per_event": total_jets / len(events),
                "mean_accepted_jets_per_event": (
                    accepted_jets / len(events)
                ),
                "accepted_jet_count_over_lund": (
                    accepted_jets / lund_jets
                ),
                "events_with_fewer_accepted_jets_than_lund": fewer,
                "events_with_same_accepted_jets_as_lund": same,
                "events_with_more_accepted_jets_than_lund": more,
                "median_nearest_axis_angle_rad": percentile(
                    nearest_angles, 0.5
                ),
                "p90_nearest_axis_angle_rad": percentile(
                    nearest_angles, 0.9
                ),
                "maximum_nearest_axis_angle_rad": (
                    max(max_angles) if max_angles else math.nan
                ),
            }
        )
    return result


def write_summary(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def make_plot(path, rows):
    import matplotlib.pyplot as plt

    labels = [row["jet_algorithm"] for row in rows]
    x = list(range(len(labels)))
    total = [row["mean_total_jets_per_event"] for row in rows]
    accepted = [
        row["mean_accepted_jets_per_event"] for row in rows
    ]
    fewer = [
        row["events_with_fewer_accepted_jets_than_lund"]
        for row in rows
    ]
    same = [
        row["events_with_same_accepted_jets_as_lund"]
        for row in rows
    ]
    more = [
        row["events_with_more_accepted_jets_than_lund"]
        for row in rows
    ]

    path.parent.mkdir(parents=True, exist_ok=True)
    figure, (count_axis, event_axis) = plt.subplots(
        1, 2, figsize=(8.5, 4.2)
    )
    width = 0.36
    count_axis.bar(
        [value - width / 2 for value in x],
        total,
        width,
        color="#999999",
        label="All clustered jets",
    )
    count_axis.bar(
        [value + width / 2 for value in x],
        accepted,
        width,
        color="#0072B2",
        label="Accepted jets",
    )
    count_axis.set(
        xticks=x,
        xticklabels=labels,
        ylabel="mean jets / LUCLUS-selected event",
    )
    count_axis.grid(axis="y", alpha=0.22)
    count_axis.legend(frameon=False, fontsize=8)

    event_axis.bar(x, fewer, color="#D55E00", label="Fewer than LUCLUS")
    event_axis.bar(
        x,
        same,
        bottom=fewer,
        color="#009E73",
        label="Same as LUCLUS",
    )
    event_axis.bar(
        x,
        more,
        bottom=[
            fewer[index] + same[index] for index in x
        ],
        color="#CC79A7",
        label="More than LUCLUS",
    )
    event_axis.set(
        xticks=x,
        xticklabels=labels,
        ylabel="events",
        ylim=(0, rows[0]["lund_selected_events"]),
    )
    event_axis.grid(axis="y", alpha=0.22)
    event_axis.legend(frameon=False, fontsize=8)
    figure.suptitle(
        "Same-event DELPHI truth-jet accounting",
        fontsize=14,
    )
    figure.tight_layout()
    figure.savefig(path, dpi=180)
    plt.close(figure)


def main():
    args = parse_args()
    events = read_events(args.inputs)
    rows = summarize(events)
    write_summary(args.summary, rows)
    if args.plot:
        make_plot(args.plot, rows)
    for row in rows:
        print(
            f"algorithm={row['jet_algorithm']} "
            f"events={row['lund_selected_events']} "
            f"accepted_events={row['accepted_events']} "
            f"accepted_jets={row['accepted_jets']} "
            f"accepted_jet_count_over_lund="
            f"{row['accepted_jet_count_over_lund']:.6f} "
            f"fewer_same_more="
            f"{row['events_with_fewer_accepted_jets_than_lund']}/"
            f"{row['events_with_same_accepted_jets_as_lund']}/"
            f"{row['events_with_more_accepted_jets_than_lund']} "
            f"median_nearest_axis_angle_rad="
            f"{row['median_nearest_axis_angle_rad']:.6g} "
            f"p90_nearest_axis_angle_rad="
            f"{row['p90_nearest_axis_angle_rad']:.6g}"
        )


if __name__ == "__main__":
    main()
