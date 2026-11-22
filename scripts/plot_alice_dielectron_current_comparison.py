#!/usr/bin/env python3
"""Compare leading-Low and Lichard Eq. (2.8) ALICE benchmark spectra."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
import sys


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--leading", type=Path, required=True)
    parser.add_argument("--lichard", type=Path, required=True)
    parser.add_argument("--output-prefix", type=Path, required=True)
    parser.add_argument(
        "--layout",
        choices=("horizontal", "stacked"),
        default="horizontal",
        help="Arrange the mass and pair-pT panels side by side or vertically",
    )
    return parser.parse_args()


def read_direct(path: Path, observable: str) -> tuple[list[dict[str, str]], str]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    for decay in ("angle_integrated", "polarized"):
        selected = sorted(
            (
                row
                for row in rows
                if row["observable"] == observable
                and row["method"] == "direct_sobol"
                and row["decay"] == decay
            ),
            key=lambda row: float(row["bin_center"]),
        )
        if selected:
            return selected, decay
    return [], ""


def values(rows: list[dict[str, str]], field: str) -> list[float]:
    return [float(row[field]) for row in rows]


def main() -> int:
    args = parse_args()
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is required for plotting", file=sys.stderr)
        return 2

    plt.rcParams.update({
        "font.size": 11.5,
        "axes.labelsize": 14,
        "xtick.labelsize": 11,
        "ytick.labelsize": 11,
        "legend.fontsize": 10,
        "figure.dpi": 150,
    })
    if args.layout == "stacked":
        figure = plt.figure(figsize=(7.2, 6.6), constrained_layout=True)
        outer = figure.add_gridspec(2, 1, hspace=0.14)
    else:
        figure = plt.figure(figsize=(10.2, 4.5), constrained_layout=True)
        outer = figure.add_gridspec(1, 2, wspace=0.16)

    specifications = (
        (
            "mass",
            r"$m_{ee}$ (GeV/$c^2$)",
            r"$d^2\sigma/(dm_{ee}\,dy)$ (mb/(GeV/$c^2$))",
            (0.14, 0.61),
        ),
        (
            "qT",
            r"$p_{T,ee}$ (GeV/$c$)",
            r"$d^2\sigma/(dp_{T,ee}\,dy)$ (mb/(GeV/$c$))",
            (0.0, 0.81),
        ),
    )

    for index, (observable, xlabel, ylabel, xlim) in enumerate(specifications):
        cell = outer[index, 0] if args.layout == "stacked" else outer[0, index]
        inner = cell.subgridspec(2, 1, height_ratios=(2.6, 1.25), hspace=0.12)
        spectrum_axis = figure.add_subplot(inner[0])
        ratio_axis = figure.add_subplot(inner[1], sharex=spectrum_axis)

        leading, leading_decay = read_direct(args.leading, observable)
        lichard, lichard_decay = read_direct(args.lichard, observable)
        if len(leading) != len(lichard):
            raise ValueError(f"incompatible {observable} binning")
        if leading_decay != lichard_decay:
            raise ValueError(
                f"incompatible {observable} decay modes: "
                f"{leading_decay} versus {lichard_decay}"
            )

        x = values(leading, "bin_center")
        leading_y = values(leading, "cross_section_mb_per_gev_per_y")
        lichard_y = values(lichard, "cross_section_mb_per_gev_per_y")
        leading_error = values(leading, "standard_error_mb_per_gev_per_y")
        lichard_error = values(lichard, "standard_error_mb_per_gev_per_y")
        ratio = [numerator / denominator for numerator, denominator in zip(lichard_y, leading_y)]

        spectrum_axis.errorbar(
            x,
            leading_y,
            yerr=leading_error,
            color="#000066",
            marker="o",
            markersize=4,
            linewidth=1.8,
            capsize=0,
            label="Leading Low (reference)",
        )
        spectrum_axis.errorbar(
            x,
            lichard_y,
            yerr=lichard_error,
            color="#cc0000",
            marker="s",
            markerfacecolor="white",
            markersize=4,
            linewidth=1.8,
            capsize=0,
            label="Lichard Eq. (2.8)",
        )
        spectrum_axis.set_yscale("log")
        spectrum_axis.set_ylabel(ylabel)
        spectrum_axis.grid(True, which="both", alpha=0.18, linewidth=0.6)
        spectrum_axis.tick_params(direction="in", which="both", top=True, right=True)
        spectrum_axis.tick_params(labelbottom=False)
        spectrum_axis.legend(
            loc="lower right" if observable == "qT" else "lower left",
            frameon=False,
        )

        ratio_axis.axhline(1.0, color="#777777", linewidth=0.9, linestyle="--")
        ratio_axis.plot(x, ratio, color="#cc0000", marker="s", markersize=4, linewidth=1.8)
        ratio_axis.set_xlabel(xlabel)
        ratio_axis.text(
            0.02,
            0.08,
            "Lichard / leading",
            transform=ratio_axis.transAxes,
            ha="left",
            va="bottom",
            fontsize=9,
            color="#4b4b4b",
        )
        ratio_axis.set_xlim(*xlim)
        ratio_axis.set_ylim(0.5, 1.05)
        ratio_axis.set_yticks((0.6, 0.8, 1.0))
        ratio_axis.grid(True, alpha=0.18, linewidth=0.6)
        ratio_axis.tick_params(direction="in", top=True, right=True)

    if args.layout == "horizontal":
        figure.text(
            0.5,
            1.015,
            "LowMC non-diffractive benchmark, polarized decay, identical events and phase-space points",
            ha="center",
            va="bottom",
            fontsize=10,
            color="#000066",
        )
    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(
        args.output_prefix.with_suffix(".pdf"),
        bbox_inches="tight",
        pad_inches=0.18,
    )
    figure.savefig(
        args.output_prefix.with_suffix(".png"),
        dpi=220,
        bbox_inches="tight",
        pad_inches=0.18,
    )
    print(f"wrote {args.output_prefix.with_suffix('.pdf')}")
    print(f"wrote {args.output_prefix.with_suffix('.png')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
