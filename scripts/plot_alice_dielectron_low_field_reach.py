#!/usr/bin/env python3
"""Visualize how factor-two lower pair scales improve the Low-current regime."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

from plot_alice_dielectron_agreement_map import read_map


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--leading", type=Path, required=True)
    parser.add_argument("--lichard", type=Path, required=True)
    parser.add_argument("--output-prefix", type=Path, required=True)
    parser.add_argument("--max-relative-error", type=float, default=0.15)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        import matplotlib.pyplot as plt
        from matplotlib import colors
        from matplotlib.patches import FancyArrowPatch, Rectangle
        import numpy as np
    except ImportError:
        print("matplotlib and numpy are required for plotting", file=sys.stderr)
        return 2

    leading = read_map(args.leading)
    lichard = read_map(args.lichard)
    if leading.keys() != lichard.keys():
        raise ValueError("leading and Lichard files use different map bins")

    mass_edges = sorted({value for key in leading for value in key[:2]})
    qt_edges = sorted({value for key in leading for value in key[2:]})
    mass_index = {
        (mass_edges[index], mass_edges[index + 1]): index
        for index in range(len(mass_edges) - 1)
    }
    qt_index = {
        (qt_edges[index], qt_edges[index + 1]): index
        for index in range(len(qt_edges) - 1)
    }
    shape = (len(mass_edges) - 1, len(qt_edges) - 1)
    difference = np.full(shape, np.nan)
    energy_softness = np.full(shape, np.nan)
    imprecise = np.ones(shape, dtype=bool)

    for key, leading_cell in leading.items():
        lichard_cell = lichard[key]
        row = mass_index[(key[0], key[1])]
        column = qt_index[(key[2], key[3])]
        leading_yield = leading_cell["yield"]
        lichard_yield = lichard_cell["yield"]
        if leading_yield <= 0.0 or lichard_yield <= 0.0:
            continue
        difference[row, column] = abs(100.0 * (lichard_yield / leading_yield - 1.0))
        energy_softness[row, column] = leading_cell["energy_fraction_p50"]
        imprecise[row, column] = max(
            leading_cell["error"] / leading_yield,
            lichard_cell["error"] / lichard_yield,
        ) > args.max_relative_error

    displayed_difference = np.ma.masked_where(imprecise, difference)
    displayed_softness = np.ma.masked_where(imprecise, energy_softness)
    mass_centers = 0.5 * (np.asarray(mass_edges[:-1]) + np.asarray(mass_edges[1:]))
    qt_centers = 0.5 * (np.asarray(qt_edges[:-1]) + np.asarray(qt_edges[1:]))
    center_qt, center_mass = np.meshgrid(qt_centers, mass_centers)

    plt.rcParams.update({
        "font.size": 11,
        "axes.labelsize": 12,
        "figure.dpi": 150,
    })
    figure, axis = plt.subplots(figsize=(11.8, 4.7), constrained_layout=True)
    colormap = plt.get_cmap("YlGnBu").copy()
    colormap.set_bad("#d9d9d9")
    image = axis.pcolormesh(
        qt_edges,
        mass_edges,
        displayed_difference,
        cmap=colormap,
        norm=colors.Normalize(vmin=0.0, vmax=50.0),
        shading="flat",
    )
    agreement_contour = axis.contour(
        center_qt,
        center_mass,
        displayed_difference,
        levels=(5.0, 10.0),
        colors=("#cc0000", "#333333"),
        linewidths=(2.4, 1.5),
    )
    axis.clabel(
        agreement_contour,
        fmt={5.0: "5% current difference", 10.0: "10%"},
        fontsize=9,
    )
    softness_contour = axis.contour(
        center_qt,
        center_mass,
        displayed_softness,
        levels=(2.0,),
        colors=("#0486ff",),
        linestyles=("--",),
        linewidths=(2.2,),
    )
    axis.clabel(
        softness_contour,
        fmt={2.0: r"median max$(q^0/E_i)=2$"},
        fontsize=8,
    )

    run2 = Rectangle(
        (0.0, 0.15),
        0.8,
        0.45,
        fill=False,
        edgecolor="#cc0000",
        linewidth=2.0,
        linestyle="--",
        label="Run 2 comparison window",
    )
    half_scale = Rectangle(
        (0.0, 0.075),
        0.4,
        0.225,
        fill=False,
        edgecolor="#0486ff",
        linewidth=2.4,
        label=r"same window at $(p_{T,ee},m_{ee})/2$",
    )
    axis.add_patch(run2)
    axis.add_patch(half_scale)
    axis.add_patch(FancyArrowPatch(
        (0.66, 0.40),
        (0.33, 0.20),
        arrowstyle="-|>",
        mutation_scale=16,
        linewidth=2.0,
        color="#0486ff",
    ))
    axis.text(
        0.48,
        0.285,
        "factor 2 lower\nin both scales",
        color="#000066",
        fontsize=10,
        weight="bold",
        ha="center",
        va="center",
    )

    axis.set_yscale("log")
    axis.set_xlim(0.0, 0.8)
    axis.set_ylim(0.05, 0.60)
    axis.set_xlabel(r"$p_{T,ee}$ (GeV/$c$)")
    axis.set_ylabel(r"$m_{ee}$ (GeV/$c^2$)")
    axis.set_yticks((0.05, 0.075, 0.10, 0.15, 0.20, 0.30, 0.40, 0.60))
    axis.set_yticklabels(("0.05", "0.075", "0.10", "0.15", "0.20", "0.30", "0.40", "0.60"))
    axis.tick_params(direction="in", which="both", top=True, right=True)
    axis.legend(loc="lower right", frameon=True, framealpha=0.9, fontsize=9)
    colorbar = figure.colorbar(image, ax=axis, pad=0.02)
    colorbar.set_label(r"$100\,|N_{\rm Lichard}/N_{\rm leading}-1|$ (%)")
    axis.text(
        0.02,
        0.02,
        "Blue dashed contours: strict energy-softness diagnostic\n"
        "Gray: event-level relative error above 15%",
        transform=axis.transAxes,
        ha="left",
        va="bottom",
        fontsize=8,
        color="#333333",
    )
    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output_prefix.with_suffix(".pdf"))
    figure.savefig(args.output_prefix.with_suffix(".png"), dpi=220)
    print(f"wrote {args.output_prefix.with_suffix('.pdf')}")
    print(f"wrote {args.output_prefix.with_suffix('.png')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
