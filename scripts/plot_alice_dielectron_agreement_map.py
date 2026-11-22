#!/usr/bin/env python3
"""Map current agreement and strict energy softness for ALICE dielectrons."""

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
        "--max-relative-error",
        type=float,
        default=0.15,
        help="Mask cells whose event-level relative error exceeds this value",
    )
    parser.add_argument(
        "--agreement-threshold",
        type=float,
        default=5.0,
        help="Maximum absolute current difference in percent for the overlap",
    )
    parser.add_argument(
        "--softness-threshold",
        type=float,
        default=0.3,
        help="Maximum median max_i(q0/E_i) for the overlap",
    )
    return parser.parse_args()


def read_map(path: Path) -> dict[tuple[float, float, float, float], dict[str, float]]:
    result: dict[tuple[float, float, float, float], dict[str, float]] = {}
    with path.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            if (
                row["observable"] != "mass_qT"
                or row["method"] != "direct_sobol"
                or row["decay"] != "polarized"
            ):
                continue
            key = (
                float(row["mass_bin_low"]),
                float(row["mass_bin_high"]),
                float(row["qT_bin_low"]),
                float(row["qT_bin_high"]),
            )
            result[key] = {
                "yield": float(row["cross_section_mb_per_gev_per_y"]),
                "error": float(row["standard_error_mb_per_gev_per_y"]),
                "energy_fraction_p50": float(row["max_energy_fraction_p50"]),
                "soft_epsilon_p50": float(row["max_soft_epsilon_p50"]),
            }
    return result


def main() -> int:
    args = parse_args()
    try:
        import matplotlib.pyplot as plt
        from matplotlib import colors
        from matplotlib.patches import Rectangle
        import numpy as np
    except ImportError:
        print("matplotlib and numpy are required for plotting", file=sys.stderr)
        return 2

    leading = read_map(args.leading)
    lichard = read_map(args.lichard)
    if leading.keys() != lichard.keys():
        raise ValueError("leading and Lichard files use different map bins")

    mass_edges = sorted({value for key in leading for value in key[:2]})
    qT_edges = sorted({value for key in leading for value in key[2:]})
    mass_index = {(mass_edges[i], mass_edges[i + 1]): i
                  for i in range(len(mass_edges) - 1)}
    qT_index = {(qT_edges[i], qT_edges[i + 1]): i
                for i in range(len(qT_edges) - 1)}
    shape = (len(mass_edges) - 1, len(qT_edges) - 1)
    signed_difference = np.full(shape, np.nan)
    energy_fraction = np.full(shape, np.nan)
    masked_for_precision = np.ones(shape, dtype=bool)

    for key, leading_cell in leading.items():
        lichard_cell = lichard[key]
        i = mass_index[(key[0], key[1])]
        j = qT_index[(key[2], key[3])]
        leading_yield = leading_cell["yield"]
        lichard_yield = lichard_cell["yield"]
        if leading_yield <= 0.0 or lichard_yield <= 0.0:
            continue
        leading_relative_error = leading_cell["error"] / leading_yield
        lichard_relative_error = lichard_cell["error"] / lichard_yield
        signed_difference[i, j] = 100.0 * (lichard_yield / leading_yield - 1.0)
        energy_fraction[i, j] = leading_cell["energy_fraction_p50"]
        masked_for_precision[i, j] = max(
            leading_relative_error, lichard_relative_error
        ) > args.max_relative_error

    display = np.ma.masked_where(masked_for_precision, signed_difference)
    mass_centers = 0.5 * (np.asarray(mass_edges[:-1]) + np.asarray(mass_edges[1:]))
    qT_centers = 0.5 * (np.asarray(qT_edges[:-1]) + np.asarray(qT_edges[1:]))
    center_qT, center_mass = np.meshgrid(qT_centers, mass_centers)

    plt.rcParams.update({
        "font.size": 11,
        "axes.labelsize": 12,
        "figure.dpi": 150,
    })
    figure, axes = plt.subplots(1, 2, figsize=(11.6, 4.8), constrained_layout=True)
    agreement_axis, softness_axis = axes
    agreement = np.ma.abs(display)
    agreement_colormap = plt.get_cmap("YlGnBu").copy()
    agreement_colormap.set_bad("#d9d9d9")
    agreement_image = agreement_axis.pcolormesh(
        qT_edges,
        mass_edges,
        agreement,
        cmap=agreement_colormap,
        norm=colors.Normalize(vmin=0.0, vmax=50.0),
        shading="flat",
    )
    contours = agreement_axis.contour(
        center_qT,
        center_mass,
        agreement,
        levels=(3.0, 5.0, 10.0),
        colors=("#0486ff", "#cc0000", "#333333"),
        linewidths=(1.8, 2.2, 1.4),
    )
    agreement_axis.clabel(
        contours, fmt={3.0: "3%", 5.0: "5%", 10.0: "10%"}, fontsize=8
    )
    agreement_colorbar = figure.colorbar(agreement_image, ax=agreement_axis, pad=0.02)
    agreement_colorbar.set_label(
        r"$100\,|N_{\rm Lichard}/N_{\rm leading}-1|$ (%)"
    )
    agreement_axis.set_title("A  Current-prescription agreement", loc="left")
    agreement_axis.text(
        0.02,
        0.02,
        "Contours: 3%, 5%, 10%\nGray: relative error above 15%",
        transform=agreement_axis.transAxes,
        ha="left",
        va="bottom",
        fontsize=8,
        color="#333333",
    )

    softness = np.ma.masked_where(masked_for_precision, energy_fraction)
    softness_colormap = plt.get_cmap("magma").copy()
    softness_colormap.set_bad("#d9d9d9")
    softness_image = softness_axis.pcolormesh(
        qT_edges,
        mass_edges,
        softness,
        cmap=softness_colormap,
        norm=colors.LogNorm(vmin=0.3, vmax=20.0),
        shading="flat",
    )
    softness_colorbar = figure.colorbar(softness_image, ax=softness_axis, pad=0.02)
    softness_colorbar.set_label(
        r"yield-weighted median of $\max_i(q^0/E_i)$"
    )
    softness_axis.set_title("B  Strict energy-softness diagnostic", loc="left")

    overlap = (
        (~masked_for_precision)
        & np.isfinite(signed_difference)
        & np.isfinite(energy_fraction)
        & (np.abs(signed_difference) < args.agreement_threshold)
        & (energy_fraction < args.softness_threshold)
    )
    if np.any(overlap):
        softness_axis.contourf(
            center_qT,
            center_mass,
            overlap.astype(float),
            levels=(0.5, 1.5),
            colors="none",
            hatches=("////",),
        )
        overlap_note = (
            f"Hatched: both <{args.agreement_threshold:g}% and "
            f"<{args.softness_threshold:g}"
        )
    else:
        overlap_note = (
            f"No bin satisfies both <{args.agreement_threshold:g}% and "
            f"median max(q0/Ei) < {args.softness_threshold:g}"
        )
    softness_axis.text(
        0.02,
        0.02,
        overlap_note,
        transform=softness_axis.transAxes,
        ha="left",
        va="bottom",
        fontsize=8,
        color="#ffffff",
        bbox={"facecolor": "#000066", "edgecolor": "none", "alpha": 0.82},
    )

    for axis in axes:
        axis.add_patch(Rectangle(
            (0.0, 0.15),
            0.8,
            0.45,
            fill=False,
            edgecolor="#cc0000",
            linewidth=1.8,
            linestyle="--",
            label="Run 2 window",
        ))
        axis.set_yscale("log")
        axis.set_xlim(0.0, 0.8)
        axis.set_ylim(mass_edges[0], mass_edges[-1])
        axis.set_xlabel(r"$p_{T,ee}$ (GeV/$c$)")
        axis.set_yticks((0.02, 0.03, 0.05, 0.10, 0.20, 0.40, 0.60))
        axis.set_yticklabels(
            ("0.02", "0.03", "0.05", "0.10", "0.20", "0.40", "0.60")
        )
        axis.tick_params(direction="in", which="both", top=True, right=True)
        axis.legend(loc="upper right", frameon=False, fontsize=8)
    agreement_axis.set_ylabel(r"$m_{ee}$ (GeV/$c^2$)")
    softness_axis.set_ylabel(r"$m_{ee}$ (GeV/$c^2$)")
    figure.suptitle(
        "Current agreement is necessary—but not sufficient",
        color="#000066",
        fontsize=13,
    )

    finite_energy = energy_fraction[~np.isnan(energy_fraction)]
    if finite_energy.size > 0:
        print(
            "median_max_energy_fraction_range="
            f"{finite_energy.min():.6g},{finite_energy.max():.6g}"
        )
    print(f"overlap_cells={int(np.count_nonzero(overlap))}")
    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output_prefix.with_suffix(".pdf"))
    figure.savefig(args.output_prefix.with_suffix(".png"), dpi=220)
    print(f"wrote {args.output_prefix.with_suffix('.pdf')}")
    print(f"wrote {args.output_prefix.with_suffix('.png')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
