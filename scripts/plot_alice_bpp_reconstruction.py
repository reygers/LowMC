#!/usr/bin/env python3
"""Compare the published ALICE bremsstrahlung curve with BPP Eq. (9a)."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

from plot_alice_dielectron_comparison import (
    balek_mass,
    balek_normalization,
    balek_qt,
    read_rows,
)


def parse_args() -> argparse.Namespace:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--published-bremsstrahlung",
        type=Path,
        default=root / "data" / "alice_fig4_published_bremsstrahlung.csv",
    )
    parser.add_argument(
        "--output-prefix",
        type=Path,
        default=Path("alice_bpp_published_vs_eq9a"),
    )
    parser.add_argument("--sigma-inelastic", type=float, default=57.8)
    parser.add_argument("--mean-charge-transfer-squared", type=float, default=1.32)
    return parser.parse_args()


def values(rows: list[dict[str, str]], observable: str) -> tuple[list[float], list[float]]:
    selected = [row for row in rows if row["observable"] == observable]
    return (
        [float(row["x"]) for row in selected],
        [float(row["cross_section_mb_per_gev_per_y"]) for row in selected],
    )


def main() -> int:
    args = parse_args()
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is required for plotting", file=sys.stderr)
        return 2

    rows = read_rows(args.published_bremsstrahlung)
    mass_x, mass_published = values(rows, "mass")
    qt_x, qt_published = values(rows, "qT")
    normalization = balek_normalization(
        args.sigma_inelastic, args.mean_charge_transfer_squared
    )
    mass_reconstruction = [balek_mass(x, normalization) for x in mass_x]
    qt_reconstruction = [balek_qt(x, normalization) for x in qt_x]
    mass_ratio = [a / b for a, b in zip(mass_published, mass_reconstruction)]
    qt_ratio = [a / b for a, b in zip(qt_published, qt_reconstruction)]

    plt.rcParams.update({
        "font.size": 11,
        "axes.labelsize": 12,
        "legend.fontsize": 9.5,
        "figure.dpi": 150,
    })
    figure = plt.figure(figsize=(10.0, 4.7), constrained_layout=True)
    grid = figure.add_gridspec(2, 2, height_ratios=(3.0, 1.25))
    mass_axis = figure.add_subplot(grid[0, 0])
    mass_ratio_axis = figure.add_subplot(grid[1, 0], sharex=mass_axis)
    qt_axis = figure.add_subplot(grid[0, 1])
    qt_ratio_axis = figure.add_subplot(grid[1, 1], sharex=qt_axis)

    published_style = {
        "color": "#2ca02c",
        "linewidth": 2.4,
        "label": "Published Run 2 curve",
    }
    reconstructed_style = {
        "color": "#000066",
        "linewidth": 2.0,
        "linestyle": "--",
        "label": "BPP Eq. (9a) reconstruction",
    }
    mass_axis.plot(mass_x, mass_published, **published_style)
    mass_axis.plot(mass_x, mass_reconstruction, **reconstructed_style)
    qt_axis.plot(qt_x, qt_published, **published_style)
    qt_axis.plot(qt_x, qt_reconstruction, **reconstructed_style)

    mass_ratio_axis.plot(mass_x, mass_ratio, color="#2ca02c", linewidth=2.2)
    qt_ratio_axis.plot(qt_x, qt_ratio, color="#2ca02c", linewidth=2.2)
    for axis in (mass_ratio_axis, qt_ratio_axis):
        axis.axhline(1.0, color="#777777", linewidth=1.0, linestyle=":")
        axis.set_ylabel("published /\nEq. (9a)")
        axis.grid(True, alpha=0.18, linewidth=0.6)
        axis.tick_params(direction="in", top=True, right=True)

    mass_axis.set_yscale("log")
    qt_axis.set_yscale("log")
    mass_axis.set_ylabel(r"$d^2\sigma/(dm_{ee}\,dy)$ (mb/(GeV/$c^2$))")
    qt_axis.set_ylabel(r"$d^2\sigma/(dp_{T,ee}\,dy)$ (mb/(GeV/$c$))")
    mass_ratio_axis.set_xlabel(r"$m_{ee}$ (GeV/$c^2$)")
    qt_ratio_axis.set_xlabel(r"$p_{T,ee}$ (GeV/$c$)")
    mass_ratio_axis.set_ylim(1.18, 1.55)
    qt_ratio_axis.set_ylim(0.60, 1.45)
    mass_axis.set_xlim(0.15, 0.60)
    qt_axis.set_xlim(0.0, 0.80)

    for axis in (mass_axis, qt_axis):
        axis.grid(True, which="both", alpha=0.18, linewidth=0.6)
        axis.tick_params(direction="in", which="both", top=True, right=True)
        axis.legend(frameon=False, loc="best")

    mass_axis.text(
        0.97,
        0.18,
        r"$p_{T,ee}<0.4$ GeV/$c$",
        transform=mass_axis.transAxes,
        ha="right",
        va="bottom",
        fontsize=10,
    )
    qt_axis.text(
        0.97,
        0.06,
        r"$0.15<m_{ee}<0.6$ GeV/$c^2$",
        transform=qt_axis.transAxes,
        ha="right",
        va="bottom",
        fontsize=10,
    )

    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output_prefix.with_suffix(".pdf"))
    figure.savefig(args.output_prefix.with_suffix(".png"), dpi=220)
    print(f"wrote {args.output_prefix.with_suffix('.pdf')}")
    print(f"wrote {args.output_prefix.with_suffix('.png')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
