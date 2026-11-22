#!/usr/bin/env python3
"""Plot matched LowMC currents against ALICE and its published baseline."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

from plot_alice_dielectron_comparison import draw_data, read_rows, select


def parse_args() -> argparse.Namespace:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lichard", type=Path, required=True)
    parser.add_argument("--leading", type=Path, required=True)
    parser.add_argument(
        "--data",
        type=Path,
        default=root / "data" / "alice_fig4_hepdata.csv",
    )
    parser.add_argument(
        "--published-bremsstrahlung",
        type=Path,
        default=root / "data" / "alice_fig4_published_bremsstrahlung.csv",
        help="vector-extracted hadronic-bremsstrahlung curve from ALICE Fig. 4",
    )
    parser.add_argument("--output-prefix", type=Path, required=True)
    return parser.parse_args()


def direct_prediction(rows: list[dict[str, str]], observable: str) -> list[dict[str, str]]:
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
            return selected
    return []


def draw_lowmc(axis, rows: list[dict[str, str]], observable: str, *,
               color: str, marker: str, linestyle: str, label: str) -> None:
    selected = direct_prediction(rows, observable)
    axis.errorbar(
        [float(row["bin_center"]) for row in selected],
        [float(row["cross_section_mb_per_gev_per_y"]) for row in selected],
        yerr=[float(row["standard_error_mb_per_gev_per_y"]) for row in selected],
        color=color,
        marker=marker,
        markerfacecolor="white" if marker == "s" else color,
        markersize=4.5,
        linewidth=2.0,
        linestyle=linestyle,
        capsize=0,
        label=label,
        zorder=4,
    )


def draw_published_bremsstrahlung(axis, rows: list[dict[str, str]],
                                  observable: str) -> None:
    selected = [row for row in rows if row["observable"] == observable]
    axis.plot(
        [float(row["x"]) for row in selected],
        [float(row["cross_section_mb_per_gev_per_y"]) for row in selected],
        color="#2ca02c",
        linestyle=":",
        linewidth=2.2,
        label="Published Run 2 bremsstrahlung curve",
    )


def main() -> int:
    args = parse_args()
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is required for plotting", file=sys.stderr)
        return 2

    data_rows = read_rows(args.data)
    published_bremsstrahlung_rows = read_rows(args.published_bremsstrahlung)
    lichard_rows = read_rows(args.lichard)
    leading_rows = read_rows(args.leading)

    plt.rcParams.update({
        "font.size": 10,
        "axes.labelsize": 11,
        "legend.fontsize": 8.5,
        "figure.dpi": 150,
    })
    figure, axes = plt.subplots(1, 2, figsize=(10.2, 4.3), constrained_layout=True)
    mass_axis, qt_axis = axes

    draw_data(mass_axis, select(data_rows, "mass"))
    draw_published_bremsstrahlung(
        mass_axis, published_bremsstrahlung_rows, "mass"
    )
    draw_lowmc(
        mass_axis,
        leading_rows,
        "mass",
        color="#000066",
        marker="o",
        linestyle="--",
        label="LowMC leading Low (soft-limit reference)",
    )
    draw_lowmc(
        mass_axis,
        lichard_rows,
        "mass",
        color="#cc0000",
        marker="s",
        linestyle="-",
        label="LowMC Lichard Eq. (2.8)",
    )
    mass_axis.set_xlabel(r"$m_{ee}$ (GeV/$c^2$)")
    mass_axis.set_ylabel(r"$d^2\sigma/(dm_{ee}\,dy)$ (mb/(GeV/$c^2$))")
    mass_axis.set_xlim(0.14, 0.61)

    draw_data(qt_axis, select(data_rows, "qT"))
    draw_published_bremsstrahlung(
        qt_axis, published_bremsstrahlung_rows, "qT"
    )
    draw_lowmc(
        qt_axis,
        leading_rows,
        "qT",
        color="#000066",
        marker="o",
        linestyle="--",
        label="LowMC leading Low (soft-limit reference)",
    )
    draw_lowmc(
        qt_axis,
        lichard_rows,
        "qT",
        color="#cc0000",
        marker="s",
        linestyle="-",
        label="LowMC Lichard Eq. (2.8)",
    )
    qt_axis.set_xlabel(r"$p_{T,ee}$ (GeV/$c$)")
    qt_axis.set_ylabel(r"$d^2\sigma/(dp_{T,ee}\,dy)$ (mb/(GeV/$c$))")
    qt_axis.set_xlim(0.0, 0.81)

    for axis in axes:
        axis.set_yscale("log")
        axis.set_ylim(1.0e-5, 0.12)
        axis.grid(True, which="both", alpha=0.18, linewidth=0.6)
        axis.tick_params(direction="in", which="both", top=True, right=True)

    mass_axis.text(
        0.97,
        0.96,
        "pp, $\\sqrt{s}=13$ TeV\n$p_{T,ee}<0.4$ GeV/$c$\n"
        "angle integrated, provisional $|y_{ee}|<0.8$",
        transform=mass_axis.transAxes,
        ha="right",
        va="top",
        fontsize=8,
    )
    qt_axis.text(
        0.97,
        0.96,
        "$0.15<m_{ee}<0.6$ GeV/$c^2$\n"
        "LowMC: non-diffractive Pythia8, angle integrated\n"
        "ALICE luminosity uncertainty: 5%",
        transform=qt_axis.transAxes,
        ha="right",
        va="top",
        fontsize=8,
    )
    mass_axis.legend(loc="lower left", frameon=False)
    qt_axis.legend(loc="lower left", frameon=False)

    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output_prefix.with_suffix(".pdf"))
    figure.savefig(args.output_prefix.with_suffix(".png"), dpi=220)
    print(f"wrote {args.output_prefix.with_suffix('.pdf')}")
    print(f"wrote {args.output_prefix.with_suffix('.png')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
