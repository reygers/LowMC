#!/usr/bin/env python3
"""Plot ALICE Fig. 4 data, the Balek baseline, and LowMC predictions."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
import sys


ALPHA = 1.0 / 137.035999084
ELECTRON_MASS_GEV = 0.00051099895000


def parse_args() -> argparse.Namespace:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--prediction",
        type=Path,
        help="CSV from lowmc_pythia8_alice_dielectron_benchmark",
    )
    parser.add_argument(
        "--data",
        type=Path,
        default=root / "data" / "alice_fig4_hepdata.csv",
        help="ALICE HEPData CSV",
    )
    parser.add_argument(
        "--output-prefix",
        type=Path,
        default=Path("alice_fig4_lowmc"),
        help="Output path without extension (writes PDF and PNG)",
    )
    parser.add_argument(
        "--sigma-inelastic",
        type=float,
        default=57.8,
        help="Balek normalization cross section in mb (default: 57.8)",
    )
    parser.add_argument(
        "--mean-charge-transfer-squared",
        type=float,
        default=1.32,
        help="Balek mean squared charge transfer",
    )
    return parser.parse_args()


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def balek_normalization(sigma_inelastic: float, charge_transfer: float) -> float:
    return sigma_inelastic * 4.0 * ALPHA**2 * charge_transfer / (3.0 * math.pi**2)


def lepton_mass_factor(mass: float) -> float:
    ratio = ELECTRON_MASS_GEV**2 / mass**2
    return (1.0 + 2.0 * ratio) * math.sqrt(max(0.0, 1.0 - 4.0 * ratio))


def balek_mass(mass: float, normalization: float, qt_max: float = 0.4) -> float:
    return (
        normalization
        * lepton_mass_factor(mass)
        * math.log1p(qt_max**2 / mass**2)
        / (2.0 * mass)
    )


def balek_qt(qt: float, normalization: float, mass_min: float = 0.15,
             mass_max: float = 0.6) -> float:
    if qt == 0.0:
        return 0.0
    qt2 = qt**2
    minimum2 = mass_min**2
    maximum2 = mass_max**2
    logarithm = math.log(
        maximum2 * (qt2 + minimum2) /
        (minimum2 * (qt2 + maximum2))
    )
    return normalization * logarithm / (2.0 * qt)


def select(rows: list[dict[str, str]], observable: str) -> list[dict[str, str]]:
    return [row for row in rows if row["observable"] == observable]


def draw_data(axis, rows: list[dict[str, str]]) -> None:
    centers = [float(row["bin_center"]) for row in rows]
    values = [float(row["value_mb_per_gev_per_y"]) for row in rows]
    stats = [float(row["stat_mb_per_gev_per_y"]) for row in rows]
    for row, center, value in zip(rows, centers, values):
        low = float(row["bin_low"])
        high = float(row["bin_high"])
        syst = float(row["syst_mb_per_gev_per_y"])
        axis.fill_between(
            [low, high],
            [max(value - syst, 1.0e-12)] * 2,
            [value + syst] * 2,
            color="black",
            alpha=0.16,
            linewidth=0,
            zorder=2,
        )
    axis.errorbar(
        centers,
        values,
        yerr=stats,
        fmt="o",
        color="black",
        markersize=4.5,
        capsize=0,
        label="ALICE excess",
        zorder=5,
    )


def draw_prediction(axis, rows: list[dict[str, str]], observable: str) -> None:
    colors = {
        "polarized": "#1f77b4",
        "isotropic": "#d62728",
        "angle_integrated": "#000066",
    }
    labels = {
        "polarized": "LowMC ND, direct (polarized)",
        "isotropic": "LowMC ND, direct (isotropic)",
        "angle_integrated": "LowMC ND, direct (angle integrated)",
    }
    markers = {"polarized": "s", "isotropic": "^", "angle_integrated": "o"}
    available_decays = {row["decay"] for row in rows}
    decays = (
        ("angle_integrated",)
        if "angle_integrated" in available_decays
        else ("polarized", "isotropic")
    )
    for decay in decays:
        selected = [
            row for row in rows
            if row["observable"] == observable and row["decay"] == decay
        ]
        sampled = sorted(
            (row for row in selected if row["method"] == "sampled"),
            key=lambda row: float(row["bin_center"]),
        )
        direct = sorted(
            (row for row in selected if row["method"] == "direct_sobol"),
            key=lambda row: float(row["bin_center"]),
        )
        if direct:
            x = [float(row["bin_center"]) for row in direct]
            y = [float(row["cross_section_mb_per_gev_per_y"]) for row in direct]
            axis.plot(x, y, color=colors[decay], linewidth=1.8, label=labels[decay])
        if sampled:
            x = [float(row["bin_center"]) for row in sampled]
            y = [float(row["cross_section_mb_per_gev_per_y"]) for row in sampled]
            error = [float(row["standard_error_mb_per_gev_per_y"]) for row in sampled]
            axis.errorbar(
                x,
                y,
                yerr=error,
                fmt=markers[decay],
                color=colors[decay],
                markerfacecolor="white",
                markersize=4,
                linewidth=1,
                label=f"LowMC ND, sampled ({decay})",
                zorder=4,
            )


def main() -> int:
    args = parse_args()
    if args.sigma_inelastic <= 0.0 or args.mean_charge_transfer_squared < 0.0:
        raise ValueError("Balek normalization inputs must be non-negative")
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is required for plotting", file=sys.stderr)
        return 2

    data_rows = read_rows(args.data)
    prediction_rows = read_rows(args.prediction) if args.prediction else []
    publication_mode = any(
        row.get("decay") == "angle_integrated" for row in prediction_rows
    )
    normalization = balek_normalization(
        args.sigma_inelastic, args.mean_charge_transfer_squared
    )

    plt.rcParams.update({
        "font.size": 10,
        "axes.labelsize": 11,
        "legend.fontsize": 8,
        "figure.dpi": 150,
    })
    figure, axes = plt.subplots(1, 2, figsize=(10.2, 4.3), constrained_layout=True)

    mass_axis, qt_axis = axes
    draw_data(mass_axis, select(data_rows, "mass"))
    mass_x = [0.15 + index * (0.60 - 0.15) / 400 for index in range(401)]
    mass_axis.plot(
        mass_x,
        [balek_mass(value, normalization) for value in mass_x],
        color="#2ca02c",
        linestyle="--",
        linewidth=2,
        label="Balek, Pišútová, and Pišút",
    )
    if prediction_rows:
        draw_prediction(mass_axis, prediction_rows, "mass")
    mass_axis.set_xlabel(r"$m_{ee}$ (GeV/$c^2$)")
    mass_axis.set_ylabel(r"$d^2\sigma/(dm_{ee}\,dy)$ (mb/(GeV/$c^2$))")
    mass_axis.set_xlim(0.14, 0.61)

    draw_data(qt_axis, select(data_rows, "qT"))
    qt_x = [0.0005 + index * (0.8 - 0.0005) / 500 for index in range(501)]
    qt_axis.plot(
        qt_x,
        [balek_qt(value, normalization) for value in qt_x],
        color="#2ca02c",
        linestyle="--",
        linewidth=2,
        label="Balek, Pišútová, and Pišút",
    )
    if prediction_rows:
        draw_prediction(qt_axis, prediction_rows, "qT")
    qt_axis.set_xlabel(r"$p_{T,ee}$ (GeV/$c$)")
    qt_axis.set_ylabel(r"$d^2\sigma/(dp_{T,ee}\,dy)$ (mb/(GeV/$c$))")
    qt_axis.set_xlim(0.0, 0.81)

    for axis in axes:
        axis.set_yscale("log")
        axis.set_ylim(1.0e-5, 0.12)
        axis.grid(True, which="both", alpha=0.18, linewidth=0.6)
        axis.tick_params(direction="in", which="both", top=True, right=True)

    acceptance_text = (
        "angle integrated, provisional $|y_{ee}|<0.8$"
        if publication_mode
        else "$|\\eta_{e^\\pm}|<0.8$, $0.075<p_{T,e^\\pm}<10$ GeV/$c$"
    )
    mass_axis.text(
        0.97,
        0.96,
        "pp, $\\sqrt{s}=13$ TeV\n$p_{T,ee}<0.4$ GeV/$c$\n"
        + acceptance_text,
        transform=mass_axis.transAxes,
        ha="right",
        va="top",
        fontsize=8,
    )
    qt_axis.text(
        0.97,
        0.96,
        "$0.15<m_{ee}<0.6$ GeV/$c^2$\nALICE luminosity uncertainty: 5%",
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
