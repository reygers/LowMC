#!/usr/bin/env python3
"""Compare fixed-grid and weighted-sample fixed-kT diagnostics."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path


DEFAULT_ALPHA = 1.0 / 137.035999084


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compare a fixed eta-phi grid integral with the mean weighted "
            "fixed-kT sample yield."
        )
    )
    parser.add_argument("grid_csv", type=Path, help="CSV file with eta,phi,S columns")
    parser.add_argument(
        "weighted_csv",
        type=Path,
        help="CSV file with weighted fixed-kT sample rows including a weight column",
    )
    parser.add_argument(
        "--alpha",
        type=float,
        default=DEFAULT_ALPHA,
        help="fine-structure constant used to convert S to photon density",
    )
    parser.add_argument(
        "--max-abs-pull",
        type=float,
        default=None,
        help="return nonzero if |grid-weighted|/stderr exceeds this value",
    )
    return parser.parse_args()


def read_grid(path: Path) -> tuple[list[float], list[float], dict[tuple[float, float], float]]:
    values: dict[tuple[float, float], float] = {}
    etas: set[float] = set()
    phis: set[float] = set()

    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"eta", "phi", "S"}
        if set(reader.fieldnames or []) != required:
            raise ValueError("grid CSV must have exactly the columns eta,phi,S")

        for row in reader:
            eta = float(row["eta"])
            phi = float(row["phi"])
            intensity = float(row["S"])
            etas.add(eta)
            phis.add(phi)
            values[(eta, phi)] = intensity

    if not values:
        raise ValueError("grid CSV is empty")

    sorted_etas = sorted(etas)
    sorted_phis = sorted(phis)
    expected = len(sorted_etas) * len(sorted_phis)
    if len(values) != expected:
        raise ValueError("grid CSV does not contain a complete rectangular grid")
    if len(sorted_etas) < 2 or len(sorted_phis) < 2:
        raise ValueError("grid CSV needs at least two eta and phi points")

    return sorted_etas, sorted_phis, values


def require_uniform_axis(axis: list[float], name: str) -> float:
    spacing = (axis[-1] - axis[0]) / (len(axis) - 1)
    tolerance = max(1.0, abs(spacing)) * 1.0e-10

    for left, right in zip(axis, axis[1:]):
        if abs((right - left) - spacing) > tolerance:
            raise ValueError(f"{name} axis is not uniformly spaced")

    return spacing


def integrate_grid(path: Path, alpha: float) -> float:
    etas, phis, values = read_grid(path)
    deta = require_uniform_axis(etas, "eta")
    dphi = require_uniform_axis(phis, "phi")

    weighted_sum = 0.0
    for eta_index, eta in enumerate(etas):
        eta_weight = 0.5 if eta_index == 0 or eta_index == len(etas) - 1 else 1.0
        for phi_index, phi in enumerate(phis):
            phi_weight = 0.5 if phi_index == 0 or phi_index == len(phis) - 1 else 1.0
            weighted_sum += eta_weight * phi_weight * values[(eta, phi)]

    return alpha / ((2.0 * math.pi) ** 2) * weighted_sum * deta * dphi


def read_weighted_sample_weights(path: Path) -> list[float]:
    weights: list[float] = []

    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if "weight" not in (reader.fieldnames or []):
            raise ValueError("weighted sample CSV must contain a weight column")

        for row in reader:
            weight = float(row["weight"])
            if not math.isfinite(weight):
                raise ValueError("weighted sample contains a non-finite weight")
            weights.append(weight)

    if not weights:
        raise ValueError("weighted sample CSV is empty")

    return weights


def sample_mean_and_error(weights: list[float]) -> tuple[float, float]:
    count = len(weights)
    mean = sum(weights) / count

    if count < 2:
        return mean, math.nan

    variance = sum((weight - mean) ** 2 for weight in weights) / (count - 1)
    return mean, math.sqrt(variance / count)


def main() -> int:
    args = parse_args()

    try:
        grid_estimate = integrate_grid(args.grid_csv, args.alpha)
        weights = read_weighted_sample_weights(args.weighted_csv)
        sample_mean, sample_error = sample_mean_and_error(weights)
    except (OSError, ValueError) as error:
        print(f"Could not compare diagnostics: {error}", file=sys.stderr)
        return 1

    difference = sample_mean - grid_estimate
    relative_difference = difference / grid_estimate if grid_estimate != 0.0 else math.nan
    pull = difference / sample_error if sample_error > 0.0 else math.nan

    print(f"grid_integral={grid_estimate:.17g}")
    print(f"weighted_mean={sample_mean:.17g}")
    print(f"weighted_standard_error={sample_error:.17g}")
    print(f"difference={difference:.17g}")
    print(f"relative_difference={relative_difference:.17g}")
    print(f"pull={pull:.17g}")
    print(f"samples={len(weights)}")

    if args.max_abs_pull is not None and math.isfinite(pull):
        if abs(pull) > args.max_abs_pull:
            return 2

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
