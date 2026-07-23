#!/usr/bin/env python3
"""Run the LowMC release-validation matrix and preserve an evidence manifest."""

from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import platform
import shlex
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Callable


@dataclass(frozen=True)
class ReleaseConfig:
    source_dir: Path
    output_dir: Path
    parallel: int
    enable_integrations: bool = False
    enable_photos: bool = False
    hepmc3_prefix: Path | None = None
    pythia_prefix: Path | None = None
    photos_root: Path | None = None


@dataclass(frozen=True)
class Step:
    name: str
    gate: str
    command: tuple[str, ...]
    discover_tests: bool = False


def parse_args() -> argparse.Namespace:
    default_source = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(
        description=(
            "Run the release-blocking LowMC validation gates in fresh build "
            "directories and write JSON and Markdown evidence."
        )
    )
    parser.add_argument(
        "--source-dir",
        type=Path,
        default=default_source,
        help="LowMC source tree; defaults to the parent of this script",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        help=(
            "fresh run directory; defaults to "
            "build/release-validation/<UTC timestamp>"
        ),
    )
    parser.add_argument(
        "--parallel",
        type=int,
        default=max(1, min(os.cpu_count() or 1, 8)),
        help="parallel build jobs",
    )
    parser.add_argument(
        "--enable-integrations",
        action="store_true",
        help="run the combined HepMC3 and Pythia8 conditional release gate",
    )
    parser.add_argument(
        "--enable-photos",
        action="store_true",
        help="build and test the optional PHOTOS benchmark",
    )
    parser.add_argument("--hepmc3-prefix", type=Path)
    parser.add_argument("--pythia-prefix", type=Path)
    parser.add_argument("--photos-root", type=Path)
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="write the planned manifest and summary without executing commands",
    )
    args = parser.parse_args()

    if args.parallel < 1:
        parser.error("--parallel must be positive")
    if args.enable_integrations:
        if args.hepmc3_prefix is None:
            parser.error("--enable-integrations requires --hepmc3-prefix")
        if args.pythia_prefix is None:
            parser.error("--enable-integrations requires --pythia-prefix")
    if args.enable_photos:
        if args.hepmc3_prefix is None:
            parser.error("--enable-photos requires --hepmc3-prefix")
        if args.photos_root is None:
            parser.error("--enable-photos requires --photos-root")
    return args


def timestamp_utc() -> str:
    return dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def timestamp_path_utc() -> str:
    return dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def prepare_output_dir(path: Path) -> None:
    if path.exists():
        if not path.is_dir():
            raise FileExistsError(f"release-validation output path is not a directory: {path}")
        if any(path.iterdir()):
            raise FileExistsError(f"release-validation output directory is not empty: {path}")
    path.mkdir(parents=True, exist_ok=True)
    (path / "logs").mkdir()
    (path / "evidence").mkdir()


def build_plan(config: ReleaseConfig) -> list[Step]:
    source = config.source_dir
    output = config.output_dir
    default_build = output / "build-default"
    default_install = output / "install-default"
    evidence = output / "evidence"
    python = sys.executable
    parallel = str(config.parallel)

    steps = [
        Step("repository-diff-check", "R1", ("git", "diff", "--check")),
        Step(
            "default-configure",
            "R1",
            (
                "cmake",
                "-S",
                str(source),
                "-B",
                str(default_build),
                "-DCMAKE_BUILD_TYPE=Release",
            ),
        ),
        Step(
            "default-build",
            "R1",
            ("cmake", "--build", str(default_build), "--parallel", parallel),
        ),
        Step(
            "default-test-discovery",
            "R1",
            ("ctest", "--test-dir", str(default_build), "--show-only=json-v1"),
            discover_tests=True,
        ),
        Step(
            "default-tests",
            "R1",
            ("ctest", "--test-dir", str(default_build), "--output-on-failure"),
        ),
        Step(
            "default-install",
            "R1",
            ("cmake", "--install", str(default_build), "--prefix", str(default_install)),
        ),
        Step(
            "count-shape",
            "R2",
            (
                python,
                str(source / "scripts" / "validate_count_shape.py"),
                str(evidence / "count-shape.csv"),
                "--aggregate-csv",
                str(evidence / "count-shape-aggregate.csv"),
                "--summary-md",
                str(evidence / "count-shape.md"),
                "--output-dir",
                str(evidence / "count-shape-runs"),
                "--diagnostic-exe",
                str(default_build / "lowmc_count_shape_diagnostic"),
                "--seeds",
                "3",
                "--seed",
                "12345",
                "--oversampling-factor",
                "100",
            ),
        ),
        Step(
            "shape-envelope",
            "R3",
            (
                python,
                str(source / "scripts" / "validate_shape_unweighting_candidates.py"),
                str(evidence / "shape-envelope.csv"),
                "--summary-md",
                str(evidence / "shape-envelope.md"),
                "--work-dir",
                str(evidence / "shape-envelope-runs"),
                "--truth-exe",
                str(default_build / "lowmc_two_particle_diagnostic"),
                "--shape-exe",
                str(default_build / "lowmc_shape_unweighting_diagnostic"),
                "--candidate-set",
                "targeted",
                "--trials",
                "30000",
                "--seed",
                "12345",
                "--seed",
                "23456",
                "--seed",
                "34567",
                "--gate-policy",
                "strict",
                "--gate-max-trial-bound-ratio",
                "0.8",
            ),
        ),
    ]

    if config.enable_integrations:
        if config.hepmc3_prefix is None or config.pythia_prefix is None:
            raise ValueError("integration plan requires HepMC3 and Pythia8 prefixes")
        build = output / "build-integrations"
        install = output / "install-integrations"
        steps.extend(
            [
                Step(
                    "integrations-configure",
                    "C1",
                    (
                        "cmake",
                        "-S",
                        str(source),
                        "-B",
                        str(build),
                        "-DCMAKE_BUILD_TYPE=Release",
                        "-DLOWMC_ENABLE_HEPMC3=ON",
                        "-DLOWMC_ENABLE_PYTHIA8=ON",
                        f"-DCMAKE_PREFIX_PATH={config.hepmc3_prefix}",
                        f"-DPYTHIA8_CONFIG_EXECUTABLE={config.pythia_prefix / 'bin' / 'pythia8-config'}",
                    ),
                ),
                Step(
                    "integrations-build",
                    "C1",
                    ("cmake", "--build", str(build), "--parallel", parallel),
                ),
                Step(
                    "integrations-test-discovery",
                    "C1",
                    ("ctest", "--test-dir", str(build), "--show-only=json-v1"),
                    discover_tests=True,
                ),
                Step(
                    "integrations-tests",
                    "C1",
                    ("ctest", "--test-dir", str(build), "--output-on-failure"),
                ),
                Step(
                    "integrations-install",
                    "C1",
                    ("cmake", "--install", str(build), "--prefix", str(install)),
                ),
            ]
        )

    if config.enable_photos:
        if config.hepmc3_prefix is None or config.photos_root is None:
            raise ValueError("PHOTOS plan requires HepMC3 and PHOTOS prefixes")
        build = output / "build-photos"
        steps.extend(
            [
                Step(
                    "photos-configure",
                    "PHOTOS",
                    (
                        "cmake",
                        "-S",
                        str(source),
                        "-B",
                        str(build),
                        "-DCMAKE_BUILD_TYPE=Release",
                        "-DLOWMC_ENABLE_PHOTOS=ON",
                        f"-DPHOTOS_ROOT={config.photos_root}",
                        f"-DCMAKE_PREFIX_PATH={config.hepmc3_prefix}",
                    ),
                ),
                Step(
                    "photos-build",
                    "PHOTOS",
                    ("cmake", "--build", str(build), "--parallel", parallel),
                ),
                Step(
                    "photos-test-discovery",
                    "PHOTOS",
                    ("ctest", "--test-dir", str(build), "--show-only=json-v1"),
                    discover_tests=True,
                ),
                Step(
                    "photos-tests",
                    "PHOTOS",
                    ("ctest", "--test-dir", str(build), "--output-on-failure"),
                ),
            ]
        )

    return steps


def probe(command: list[str], cwd: Path) -> dict[str, object]:
    try:
        result = subprocess.run(command, cwd=cwd, capture_output=True, text=True)
    except OSError as error:
        return {"command": command, "returncode": None, "output": str(error)}
    output = (result.stdout or result.stderr).strip()
    return {
        "command": command,
        "returncode": result.returncode,
        "output": output,
    }


def first_line(probe_result: dict[str, object]) -> str | None:
    output = str(probe_result.get("output", "")).strip()
    return output.splitlines()[0] if output else None


def dependency_metadata(config: ReleaseConfig) -> dict[str, object]:
    result: dict[str, object] = {
        "hepmc3": {
            "enabled": config.enable_integrations or config.enable_photos,
            "prefix": str(config.hepmc3_prefix) if config.hepmc3_prefix else None,
            "version": None,
        },
        "pythia8": {
            "enabled": config.enable_integrations,
            "prefix": str(config.pythia_prefix) if config.pythia_prefix else None,
            "version": None,
        },
        "photos": {
            "enabled": config.enable_photos,
            "prefix": str(config.photos_root) if config.photos_root else None,
            "version": None,
        },
    }
    if config.hepmc3_prefix:
        version = probe(
            [str(config.hepmc3_prefix / "bin" / "HepMC3-config"), "--version"],
            config.source_dir,
        )
        result["hepmc3"]["version"] = first_line(version)  # type: ignore[index]
    if config.pythia_prefix:
        version = probe(
            [str(config.pythia_prefix / "bin" / "pythia8-config"), "--version"],
            config.source_dir,
        )
        result["pythia8"]["version"] = first_line(version)  # type: ignore[index]
    if config.photos_root:
        candidates = sorted((config.photos_root / "lib").glob("libPhotospp.*"))
        result["photos"]["libraries"] = [str(path) for path in candidates]  # type: ignore[index]
    return result


def collect_metadata(config: ReleaseConfig, steps: list[Step]) -> dict[str, object]:
    revision = probe(["git", "rev-parse", "HEAD"], config.source_dir)
    status = probe(["git", "status", "--porcelain"], config.source_dir)
    cmake = probe(["cmake", "--version"], config.source_dir)
    compiler_command = os.environ.get("CXX", "c++")
    compiler = probe([compiler_command, "--version"], config.source_dir)
    return {
        "schema_version": 1,
        "status": "planned",
        "started_at_utc": timestamp_utc(),
        "finished_at_utc": None,
        "source_dir": str(config.source_dir),
        "output_dir": str(config.output_dir),
        "revision": first_line(revision),
        "dirty_tree": bool(str(status.get("output", "")).strip()),
        "dirty_paths": str(status.get("output", "")).splitlines(),
        "platform": platform.platform(),
        "tools": {
            "python": sys.version.splitlines()[0],
            "cmake": first_line(cmake),
            "compiler_command": compiler_command,
            "compiler": first_line(compiler),
        },
        "configuration": {
            "build_type": "Release",
            "parallel": config.parallel,
            "enable_integrations": config.enable_integrations,
            "enable_photos": config.enable_photos,
        },
        "dependencies": dependency_metadata(config),
        "planned_steps": [
            {
                "name": step.name,
                "gate": step.gate,
                "command": list(step.command),
            }
            for step in steps
        ],
        "steps": [],
    }


def parse_discovered_test_count(stdout: str) -> int:
    payload = json.loads(stdout)
    tests = payload.get("tests")
    if not isinstance(tests, list):
        raise ValueError("CTest JSON output does not contain a tests list")
    return len(tests)


def execute_steps(
    steps: list[Step],
    source_dir: Path,
    output_dir: Path,
    environment: dict[str, str] | None = None,
    on_result: Callable[[dict[str, object]], None] | None = None,
    verbose: bool = True,
) -> tuple[bool, list[dict[str, object]]]:
    logs = output_dir / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    results: list[dict[str, object]] = []
    env = os.environ.copy()
    if environment:
        env.update(environment)

    for index, step in enumerate(steps, start=1):
        log_path = logs / f"{index:02d}-{step.name}.log"
        command = list(step.command)
        if verbose:
            print(f"[{step.gate}] {step.name}: {shlex.join(command)}", flush=True)
        started = time.perf_counter()
        try:
            process = subprocess.run(
                command,
                cwd=source_dir,
                env=env,
                capture_output=True,
                text=True,
            )
            returncode: int | None = process.returncode
            stdout = process.stdout
            stderr = process.stderr
            error = None
        except OSError as caught:
            returncode = None
            stdout = ""
            stderr = ""
            error = str(caught)
        duration = time.perf_counter() - started
        log_text = f"$ {shlex.join(command)}\n\n[stdout]\n{stdout}\n[stderr]\n{stderr}"
        if error:
            log_text += f"\n[execution error]\n{error}\n"
        log_path.write_text(log_text, encoding="utf-8")

        result: dict[str, object] = {
            "name": step.name,
            "gate": step.gate,
            "command": command,
            "status": "passed" if returncode == 0 else "failed",
            "returncode": returncode,
            "duration_seconds": duration,
            "log": str(log_path.relative_to(output_dir)),
        }
        if error:
            result["error"] = error
        if returncode == 0 and step.discover_tests:
            try:
                result["discovered_test_count"] = parse_discovered_test_count(stdout)
            except (json.JSONDecodeError, ValueError) as caught:
                result["status"] = "failed"
                result["error"] = f"could not parse CTest discovery output: {caught}"

        results.append(result)
        if on_result:
            on_result(result)
        if result["status"] != "passed":
            if verbose:
                print(f"FAILED: {step.name}; see {log_path}", file=sys.stderr)
            return False, results
        if verbose:
            print(f"passed in {duration:.2f} s", flush=True)

    return True, results


def write_manifest(path: Path, manifest: dict[str, object]) -> None:
    path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


def write_summary(path: Path, manifest: dict[str, object]) -> None:
    configuration = manifest["configuration"]
    dependencies = manifest["dependencies"]
    steps = manifest["steps"]
    lines = [
        "# LowMC Release Validation",
        "",
        f"- Status: `{manifest['status']}`",
        f"- Revision: `{manifest.get('revision') or 'unknown'}`",
        f"- Dirty source tree: `{'yes' if manifest['dirty_tree'] else 'no'}`",
        f"- Started: `{manifest['started_at_utc']}`",
        f"- Finished: `{manifest.get('finished_at_utc') or 'not finished'}`",
        f"- Integrations enabled: `{configuration['enable_integrations']}`",
        f"- PHOTOS enabled: `{configuration['enable_photos']}`",
        "",
        "## Dependencies",
        "",
        "| Component | Enabled | Prefix | Version |",
        "|---|---:|---|---|",
    ]
    for key, label in (("hepmc3", "HepMC3"), ("pythia8", "Pythia8"), ("photos", "PHOTOS")):
        item = dependencies[key]
        lines.append(
            f"| {label} | {item['enabled']} | {item.get('prefix') or '-'} | "
            f"{item.get('version') or '-'} |"
        )
    lines.extend(
        [
            "",
            "## Steps",
            "",
            "| Gate | Step | Status | Seconds | Tests | Log |",
            "|---|---|---|---:|---:|---|",
        ]
    )
    for item in steps:
        tests = item.get("discovered_test_count", "")
        lines.append(
            f"| {item['gate']} | {item['name']} | {item['status']} | "
            f"{item['duration_seconds']:.2f} | {tests} | `{item['log']}` |"
        )
    if manifest["status"] == "planned":
        lines.extend(
            [
                "",
                "Dry run only. Commands are recorded in `manifest.json` and were not executed.",
            ]
        )
    lines.extend(
        [
            "",
            "DELPHI and PHOTOS production-statistics studies are deliberately not run by "
            "this driver. Their trigger conditions remain in "
            "`docs/RELEASE_VALIDATION_MATRIX.md`.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    args = parse_args()
    source_dir = args.source_dir.resolve()
    if not (source_dir / "CMakeLists.txt").is_file():
        print(f"not a LowMC source tree: {source_dir}", file=sys.stderr)
        return 2
    output_dir = (
        args.output_dir.resolve()
        if args.output_dir
        else source_dir / "build" / "release-validation" / timestamp_path_utc()
    )
    config = ReleaseConfig(
        source_dir=source_dir,
        output_dir=output_dir,
        parallel=args.parallel,
        enable_integrations=args.enable_integrations,
        enable_photos=args.enable_photos,
        hepmc3_prefix=args.hepmc3_prefix.resolve() if args.hepmc3_prefix else None,
        pythia_prefix=args.pythia_prefix.resolve() if args.pythia_prefix else None,
        photos_root=args.photos_root.resolve() if args.photos_root else None,
    )

    try:
        prepare_output_dir(output_dir)
        steps = build_plan(config)
    except (FileExistsError, ValueError) as error:
        print(error, file=sys.stderr)
        return 2

    manifest = collect_metadata(config, steps)
    manifest_path = output_dir / "manifest.json"
    summary_path = output_dir / "summary.md"
    write_manifest(manifest_path, manifest)
    write_summary(summary_path, manifest)

    if args.dry_run:
        print(f"planned {len(steps)} steps")
        print(f"manifest: {manifest_path}")
        print(f"summary: {summary_path}")
        return 0

    def update_manifest(result: dict[str, object]) -> None:
        manifest["steps"].append(result)  # type: ignore[union-attr]
        manifest["status"] = "running"
        write_manifest(manifest_path, manifest)
        write_summary(summary_path, manifest)

    environment: dict[str, str] = {}
    if config.pythia_prefix:
        environment["PATH"] = (
            str(config.pythia_prefix / "bin") + os.pathsep + os.environ.get("PATH", "")
        )
    passed, _ = execute_steps(
        steps,
        source_dir,
        output_dir,
        environment=environment,
        on_result=update_manifest,
    )
    manifest["status"] = "passed" if passed else "failed"
    manifest["finished_at_utc"] = timestamp_utc()
    write_manifest(manifest_path, manifest)
    write_summary(summary_path, manifest)
    print(f"manifest: {manifest_path}")
    print(f"summary: {summary_path}")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
