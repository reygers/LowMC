#!/usr/bin/env python3
"""Dependency-free tests for the release-validation driver."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import validate_release_candidate as release


class ReleaseValidationTests(unittest.TestCase):
    def config(self, root: Path, **kwargs: object) -> release.ReleaseConfig:
        return release.ReleaseConfig(
            source_dir=Path(__file__).resolve().parents[1],
            output_dir=root,
            parallel=2,
            **kwargs,
        )

    def test_default_plan_contains_only_required_gates(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            steps = release.build_plan(self.config(Path(directory)))
        self.assertEqual(
            [step.name for step in steps],
            [
                "repository-diff-check",
                "default-configure",
                "default-build",
                "default-test-discovery",
                "default-tests",
                "default-install",
                "count-shape",
                "shape-envelope",
            ],
        )
        self.assertEqual({step.gate for step in steps}, {"R1", "R2", "R3"})

    def test_optional_plans_include_explicit_prefixes(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            config = self.config(
                Path(directory),
                enable_integrations=True,
                enable_photos=True,
                hepmc3_prefix=Path("/opt/hepmc3"),
                pythia_prefix=Path("/opt/pythia8"),
                photos_root=Path("/opt/photos"),
            )
            steps = release.build_plan(config)
        commands = {step.name: list(step.command) for step in steps}
        self.assertIn("-DCMAKE_PREFIX_PATH=/opt/hepmc3", commands["integrations-configure"])
        self.assertIn(
            "-DPYTHIA8_CONFIG_EXECUTABLE=/opt/pythia8/bin/pythia8-config",
            commands["integrations-configure"],
        )
        self.assertIn("-DPHOTOS_ROOT=/opt/photos", commands["photos-configure"])
        self.assertIn("integrations-tests", commands)
        self.assertIn("photos-tests", commands)

    def test_failure_stops_plan_and_preserves_log(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            marker = output / "should-not-exist"
            steps = [
                release.Step(
                    "expected-failure",
                    "TEST",
                    (sys.executable, "-c", "raise SystemExit(7)"),
                ),
                release.Step(
                    "must-not-run",
                    "TEST",
                    (sys.executable, "-c", f"open({str(marker)!r}, 'w').close()"),
                ),
            ]
            passed, results = release.execute_steps(
                steps,
                Path(__file__).resolve().parents[1],
                output,
                verbose=False,
            )
            self.assertFalse(passed)
            self.assertEqual(len(results), 1)
            self.assertEqual(results[0]["returncode"], 7)
            self.assertFalse(marker.exists())
            self.assertTrue((output / str(results[0]["log"])).is_file())

    def test_ctest_count_parser(self) -> None:
        payload = '{"tests": [{"name": "a"}, {"name": "b"}]}'
        self.assertEqual(release.parse_discovered_test_count(payload), 2)


if __name__ == "__main__":
    unittest.main()
