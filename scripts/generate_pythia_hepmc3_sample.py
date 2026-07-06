#!/usr/bin/env python3
"""Generate a small optional Pythia-produced HepMC3 sample for LowMC validation."""

from __future__ import annotations

import argparse
import shlex
import subprocess
import sys
from pathlib import Path


GENERATOR_SOURCE = r"""
#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC3.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/WriterAscii.h"

#include <cstdlib>
#include <iostream>
#include <string>

int main(int argc, char** argv) {
  if (argc != 5) {
    std::cerr << "Usage: lowmc_pythia_hepmc3_sample <output.hepmc3> <events> <seed> <eCMGeV>\n";
    return EXIT_FAILURE;
  }

  const std::string outputPath = argv[1];
  const int requestedEvents = std::stoi(argv[2]);
  const int seed = std::stoi(argv[3]);
  const double eCM = std::stod(argv[4]);
  if (requestedEvents <= 0 || seed <= 0 || eCM <= 0.0) {
    std::cerr << "FAILED: events, seed, and eCM must be positive\n";
    return EXIT_FAILURE;
  }

  Pythia8::Pythia pythia;
  pythia.readString("Print:quiet = on");
  pythia.readString("Beams:idA = 2212");
  pythia.readString("Beams:idB = 2212");
  pythia.readString("Beams:eCM = " + std::to_string(eCM));
  pythia.readString("SoftQCD:nonDiffractive = on");
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = " + std::to_string(seed));
  if (!pythia.init()) {
    std::cerr << "FAILED: Pythia initialization failed\n";
    return EXIT_FAILURE;
  }

  HepMC3::WriterAscii writer(outputPath);
  if (writer.failed()) {
    std::cerr << "FAILED: could not open HepMC3 output file: " << outputPath << "\n";
    return EXIT_FAILURE;
  }
  HepMC3::Pythia8ToHepMC3 converter;

  int attemptedEvents = 0;
  int writtenEvents = 0;
  while (writtenEvents < requestedEvents) {
    ++attemptedEvents;
    if (!pythia.next()) {
      if (attemptedEvents > requestedEvents * 20) {
        std::cerr << "FAILED: too many Pythia event-generation failures\n";
        return EXIT_FAILURE;
      }
      continue;
    }

    HepMC3::GenEvent event;
    if (!converter.fill_next_event(pythia, &event, writtenEvents)) {
      std::cerr << "FAILED: Pythia-to-HepMC3 conversion failed\n";
      return EXIT_FAILURE;
    }
    writer.write_event(event);
    ++writtenEvents;
  }

  pythia.stat();
  std::cout << "pythia_hepmc3_sample events=" << writtenEvents
            << " attempted_events=" << attemptedEvents
            << " seed=" << seed
            << " eCM_GeV=" << eCM
            << " output=" << outputPath << "\n";
  return EXIT_SUCCESS;
}
"""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compile and run a tiny local Pythia8 -> HepMC3 generator. "
            "This is an optional/manual validation helper and does not add a "
            "Pythia dependency to LowMC."
        )
    )
    parser.add_argument("--pythia-prefix", type=Path, default=Path.home() / "pythia8317")
    parser.add_argument("--hepmc3-prefix", type=Path, default=Path.home() / "hepmc3")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, default=Path("build/pythia_hepmc3_sample"))
    parser.add_argument("--events", type=int, default=10)
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--ecm", type=float, default=13000.0)
    parser.add_argument("--cxx", default="c++")
    parser.add_argument(
        "--keep-source",
        action="store_true",
        help="keep the generated C++ source in the work directory",
    )
    return parser.parse_args()


def run_command(command: list[str]) -> str:
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        if result.stdout:
            print(result.stdout, file=sys.stdout, end="")
        if result.stderr:
            print(result.stderr, file=sys.stderr, end="")
        raise subprocess.CalledProcessError(result.returncode, command)
    if result.stderr:
        print(result.stderr, file=sys.stderr, end="")
    return result.stdout


def require_file(path: Path, description: str) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"{description} not found: {path}")


def split_flags(text: str) -> list[str]:
    return shlex.split(text.strip())


def main() -> int:
    args = parse_args()
    pythia_config = args.pythia_prefix / "bin" / "pythia8-config"
    hepmc3_config = args.hepmc3_prefix / "bin" / "HepMC3-config"
    require_file(pythia_config, "pythia8-config")
    require_file(hepmc3_config, "HepMC3-config")
    if args.events <= 0:
        raise ValueError("--events must be positive")
    if args.seed <= 0:
        raise ValueError("--seed must be positive")
    if args.ecm <= 0.0:
        raise ValueError("--ecm must be positive")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.work_dir.mkdir(parents=True, exist_ok=True)
    source_path = args.work_dir / "lowmc_pythia_hepmc3_sample.cpp"
    executable_path = args.work_dir / "lowmc_pythia_hepmc3_sample"
    source_path.write_text(GENERATOR_SOURCE, encoding="utf-8")

    pythia_flags = split_flags(run_command([str(pythia_config), "--cxxflags", "--libs"]))
    hepmc3_flags = split_flags(run_command([str(hepmc3_config), "--cflags", "--libs"]))
    hepmc3_libdir = run_command([str(hepmc3_config), "--libdir"]).strip()
    compile_command = [
        args.cxx,
        "-std=c++17",
        str(source_path),
        "-o",
        str(executable_path),
        *pythia_flags,
        *hepmc3_flags,
        f"-Wl,-rpath,{hepmc3_libdir}",
    ]
    run_command(compile_command)

    generator_output = run_command(
        [
            str(executable_path),
            str(args.output),
            str(args.events),
            str(args.seed),
            f"{args.ecm:.17g}",
        ]
    )
    print(generator_output, end="")
    if not args.keep_source:
        source_path.unlink(missing_ok=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
