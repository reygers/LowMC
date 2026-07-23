# LowMC 0.1.0 Release Candidate

## Status

LowMC 0.1.0 is the first public research-prototype release candidate. It is
intended for controlled soft-photon studies and validation, not yet as a
production event-generator afterburner.

The generator-independent C++ core is the primary public interface. HepMC3 and
Pythia8 are optional event-interface integrations. PHOTOS is on a different
level: it is used only as an external comparison in the fixed
`Z -> mu+ mu-` benchmark and is not a LowMC event interface. The HepMC3
programs deliberately keep their `*_diagnostic` names in this release.

## Included

- coherent event-level leading-Low current and photon density;
- weighted fixed- and full-phase-space photon sampling;
- adaptive event-yield integration and Poisson multiplicity sampling;
- corrected event-specific shape unweighting with explicit envelope failure
  diagnostics;
- direct `ChargedParticle` C++ interface;
- a six-example learning path with descriptive filenames, split into
  dependency-free basic examples and advanced shape/HepMC3/Pythia8 examples;
- installed static core library, public headers, and CMake package target
  `LowMC::core`;
- optional HepMC3 selection, weighted diagnostics, and provisional physical or
  enhanced unweighted-photon output;
- optional direct Pythia8 selection and weighted diagnostics, including a
  direct-versus-HepMC3 cross-check;
- opt-in external PHOTOS comparison for the fixed `Z -> mu+ mu-` benchmark;
- DELPHI hadronic inner-bremsstrahlung benchmark;
- reproducible release-validation driver with JSON, Markdown, and per-command
  log output.

## Installation

The dependency-free core requires CMake 3.16 or newer and a C++17 compiler.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
cmake --install build --prefix /path/to/lowmc
```

A downstream CMake project can then use:

```cmake
find_package(LowMC 0.1 CONFIG REQUIRED)
target_link_libraries(my_target PRIVATE LowMC::core)
```

HepMC3 and Pythia8 remain opt-in:

```sh
cmake -S . -B build-integrations \
  -DCMAKE_BUILD_TYPE=Release \
  -DLOWMC_ENABLE_HEPMC3=ON \
  -DLOWMC_ENABLE_PYTHIA8=ON \
  -DCMAKE_PREFIX_PATH=/path/to/hepmc3 \
  -DPYTHIA8_CONFIG_EXECUTABLE=/path/to/pythia8/bin/pythia8-config
cmake --build build-integrations --parallel
ctest --test-dir build-integrations --output-on-failure
```

The installed HepMC3 afterburner example is documented in
[USAGE.md](USAGE.md#4-hepmc3-afterburner-output). Configuration of the
PHOTOS comparison benchmark and the full fresh-checkout guide are in
[INSTALL.md](INSTALL.md).

The complete example index and suggested learning order are in
[examples/README.md](../examples/README.md).

## Release Validation

The release candidate is validated with:

```sh
python3 scripts/validate_release_candidate.py \
  --enable-integrations \
  --hepmc3-prefix /path/to/hepmc3 \
  --pythia-prefix /path/to/pythia8 \
  --enable-photos \
  --photos-root /path/to/photospp
```

The retained evidence directory contains `manifest.json`, `summary.md`, one log
per command, and the `R2`/`R3` CSV and Markdown evidence. A valid release run
must report `dirty_tree=false`, pass the installed-core consumer, and pass all
applicable gates. The exact gate policy is in
[RELEASE_VALIDATION_MATRIX.md](RELEASE_VALIDATION_MATRIX.md).

The expected current test inventory is:

- `98` dependency-free tests;
- the `18`-row count-plus-shape gate;
- the strict shape-envelope gate;
- `192` tests with HepMC3 and Pythia8 enabled;
- `101` tests in the PHOTOS benchmark build.

Counts are informative; success requires every discovered test to pass,
including tests added after these notes were written.

## Scientific Validation

The frozen DELPHI comparison uses the complete charge-conserving event current
and targets DELPHI's calculated inner-bremsstrahlung reference, not the photon
excess. In the combined `3000`-event, `7523`-jet result,
`LowMC/direct = 1.00347`, `direct/DELPHI = 1.03218`, and the diagonal
`chi2/ndf = 12.08/10`. The primary `pT_gamma` comparison uses
`pT_gamma < 80 MeV`; the angular comparison follows the documented Figure 4
projection conventions. Details are in
[DELPHI_HADRONIC_BENCHMARK_PLAN.md](DELPHI_HADRONIC_BENCHMARK_PLAN.md).

For the fixed `Z -> mu+ mu-` benchmark, all six PHOTOS soft-window rate ratios
are statistically compatible with unity. LowMC/PHOTOS is
`0.995 +/- 0.019`, `1.010 +/- 0.013`, `1.013 +/- 0.010`,
`1.007 +/- 0.008`, `1.009 +/- 0.007`, and `1.009 +/- 0.007` as the upper
photon-energy cut increases from `0.06` to `0.50 GeV`. Details are in
[PHOTOS_BENCHMARK_PLAN.md](PHOTOS_BENCHMARK_PLAN.md).

These production-statistics studies are physics-change gates, not tests rerun
for every code-only release.

## Known Limitations

- LowMC 0.1.0 implements the leading real soft-photon approximation. It does
  not provide virtual corrections, exponentiation, or dielectron production.
- HepMC3 output adds stable photons without recoil and preserves the original
  event record. This is an explicit convention, not exact event-level
  four-momentum conservation after photon insertion.
- The default failed-event policy skips events with invalid selections,
  non-conserved selected charge, failed yields or envelopes, or insufficient
  candidates. Strict and atomic diagnostic modes are available.
- An event-specific sampled envelope is validated and failures are reported,
  but it is not a mathematical production guarantee for every possible event.
- Extremely boosted and high-multiplicity topologies can be inefficient. The
  `p = 1000 GeV` yield case remains a separately labelled stress case.
- The PHOTOS interface is a validation benchmark, not a general decay
  afterburner replacement. A VinciaQED comparison remains future independent
  validation.
- There is no stable production alias for the HepMC3 afterburner in this
  release.

## License

LowMC is distributed under the MIT License. See [LICENSE](../LICENSE).
