# LowMC

LowMC is an early-development C++ Monte Carlo library for coherent soft-photon
and soft-dielectron emission based on the leading Low theorem.

The project is designed as a generator-independent core with optional event
generator adapters. The core works on selected charged particles and returns
Low-photon samples, weights, and diagnostics. Optional HepMC3 and Pythia8 layers
translate event records into the core input format.

## Current Status

LowMC is a research prototype, not yet a production release.

Currently implemented:

- generator-independent soft-photon core
- coherent event-level Low current
- weighted photon sampling
- weighted soft-dielectron sampling through virtual photons and Kroll--Wada conversion
- fixed- and full-phase-space diagnostics
- HepMC3 adapter and provisional afterburner diagnostic
- Pythia8 adapter diagnostics
- direct charged-particle example

Still under validation:

- production unweighted-photon output
- final HepMC3 command-line interface
- broad physics and performance validation
- recoil and finite-mass corrections beyond the leading virtual Low current

## Repository Layout

```text
include/       public C++ headers
src/           generator-independent core implementation
tests/         core tests
diagnostics/   core diagnostic executables
adapters/      optional HepMC3 and Pythia8 integrations
fixtures/      small dependency-free input fixtures
examples/      user-facing examples
docs/          design, physics, usage, and development notes
scripts/       plotting, validation, and helper scripts
```

## Build And Test

The default build has no HepMC3 or Pythia8 dependency.

```sh
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
cmake --install build --prefix /path/to/lowmc
```

Installed CMake consumers use `find_package(LowMC 0.1 CONFIG REQUIRED)` and
link `LowMC::core`.

Build the optional HepMC3 adapter:

```sh
cmake -S . -B build-hepmc3 \
  -DLOWMC_ENABLE_HEPMC3=ON \
  -DCMAKE_PREFIX_PATH=/path/to/hepmc3
cmake --build build-hepmc3
```

Build the optional Pythia8 adapter:

```sh
cmake -S . -B build-pythia8 -DLOWMC_ENABLE_PYTHIA8=ON
cmake --build build-pythia8
```

## Examples

The repository contains a seven-example learning path under
[examples/](examples/README.md). Descriptive filenames are grouped into
`basic/` core examples and `advanced/` shape-unweighting, HepMC3, and Pythia8
examples. The basic executables are built by default; optional examples are
enabled with their corresponding adapters.

## Basic Usage

The most stable interface is the direct core workflow:

```text
selected charged particles
-> LowMC core
-> weighted Low-photon samples
-> user histograms or diagnostics
```

For stored event records, the optional HepMC3 layer can select charged
particles and run the same core. A provisional HepMC3 afterburner diagnostic can
also write output files containing the original event record plus generated Low
photons.

See [docs/USAGE.md](docs/USAGE.md) and
[docs/API_EXAMPLE.md](docs/API_EXAMPLE.md) for concrete examples.

## Physics Model

LowMC implements the leading soft-photon factor

```text
J^mu(k) = sum_i eta_i Q_i p_i^mu / (p_i . k)
S(k)    = -J_mu(k) J^mu(k)
```

where `Q_i` is the charge in units of the positron charge and `eta_i` encodes
the incoming/outgoing sign convention. The event-generator ensemble supplies
the hadronic event distribution; LowMC evaluates the event-level coherent
soft-photon factor.

See [docs/PHYSICS.md](docs/PHYSICS.md) for the full convention.

## Documentation

- [docs/DESIGN.md](docs/DESIGN.md): project design and boundaries
- [docs/PHYSICS.md](docs/PHYSICS.md): physics formulas and conventions
- [docs/INSTALL.md](docs/INSTALL.md): fresh-clone build and smoke-test guide
- [docs/USAGE.md](docs/USAGE.md): intended user workflows
- [docs/ALICE_DIELECTRON_BENCHMARK.md](docs/ALICE_DIELECTRON_BENCHMARK.md):
  ALICE Figure 4 data, Balek baseline, and LowMC comparison workflow
- [docs/ALGORITHM.md](docs/ALGORITHM.md): sampling and algorithm notes
- [docs/DECISIONS.md](docs/DECISIONS.md): recorded design decisions
- [docs/STATUS.md](docs/STATUS.md): current implementation status
- [docs/RELEASE_VALIDATION_MATRIX.md](docs/RELEASE_VALIDATION_MATRIX.md):
  release-blocking and scientific validation gates
- [docs/REPRODUCIBILITY.md](docs/REPRODUCIBILITY.md): software environment,
  random-number portability, units, concurrency, and benchmark settings
- [docs/RELEASE_NOTES_0.1.0.md](docs/RELEASE_NOTES_0.1.0.md): first public
  research-prototype release notes
- [examples/README.md](examples/README.md): basic and advanced examples
- [docs/NEXT_TASK.md](docs/NEXT_TASK.md): next planned development step

## License

LowMC is released under the MIT License. See [LICENSE](LICENSE).
