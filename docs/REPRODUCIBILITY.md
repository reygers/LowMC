# Reproducibility and portability

This note records details that affect numerical reproduction of LowMC results.
It complements the release gates in
[`RELEASE_VALIDATION_MATRIX.md`](RELEASE_VALIDATION_MATRIX.md).

## Reference software environment

The paper benchmarks were produced with:

- LowMC configured as a C++17 Release build;
- Pythia 8.317;
- HepMC3 3.03.01;
- PHOTOS 3.64 for the independent fixed-Born comparison only;
- Apple clang 21.0.0 and CMake/CTest 3.31.3 on macOS.

The dependency-free core is also built and tested with GCC on Ubuntu by the
checked-in GitHub Actions workflow. Pythia8, HepMC3, and PHOTOS are optional;
the default core build does not need them.

## Random-number generation

Pseudorandom LowMC paths instantiate a local `std::mt19937_64` engine from the
seed stored in the relevant configuration. Uniform draws use
`std::generate_canonical<double, 53>`. Photon multiplicities use
`std::poisson_distribution<int>`.

A fixed seed reproduces a run made with the same build and C++ standard-library
implementation. The C++ standard does not require standard distributions, or
all helper mappings from engine output to floating point, to produce identical
sequences in different library implementations. Seeds are therefore important
provenance, but are not a promise of bitwise identity across libstdc++ and
libc++. Physics results should be compared with numerical tolerances or
statistical uncertainties.

## Units and conventions

- four-vector components, masses, transverse momenta, and energies: GeV;
- angles: radians;
- pseudorapidity and rapidity: dimensionless;
- charges: units of the positron charge `e`;
- cross sections: not part of the generator-independent core; benchmark
  adapters state their units explicitly, normally mb.

The component order is `(E, px, py, pz)`. Incoming and outgoing signs are
stored separately from electric charge.

## Concurrency and memory

The core has no hidden mutable event state and creates its random engines
locally. Independent calls can therefore run concurrently when their input and
result objects are not shared. LowMC does not currently provide an internal
thread pool.

Weighted generators retain the requested samples in their result objects, so
their memory use grows linearly with the stored sample count. Adaptive
integration and shape generation additionally store event-specific grid and
envelope data; their memory use grows with the configured number of cells.
Users processing large samples should consume and release one event result at
a time unless they intentionally need to retain it.

## Testing scope

The default CTest suite covers the current and density, phase-space
conventions, weighted photon and dielectron sampling, integration,
multiplicity, shape-envelope diagnostics, install/export behavior, and basic
script interfaces. Conditional suites add Pythia8 and HepMC3 adapter tests;
PHOTOS is a benchmark dependency rather than a runtime dependency. The release
matrix records the exact commands and physics-change gates. No line- or
branch-coverage percentage is currently claimed.

## Benchmark generator settings

The DELPHI hadronic-Z diagnostic records every explicit Pythia setting in its
source. It uses Pythia 8.317 defaults except for the following changes:

```text
PDF:lepton = off
WeakSingleBoson:ffbar2gmZ = on
23:onMode = off
23:onIfAny = 1 2 3 4 5
Beams:idA = 11
Beams:idB = -11
Beams:eCM = mZ from the Pythia particle data
TimeShower:QEDshowerByQ = off
TimeShower:QEDshowerByL = off
TimeShower:QEDshowerByOther = off
TimeShower:QEDshowerByGamma = off
SpaceShower:QEDshowerByQ = off
SpaceShower:QEDshowerByL = off
HadronLevel:QED = off
```

The generator-level approximation requires at least five final charged
particles, visible energy above `0.2 Ecm`, and a thrust axis between 30 and
150 degrees. LUCLUS-like truth jets use a 3 GeV joining scale; accepted jets
have momentum at least 5 GeV, lie between 30 and 150 degrees, and contain no
final electron or positron. The photon cuts and normalization are documented
in [`DELPHI_HADRONIC_BENCHMARK_PLAN.md`](DELPHI_HADRONIC_BENCHMARK_PLAN.md).
