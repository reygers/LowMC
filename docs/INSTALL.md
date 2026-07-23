# INSTALL

This page is a short smoke-test guide for a fresh LowMC checkout.

LowMC is currently a research prototype. The default build intentionally avoids
optional event-generator dependencies.

## Requirements

Default build:

- CMake 3.16 or newer
- C++17 compiler
- Python 3, for helper-script tests

Optional adapters:

- HepMC3, only when `LOWMC_ENABLE_HEPMC3=ON`
- Pythia8, only when `LOWMC_ENABLE_PYTHIA8=ON`
- PHOTOS with its HepMC3 interface, only when `LOWMC_ENABLE_PHOTOS=ON`

## Default Build

From the repository root:

```sh
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
```

This should build the generator-independent core, core diagnostics, examples,
and default tests.

## HepMC3 Adapter

If HepMC3 is installed:

```sh
cmake -S . -B build-hepmc3 \
  -DLOWMC_ENABLE_HEPMC3=ON \
  -DCMAKE_PREFIX_PATH=/path/to/hepmc3
cmake --build build-hepmc3
ctest --test-dir build-hepmc3 -R hepmc3 --output-on-failure
```

For a local install under `$HOME/hepmc3`, use:

```sh
-DCMAKE_PREFIX_PATH=$HOME/hepmc3
```

## Pythia8 Adapter

If Pythia8 is installed and `pythia8-config` is available:

```sh
cmake -S . -B build-pythia8 -DLOWMC_ENABLE_PYTHIA8=ON
cmake --build build-pythia8
ctest --test-dir build-pythia8 -R pythia8 --output-on-failure
```

LowMC also searches the common local prefix `$HOME/pythia8317/bin` for
`pythia8-config`.

The opt-in Pythia residual-tail diagnostic compares the analytic leg-centred
proposal with the tuned grid on generated events:

```sh
./build-pythia8/adapters/pythia8/lowmc_pythia8_leg_proposal_tail_diagnostic \
  build-pythia8/leg_tail_ee_seed24681.csv \
  --mode ee-hadrons \
  --events 20 \
  --samples 20000 \
  --seed 24681
```

Use `--mode softqcd` for 13 TeV non-diffractive `pp` events. After running
multiple seeds, create the paired summary and diagnostic plot with:

```sh
python3 scripts/summarize_pythia_leg_tail_study.py \
  build-pythia8/leg_tail_*.csv \
  --summary build-pythia8/leg_tail_summary.csv \
  --paired build-pythia8/leg_tail_paired.csv \
  --plot build-pythia8/leg_tail_study.png
```

This remains a weighted diagnostic and does not construct or select a
production envelope.

The separate strict envelope diagnostic calibrates an empirical per-event
bound and checks it on independent proposal draws:

```sh
./build-pythia8/adapters/pythia8/lowmc_pythia8_leg_envelope_diagnostic \
  build-pythia8/leg_envelope_ee_seed35791.csv \
  --mode ee-hadrons \
  --events 20 \
  --calibration-samples 20000 \
  --trial-samples 100000 \
  --seed 35791
```

The default scan writes rows for safety factors `1.2`, `1.5`, and `2.0`.
Repeat `--safety-factor` to select a custom list. Summarize multiple seeds and
both event modes with:

```sh
python3 scripts/summarize_pythia_leg_envelope.py \
  build-pythia8/leg_envelope_*.csv \
  --summary build-pythia8/leg_envelope_summary.csv \
  --plot build-pythia8/leg_envelope_study.png
```

Any independent trial overshoot sets `envelope_ok=0`; it must be treated as a
failure. This command does not generate photons or replace the production
shape-unweighting envelope.

### macOS Pythia8 Troubleshooting

If the Pythia8 tests segfault very early, including in simple operations such
as `Pythia8::Event::reset()` or `append()`, suspect an ABI or toolchain
mismatch before assuming a LowMC physics bug.

On macOS, the most reproducible validation path so far is:

- build Pythia8 from source, for example Pythia 8.317 under `$HOME/pythia8317`
- use the same compiler family for Pythia8 and LowMC
- prefer `/usr/bin/clang++` for source builds on macOS, because it is the
  SDK-aware Apple compiler driver
- avoid mixing a package-manager Pythia8 build with a LowMC build that uses a
  different compiler, standard library, or SDK setup

MacPorts or Homebrew installations may work, but if adapter tests fail with
early crashes, first repeat the check with a source-built Pythia8 using the
same toolchain as LowMC.

## Combined Pythia8 And HepMC3 Check

If both optional dependencies are available:

```sh
cmake -S . -B build-pythia8-hepmc3 \
  -DLOWMC_ENABLE_PYTHIA8=ON \
  -DLOWMC_ENABLE_HEPMC3=ON \
  -DCMAKE_PREFIX_PATH=/path/to/hepmc3
cmake --build build-pythia8-hepmc3
ctest --test-dir build-pythia8-hepmc3 -R "pythia8|hepmc3" --output-on-failure
```

## Analytic Leg-Centred Proposal Diagnostic

The PHOTOS-inspired leg-centred proposal is part of the default build and does
not require PHOTOS. Write a fixed-`kT` weighted sample with:

```sh
./build/lowmc_fixed_kt_weighted_samples \
  build/leg_samples.csv \
  --proposal leg \
  --samples 100000 \
  --leg-uniform-fraction 0.1
```

Run the predeclared three-seed comparison against the tuned grid proposal on
the repository particle fixtures with:

```sh
python3 scripts/compare_leg_centered_proposal.py \
  build/leg_centered_proposal_comparison.csv \
  --benchmark build/lowmc_fixed_kt_weighted_samples \
  --samples 100000
```

This is a weighted diagnostic. It does not select the leg-centred proposal as
a production default or construct a rejection envelope.

## DELPHI Hadronic Bremsstrahlung Benchmark

With the optional Pythia8 adapter configured, build and run the paired direct
integral/LowMC diagnostic:

```sh
cmake --build build-pythia8 \
  --target lowmc_pythia8_delphi_hadronic_benchmark
./build-pythia8/adapters/pythia8/lowmc_pythia8_delphi_hadronic_benchmark \
  build-pythia8/delphi_hadronic_seed24680.csv \
  --events 20 \
  --seed 24680 \
  --mc-samples 300000 \
  --direct-loge-bins 48 \
  --direct-eta-bins 32 \
  --direct-phi-bins 48 \
  --direct-refinement-levels 3 \
  --direct-subbins 2
```

Aggregate event-level uncertainties and make the unscaled comparison to
DELPHI Table 4 and Figure 4b with:

```sh
MPLCONFIGDIR=build-pythia8/matplotlib \
python3 scripts/summarize_delphi_hadronic_benchmark.py \
  build-pythia8/delphi_hadronic_seed24680.csv \
  --summary build-pythia8/delphi_hadronic_summary.csv \
  --plot build-pythia8/delphi_hadronic.png \
  --final-plot build-pythia8/delphi_hadronic_final.png \
  --joint-inputs build-pythia8/delphi_hadronic_seed24680_joint.csv \
  --theta-summary build-pythia8/delphi_hadronic_theta_summary.csv \
  --theta-plot build-pythia8/delphi_hadronic_theta.png
```

Multiple independent input CSVs may be passed to the summarizer. It treats
Pythia events as the statistical units, uses a delete-one jackknife for the
per-accepted-jet ratio, and reports the LowMC photon-integration error
separately. The benchmark executable automatically creates the sibling
`_joint.csv`; use `--joint-output` to select another path. The Figure 4b
reference is vector-digitized and has no published covariance, so the angular
comparison is not assigned a chi-square. Figure 4b is the angular projection
over the displayed `0 < pT < 0.4 GeV` range. The theta summary also records a
separately labelled `pT < 80 MeV` projection and verifies that its integral
matches the primary Table 4-phase-space result. The diagnostic does not emit
photons or modify production defaults.

The default jet algorithm is LUCLUS (`--jet-algorithm Lund`) with
`djoin = 3 GeV`. The bounded DELPHI jet variations use
`--jet-algorithm Durham` or `--jet-algorithm Jade`; both select
`ycut = 0.01`. The chosen algorithm and scale are written to both CSV files.

For a fast same-event comparison of jet counts and axes, without evaluating
the Low current, run:

```sh
./build-pythia8/adapters/pythia8/lowmc_pythia8_delphi_hadronic_benchmark \
  build-pythia8/delphi_jet_accounting_seed24680.csv \
  --events 200 \
  --seed 24680 \
  --jet-accounting-only

python3 scripts/summarize_delphi_jet_accounting.py \
  build-pythia8/delphi_jet_accounting_seed24680.csv \
  --summary build-pythia8/delphi_jet_accounting_summary.csv \
  --plot build-pythia8/delphi_jet_accounting.png
```

The accounting mode selects events with the nominal LUCLUS-like Lund
configuration and records Lund, Durham, and JADE results for each identical
generated event.

## PHOTOS Benchmark

PHOTOS remains an opt-in benchmark dependency and is not linked into the
generator-independent core. With PHOTOS and HepMC3 installed, configure:

```sh
cmake -S . -B build-photos \
  -DLOWMC_ENABLE_PHOTOS=ON \
  -DPHOTOS_ROOT=/path/to/photospp \
  -DCMAKE_PREFIX_PATH=/path/to/hepmc3
cmake --build build-photos --target lowmc_photos_zmumu_benchmark
ctest --test-dir build-photos -R photos_zmumu --output-on-failure
```

Run the default fixed `Z -> mu+ mu-` pilot with:

```sh
./build-photos/lowmc_photos_zmumu_benchmark \
  build-photos/photos_zmumu_photons.csv \
  build-photos/photos_zmumu_summary.csv \
  --lowmc-output build-photos/lowmc_zmumu_photons.csv
```

Create the normalized shape and ratio comparison with:

```sh
python3 scripts/plot_photos_zmumu_comparison.py \
  build-photos/photos_zmumu_photons.csv \
  build-photos/lowmc_zmumu_photons.csv \
  build-photos/photos_zmumu_shapes.png \
  --bins-csv build-photos/photos_zmumu_shape_bins.csv \
  --summary-csv build-photos/photos_zmumu_shape_summary.csv
```

Run the predeclared three-seed soft-energy scan with:

```sh
python3 scripts/run_photos_zmumu_soft_scan.py \
  build-photos/photos_zmumu_soft_scan.csv \
  build-photos/photos_zmumu_soft_scan.png \
  --benchmark build-photos/lowmc_photos_zmumu_benchmark \
  --events 100000
```

The default scan uses six nested values
`E_gamma_max = 0.06, 0.075, 0.10, 0.15, 0.25, 0.50 GeV` and three independent
PHOTOS seed pairs. The resulting points are correlated because the same seed
pairs and nested selections are intentionally used at every energy maximum.

The PHOTOS 3.64 fixed-order mode used here is not stable at arbitrarily small
generation cutoffs. The benchmark default uses a fractional cutoff of
`3e-4`, corresponding to about `27.4 MeV` for the chosen `Z` mass, and applies
the common analysis only above `kT = 50 MeV`.

## What To Report

For installation feedback, please include:

- operating system and version
- compiler and version
- CMake version
- whether the default build passed
- whether optional HepMC3 or Pythia8 builds were attempted
- the first failing command and its full error output

## Private Material

Large literature files, slides, and local notes are not part of the public
repository. The repository should be buildable from the tracked files alone.
