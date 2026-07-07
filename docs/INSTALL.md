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
