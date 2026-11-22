# RELEASE VALIDATION MATRIX

## Scope

This document defines the validation contract for the first public LowMC
research-prototype release. It does not promote the provisional HepMC3
diagnostic to a production afterburner, change physics cuts, or change sampling
defaults.

The release is component based:

- the dependency-free core gates every release;
- HepMC3 and Pythia8 gates are mandatory when those integrations are advertised
  or shipped;
- PHOTOS is an optional external benchmark dependency, not a LowMC runtime
  dependency;
- DELPHI, PHOTOS, analytic-bound, rapidity-center, proposal, and performance
  studies are retained as scientific evidence rather than being rerun for every
  code-only release.

The checked-in GitHub Actions workflow runs the dependency-free R1 core and
installed-consumer gate with GCC on Ubuntu for every push and pull request.
The optional integration and scientific gates below remain explicit release
checks because they require external packages or frozen benchmark inputs.

## Gate Classes

| Class | Meaning | Failure consequence |
|---|---|---|
| `R` | Release-blocking for every release | Do not release |
| `C` | Release-blocking when the named optional component is advertised or changed | Do not advertise that component as validated |
| `P` | Physics-change gate | Do not release a change to the affected physics path |
| `S` | Scientific benchmark or stress diagnostic | Record and understand changes; not automatically release-blocking |

## Clean Release Candidate

Use new build directories and a Release build. Replace the example dependency
prefixes with the actual installations.

### R1: Repository and dependency-free core

```sh
git diff --check
cmake -S . -B build-release-default -DCMAKE_BUILD_TYPE=Release
cmake --build build-release-default --parallel
ctest --test-dir build-release-default --output-on-failure
cmake --install build-release-default \
  --prefix build-release-default/install
cmake -S tests/install_consumer \
  -B build-release-default/install-consumer \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH=$PWD/build-release-default/install
cmake --build build-release-default/install-consumer --parallel
ctest --test-dir build-release-default/install-consumer --output-on-failure
```

Pass criteria:

- `git diff --check`, configure, build, test, and install all exit successfully;
- every discovered CTest passes;
- an external CMake consumer finds `LowMC 0.1`, links `LowMC::core` using only
  the install prefix, includes `<lowmc/LowMC.h>`, and runs successfully;
- no optional event-generator dependency is required by this configuration.

This covers the Low current and density, phase-space conventions, weighted
sampling, adaptive yield, Poisson count sampling, corrected shape unweighting,
the analytic audit implementation, diagnostics, and the direct example.
Script `--help` tests are interface smoke tests, not substitutes for the
physics-change gates below.

### R2: Count-plus-shape integration

```sh
python3 scripts/validate_count_shape.py \
  build-release-evidence/count-shape.csv \
  --aggregate-csv build-release-evidence/count-shape-aggregate.csv \
  --summary-md build-release-evidence/count-shape.md \
  --output-dir build-release-evidence/count-shape-runs \
  --diagnostic-exe build-release-default/lowmc_count_shape_diagnostic \
  --seeds 3 \
  --seed 12345 \
  --oversampling-factor 100
```

Pass criteria:

- all 18 rows report `status=ok`, `envelope_ok=1`, and
  `enough_candidates=1`;
- requested and retained photon counts agree;
- there are zero envelope failures.

This is the central end-to-end core gate for the unweighted-photon path.
The exact retained counts are seed-dependent evidence, not fixed acceptance
criteria.

### R3: Strict shape-envelope suite

```sh
python3 scripts/validate_shape_unweighting_candidates.py \
  build-release-evidence/shape-envelope.csv \
  --summary-md build-release-evidence/shape-envelope.md \
  --work-dir build-release-evidence/shape-envelope-runs \
  --truth-exe build-release-default/lowmc_two_particle_diagnostic \
  --shape-exe build-release-default/lowmc_shape_unweighting_diagnostic \
  --candidate-set targeted \
  --trials 30000 \
  --seed 12345 \
  --seed 23456 \
  --seed 34567 \
  --gate-policy strict \
  --gate-max-trial-bound-ratio 0.8
```

Pass criteria are encoded by strict mode: zero envelope failures, zero fallback
uses, sufficient accepted candidates, and an aggregate maximum
trial-to-bound ratio below `0.8`. Shape distances remain informative because
some fixed-grid truth histograms are resolution limited.

### C1: HepMC3 and Pythia8 integrations

The first public release currently documents both integrations, so this
conditional gate is required for that release.

```sh
cmake -S . -B build-release-integrations \
  -DCMAKE_BUILD_TYPE=Release \
  -DLOWMC_ENABLE_HEPMC3=ON \
  -DLOWMC_ENABLE_PYTHIA8=ON \
  -DCMAKE_PREFIX_PATH=/path/to/hepmc3
cmake --build build-release-integrations --parallel
ctest --test-dir build-release-integrations --output-on-failure
cmake --install build-release-integrations \
  --prefix build-release-integrations/install
```

Pass criteria:

- configure, build, all discovered tests, and installation succeed;
- HepMC3 selection, units, charge lookup, weighted sampling, physical and
  enhanced output, metadata round-trip, preservation, skip/fail behavior,
  atomic output, generated samples, count/yield checks, and build-tree and
  installed afterburner examples pass;
- Pythia toy and SoftQCD selection and weighted diagnostics pass;
- the direct-Pythia versus Pythia-to-HepMC3 cross-check passes.

Missing HepMC3 or Pythia8 is not a failure of the dependency-free core, but a
release made without this gate must not claim that the omitted integration was
validated for that release.

## Physics-Change Gates

These gates are required only when the listed physics implementation,
normalization, selection, or default configuration changes.

| ID | Trigger | Command or reference | Required outcome |
|---|---|---|---|
| `P1` | Low density, phase-space measure, weighted sampler, or yield estimator | `python3 scripts/validate_event_yield.py build-release-evidence/event-yield.csv --summary-md build-release-evidence/event-yield.md --output-dir build-release-evidence/event-yield-runs --yield-exe build-release-default/lowmc_event_yield_diagnostic --weighted-exe build-release-default/lowmc_full_weighted_sampler_study --samples 20000 --seed 12345` | Ordinary cases retain the documented adaptive-yield agreement; the `p=1000 GeV` case remains a separately labelled stress case |
| `P2` | Low current, signs, normalization, DELPHI selection, jet association, or per-jet normalization | Follow `docs/DELPHI_HADRONIC_BENCHMARK_PLAN.md` | Preserve or explain the recorded agreement with DELPHI inner bremsstrahlung; do not compare against the anomalous excess as a LowMC pass criterion |
| `P3` | Low current, decay normalization, finite-mass treatment, PHOTOS benchmark implementation, or common acceptance | Follow `docs/PHOTOS_BENCHMARK_PLAN.md` | Preserve or explain soft-window compatibility with unity; harder-window differences remain model diagnostics |

`P2` and `P3` use statistical criteria and recorded configurations, not
bit-for-bit output. A statistically significant change blocks the affected
physics change until it is classified and documented.

## Optional External PHOTOS Build

When the PHOTOS benchmark implementation is changed, also run:

```sh
cmake -S . -B build-release-photos \
  -DCMAKE_BUILD_TYPE=Release \
  -DLOWMC_ENABLE_PHOTOS=ON \
  -DPHOTOS_ROOT=/path/to/photospp \
  -DCMAKE_PREFIX_PATH=/path/to/hepmc3
cmake --build build-release-photos --parallel
ctest --test-dir build-release-photos --output-on-failure
```

All discovered tests must pass. This checks that the optional benchmark still
builds and its fixed-configuration smoke tests run; it does not replace the
production-statistics comparison in `P3`.

## Scientific Evidence, Not Routine Release Gates

The following remain valuable but do not block an unrelated release:

- the DELPHI full-statistics reproduction after its result is frozen;
- the PHOTOS production soft-window scan after its result is frozen;
- VinciaQED comparison, which is still future independent validation;
- Stefan's analytic global/cell bounds, which passed safety but are too loose
  to replace the current envelope;
- true-rapidity targeting, which did not justify a default change;
- PHOTOS-inspired leg-centered and possible Vincia-inspired pair proposals;
- reduced-probe reviews, including the diagnostic-only `5 x 5` study;
- extreme boosted cases, proposal comparisons, and performance scans.

These studies become `P` gates if a proposed change directly affects their
conclusion. Performance is currently recorded without a hardware-independent
hard threshold; a correctness failure always blocks, while a runtime change
requires investigation and documentation.

## Failure Handling

- Any deterministic `R` or applicable `C` failure blocks release.
- Fixed-seed stochastic gates are reproducible gates. Do not convert a failing
  run into a pass by rerunning until it succeeds.
- An unavailable optional dependency means “not validated,” not “failed core.”
- The documented limitations of the `p=1000 GeV` yield stress case are not a
  release failure unless a change worsens or misrepresents them.
- Expected negative tests pass only when they fail in the expected way and
  preserve the documented output contract.

## Current Evidence Snapshot

On 2026-07-23, the complete driver implementation check reported:

- the dependency-free Release build passed all `98` tests and installed;
- the `R2` count-plus-shape gate passed all `18` rows;
- the `R3` strict shape-envelope gate passed;
- the combined HepMC3 plus Pythia8 Release build passed all `192` tests;
- the PHOTOS 3.64 Release build passed all `101` tests.

These counts are evidence, not permanent pass criteria; newly added tests must
also pass. This implementation check used an uncommitted development tree; the
tagged release candidate requires a retained clean-tree manifest.

## Validation Driver

`scripts/validate_release_candidate.py` implements this matrix without changing
its physics criteria. By default it runs `R1`, `R2`, and `R3` in fresh
timestamped build directories. Explicit switches add the combined `C1` gate
and the opt-in external PHOTOS benchmark build. Every run preserves command logs,
`manifest.json`, and `summary.md`; it stops on the first blocking failure.

The driver deliberately does not run the DELPHI or PHOTOS
production-statistics studies. Their trigger conditions remain the `P2` and
`P3` physics-change gates above.
