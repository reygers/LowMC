# USAGE

This document describes how users are expected to use LowMC.

LowMC is primarily a generator-independent C++ library. The core takes selected
charged particles, a photon phase space, and sampler configuration, then returns
Low-photon samples and diagnostics. Event-record adapters and command-line
programs are layers around that core.

## Main Workflows

### 1. Direct Core API

Use this path for toy events, tests, and custom analyses that already know which
charged particles should enter the coherent Low current.

```text
user code
-> construct selected ChargedParticle objects
-> configure PhotonPhaseSpace
-> call LowMC core
-> fill histograms or inspect diagnostics
```

The direct API is the simplest and most stable interface for controlled
studies. See `docs/API_EXAMPLE.md` for a concrete C++ example.

### 2. HepMC3 Input And Histograms

Use this path for stored event samples or generator output written as HepMC3.

```text
HepMC3 file
-> optional LowMC HepMC3 adapter selects charged particles
-> LowMC core generates weighted photons or diagnostics
-> user analysis fills ROOT, Python, or custom histograms
```

In this mode LowMC does not need to write a new HepMC3 file. The current
HepMC3 diagnostic reads events, reports adapter diagnostics, and can call the
weighted Low-photon path. Weighted photons are analysis samples: histograms
must be filled with their weights.

### 3. Pythia Workflows

Pythia can be supported in two ways:

```text
Pythia
-> HepMC3 output
-> LowMC HepMC3 workflow
```

or, later:

```text
Pythia event in memory
-> optional LowMC Pythia adapter
-> LowMC core
-> user analysis
```

The direct Pythia adapter should remain a thin optional layer. It must not
introduce a Pythia dependency into `lowmc_core`. Its boundary is planned as:

```text
Pythia8::Event
-> AdapterParticleCandidate records
-> dependency-free LowMC selection helper
-> selected ChargedParticle list
-> weighted LowMC photons or diagnostics
```

The first direct Pythia diagnostic should focus on selected-particle
diagnostics. Weighted photon generation can follow after the Pythia
incoming/final-state mapping is validated. HepMC3 output writing is not part of
the direct Pythia workflow.

Build the optional direct Pythia adapter with:

```sh
cmake -S . -B build-pythia8 -DLOWMC_ENABLE_PYTHIA8=ON
cmake --build build-pythia8 --target lowmc_pythia8_adapter_diagnostic
```

The optional build uses `pythia8-config`; keep Pythia out of the default build
unless this adapter is explicitly requested.

Run the deterministic toy selected-particle diagnostic:

```sh
./build-pythia8/adapters/pythia8/lowmc_pythia8_adapter_diagnostic --mode toy
```

Run one generated SoftQCD event in memory:

```sh
./build-pythia8/adapters/pythia8/lowmc_pythia8_adapter_diagnostic \
  --mode softqcd --events 1 --seed 12345
```

To call the weighted Low-photon path directly on the selected Pythia event,
build and run:

```sh
cmake --build build-pythia8 --target lowmc_pythia8_weighted_photon_diagnostic
./build-pythia8/adapters/pythia8/lowmc_pythia8_weighted_photon_diagnostic \
  --mode softqcd --events 1 --samples 4 --seed 12345
```

This prints weighted photon summaries only. It does not write HepMC3 output.

To cross-check the direct Pythia selected-particle path against the
Pythia-to-HepMC3 path for the same generated event, enable both optional
adapters:

```sh
cmake -S . -B build-pythia8-hepmc3 \
  -DLOWMC_ENABLE_PYTHIA8=ON \
  -DLOWMC_ENABLE_HEPMC3=ON \
  -DCMAKE_PREFIX_PATH=$HOME/hepmc3
cmake --build build-pythia8-hepmc3 --target lowmc_pythia8_hepmc3_crosscheck
./build-pythia8-hepmc3/adapters/pythia8/lowmc_pythia8_hepmc3_crosscheck \
  --events 3 --samples 4 --seed 12345
```

This diagnostic also compares weighted photon summaries for the same generated
events using the same sampling seed and phase-space configuration. It compares
selected-particle diagnostics and weighted summary fields, not event-record
particle ordering, because the Pythia-to-HepMC3 conversion may reorder records.

For manual validation with a local Pythia installation, LowMC provides an
optional helper:

```sh
python3 scripts/generate_pythia_hepmc3_sample.py \
  --output build-hepmc3/pythia_softqcd_10.hepmc3 \
  --events 10 \
  --seed 12345
```

The helper compiles a tiny temporary Pythia8 + HepMC3 generator in its work
directory, using local `pythia8-config` and `HepMC3-config` paths. It is not
part of the default build or test suite and does not make Pythia a LowMC
dependency.

### 4. HepMC3 Afterburner Output

Use this path when the desired output is a new HepMC3 file containing the
original event record plus added unweighted LowMC photons.

```text
input.hepmc3
-> LowMC selects charged particles
-> LowMC samples physical unweighted Low photons
-> output.hepmc3 with added stable photons
```

This path is currently exposed through the optional diagnostic executable
`lowmc_hepmc3_afterburner_diagnostic`. The executable name is intentionally
conservative: the file format and basic gates are in place, but broader
production validation is still ongoing.

For now this is the only HepMC3 afterburner command name. A shorter production
alias, such as `lowmc_hepmc3_afterburner`, should wait until validation gates
and output-mode conventions are more mature.

Build with the optional HepMC3 adapter enabled:

```sh
cmake -S . -B build-hepmc3 -DLOWMC_ENABLE_HEPMC3=ON -DCMAKE_PREFIX_PATH=/path/to/hepmc3
cmake --build build-hepmc3 --target lowmc_hepmc3_afterburner_diagnostic
```

To install the optional HepMC3 command-line tools:

```sh
cmake --install build-hepmc3 --prefix /path/to/lowmc-install
```

This currently installs the provisional diagnostic names into `bin/`:

- `lowmc_hepmc3_adapter_diagnostic`
- `lowmc_hepmc3_selected_particle_text`
- `lowmc_hepmc3_weighted_photon_diagnostic`
- `lowmc_hepmc3_afterburner_diagnostic`
- `lowmc_hepmc3_output_summary_diagnostic`

Stable shorter aliases should wait until the HepMC3 output contract and
validation gates are less provisional.

For a minimal end-to-end check using the checked-in charge-conserving fixture,
run:

```sh
python3 scripts/hepmc3_quickstart.py \
  --tool-dir build-hepmc3/adapters/hepmc3 \
  --output-dir build/hepmc3_quickstart
```

The script reads the fixture, writes a new HepMC3 output file under the chosen
output directory, and then runs the output-summary diagnostic on that file. For
installed tools, pass the installed `bin/` directory as `--tool-dir`. The
quickstart uses the afterburner's `--atomic-output` mode by default, because it
is a user-facing repository example rather than a low-level diagnostic
streaming test. To exercise the diagnostic streaming mode explicitly, pass
`--no-atomic-output`.

The repository quickstart remains useful for build-tree checks. Installed
packages also provide a self-contained HepMC3 example directory:

```text
<prefix>/share/lowmc/examples/hepmc3/
  README.md
  hepmc3_quickstart.py
  charge_conserving_event.hepmc3
```

After installation, run the installed example with:

```sh
cd /path/to/lowmc-install/share/lowmc/examples/hepmc3
python3 hepmc3_quickstart.py --output-dir quickstart-output
```

The installed quickstart uses the installed `bin/` tools from the same prefix
by default and writes output only to the selected output directory.

Run on a HepMC3 ASCII input file:

```sh
./build-hepmc3/adapters/hepmc3/lowmc_hepmc3_afterburner_diagnostic \
  input.hepmc3 output-with-lowmc.hepmc3
```

The output file preserves the original HepMC3 event record and appends accepted
LowMC photons. The first HepMC3 topology is:

- stable photons with PDG id `22`
- status code `1`
- no recoil applied to original particles
- original event vertices and event weights unchanged
- `LowMC:` provenance attributes stored directly on each added photon
- no new afterburner production vertex in the first durable convention

The first topology uses attributed stable root photons because this
representation round-trips cleanly through HepMC3 ASCII and avoids implying a
physical production vertex or recoil model. A dedicated afterburner-vertex mode
may be revisited later if a portable HepMC3 convention is found and tested.

The default failed-event policy is `skip`: events with invalid selection,
non-conserved signed charge, failed yield gates, failed shape envelopes, or
insufficient accepted shape candidates are written unchanged and reported in
the text diagnostics. For strict runs, use:

```sh
./build-hepmc3/adapters/hepmc3/lowmc_hepmc3_afterburner_diagnostic \
  --failure-policy fail input.hepmc3 output-with-lowmc.hepmc3
```

The current diagnostic strict mode is not an atomic output-file contract. If a
later event fails after earlier events were written, the process exits nonzero
and the output path may contain only the already-written events. Treat such a
file as diagnostic partial output, not as a completed LowMC event sample.

For an explicit atomic diagnostic run, use:

```sh
./build-hepmc3/adapters/hepmc3/lowmc_hepmc3_afterburner_diagnostic \
  --atomic-output input.hepmc3 output-with-lowmc.hepmc3
```

With `--atomic-output`, LowMC writes to a deterministic sibling temporary path
`<output>.tmp` and publishes the requested output path only after the full input
file succeeds. If the run fails, no new requested output is published and the
temporary path is removed by the diagnostic.

The default physical HepMC3 output uses enhancement factor `1`:
`N_gamma ~ Poisson(lambda_event)`. For larger Low-photon statistics in HepMC3
samples, use the explicit enhanced mode:

```sh
./build-hepmc3/adapters/hepmc3/lowmc_hepmc3_afterburner_diagnostic \
  --enhancement-factor 10 input.hepmc3 output-enhanced-lowmc.hepmc3
```

Enhanced mode draws
`N_gamma ~ Poisson(enhancementFactor * lambda_event)`. The first implementation
stores the enhancement metadata directly on each added Low photon and duplicates
the event-level mode metadata on events that receive LowMC photons:
`LowMC:photon_mode = "enhanced"`,
`LowMC:enhancement_factor = F`, and `LowMC:weight = 1 / F`.
Skipped events are written unchanged and do not receive LowMC event metadata.

In enhanced mode, direct photon counting gives the enhanced sample yield.
Filling histograms with `LowMC:weight` recovers the physical Low-photon yield.
Weighted diagnostic samples remain a separate concept and should be marked as
diagnostic rather than ordinary enhanced event-record photons.

To inspect an output file without modifying it, use the summary diagnostic:

```sh
./build-hepmc3/adapters/hepmc3/lowmc_hepmc3_output_summary_diagnostic \
  output-enhanced-lowmc.hepmc3
```

It prints one compact line per event with the number of LowMC photons, photon
mode, enhancement factor, summed photon weight, and whether event-level LowMC
metadata is present. Events skipped by the afterburner should appear with zero
LowMC photons and no LowMC event metadata.

This workflow depends on the current validated pieces:

- shape-unweighted photon sampling
- event-level multiplicity and normalization

The accepted first multiplicity convention is
`N_gamma ~ Poisson(lambda_event)`, where `lambda_event` is the integrated Low
yield in the configured `PhotonPhaseSpace`.

### HepMC3 Prototype Readiness

The HepMC3 afterburner is ready to describe as a first controlled usable
prototype: it is suitable for small HepMC3 files, installed-example runs,
manual validation studies, and user experiments where the diagnostic status of
the command is understood.

Safe for users to try:

- build-tree and installed diagnostic executables
- repository and installed quickstarts
- physical output with enhancement factor `1`
- explicitly enhanced output with `--enhancement-factor F` and documented
  `LowMC:weight = 1 / F`
- explicit atomic output with `--atomic-output`
- Pythia-produced HepMC3 files for small, bounded validation runs

Still diagnostic or provisional:

- the command name remains `lowmc_hepmc3_afterburner_diagnostic`
- the stable-root-photon topology is the first durable convention, not a final
  claim about the best long-term HepMC3 topology
- larger Pythia-produced multi-event afterburner runs can be slow and should
  remain bounded/manual; use `scripts/validate_pythia_afterburner.py` for a
  repeatable opt-in validation summary
- naturally occurring shape-envelope failures in realistic files have not yet
  been characterized
- direct Pythia in-memory selected-particle and weighted-photon diagnostics
  exist, and the bounded three-event direct-vs-HepMC3 cross-check matches for
  both selected-particle diagnostics and weighted photon summaries with seed
  `12345`

Therefore this is not yet a production afterburner. A production-style command
or alias should wait for broader runtime validation, natural failure-mode
evidence, and final command naming.

### HepMC3 Prototype Validation Status

The current HepMC3 afterburner path is validated as a first controlled
prototype for small HepMC3 ASCII fixtures.

Implemented checks cover:

- single-event physical output with one attributed LowMC photon
- enhanced output with `--enhancement-factor 10` and per-photon
  `LowMC:weight = 0.1`
- multi-event output where a charge-conserving event receives LowMC photons
  and a non-charge-conserving event is written unchanged
- a three-event, inspectable HepMC3 fixture with two charge-conserving events
  that receive enhanced LowMC photons and one skipped non-conserved event
- a generated ten-event HepMC3 sample with varied pp-like final states: nine
  charge-conserving events receive enhanced LowMC photons and one deliberate
  non-conserved event is skipped
- strict `--failure-policy fail` coverage for the same three-event fixture,
  including the current partial diagnostic output behavior
- explicit `--atomic-output` coverage for successful output and strict
  failure without publishing the requested output file
- invalid or empty charged-particle selection using a neutral-only event:
  default `skip` output is unchanged, and `--failure-policy fail` returns
  nonzero
- insufficient accepted shape candidates using a deliberately low-trial
  charge-conserving event: default `skip` output is unchanged, and
  `--failure-policy fail` returns nonzero
- shape-envelope failure using deliberately unsafe diagnostic safety factors:
  default `skip` output is unchanged, and `--failure-policy fail` returns
  nonzero
- a small deterministic count/yield smoke scan over physical and enhanced
  output for the single-event and realistic three-event fixtures, checking
  requested/added photon counts, summary counts, skipped events, and reported
  `lambda_event` values at loose Poisson precision
- a tighter generated-sample count/yield scan over two generated events, two
  seeds, and enhancement factors `10` and `25`, checking the same
  requested/added/summary consistency with a two-sigma aggregate Poisson band
- an opt-in high-statistics count/yield validation script,
  `scripts/hepmc3_count_yield_highstat.py`, which generates a reproducible
  HepMC3 sample by default and prints an aggregate validation line suitable for
  manual logs
- a manual Pythia8-to-HepMC3 validation path: SoftQCD samples can be generated
  with the local Pythia installation, adapter diagnostics read the generated
  events, and a bounded five-event afterburner validation with enhancement
  factor `10` has passed; all five events produced enhanced LowMC photons, all
  shape diagnostics were `ok`, no events were skipped, and the read-only
  output-summary diagnostic matched the afterburner text diagnostics
- `scripts/validate_pythia_afterburner.py`, an opt-in helper that wraps sample
  generation, the HepMC3 afterburner, and the output-summary diagnostic, then
  writes logs plus `validation_summary.csv`
- a first three-seed helper scan with two Pythia-produced events per seed:
  all six events passed with no skips, all shape statuses `ok`, matching
  output-summary diagnostics, and visible afterburner runtime spread from
  about `15 s` to `43 s` per two-event sample
- a boosted transverse `p = 100 GeV` HepMC3 fixture that produces enhanced
  LowMC photons with valid metadata and successful shape diagnostics
- a forward pion fixture with `p_T = 1 GeV` and `|eta| = 5` that produces
  enhanced LowMC photons with valid metadata and successful shape diagnostics
- an explicit transverse `p = 1000 GeV` pion edge fixture that produces
  enhanced LowMC photons with valid metadata and successful shape diagnostics
- preservation of original vertices, original event weights, and no-recoil
  behavior in the checked fixtures
- per-photon and event-level LowMC metadata round-tripping through HepMC3 ASCII
- the read-only output-summary diagnostic
- installation of the optional HepMC3 diagnostic executables
- the repository quickstart helper using the checked-in fixture
- the installed HepMC3 quickstart example

Remaining validation gaps before calling the afterburner a first usable
prototype:

- better runtime/complexity reporting for bounded Pythia-produced afterburner
  scans
- validation of naturally occurring shape-envelope failures in realistic event
  files, beyond the current forced diagnostic-safety test
- broader boosted-particle HepMC3 stress tests beyond the current single-event
  transverse and forward edge fixtures
- decision on whether the current diagnostic command names should remain or
  receive stable production aliases

The minimum checklist for promoting the diagnostic behavior to a future
production-style HepMC3 command or alias is maintained in `docs/ROADMAP.md`
under "HepMC3 Production-Afterburner Readiness".

## Weighted Versus Physical Photons

Weighted photons are for analysis and validation. Their `Photon::weight`
contributes to histograms. Averaging weighted histograms over input events
estimates inclusive Low-photon spectra.

Physical event-record photons are different. They should be ordinary stable
photons with no analysis weight, sampled from the normalized Low shape and from
the physical event multiplicity model.

Diagnostic oversampling may use
`N_sample ~ Poisson(oversamplingFactor * lambda_event)` with per-sample weight
`1 / oversamplingFactor`. This improves histogram statistics but does not
change the physical photon yield.

## Boosted Particles And Validity Diagnostics

LowMC should not impose a fixed upper momentum cutoff on charged particles. The
Low formula itself has no such numerical cutoff, but the angular structures near
charged-particle directions become narrower as the scale `m / E` decreases.
Highly boosted particles can therefore make unweighting inefficient or expose an
insufficient event-specific envelope.

The intended user-facing contract is:

- weighted generation remains the robust fallback for difficult events,
- unweighted event-record photons require a validated event-specific envelope,
- LowMC reports envelope failures and efficiency diagnostics rather than
  silently writing untrusted unweighted photons,
- papers and manuals should quote benchmarked and stress-tested configurations
  instead of advertising a hard allowed momentum range.

The transverse `p = 1000 GeV` two-pion toy is currently treated as an extreme
stress test. It is useful for finding boosted-tail failures, but it should not
alone define the default production strategy.

## Package Shape

The intended package structure is:

- `lowmc_core`: generator-independent physics and sampling library
- optional adapters, such as HepMC3 and later Pythia8
- small reference executables and diagnostics
- optional plotting scripts for CSV or text outputs

The core should remain independent of HepMC3, Pythia8, ROOT, and plotting
libraries.
