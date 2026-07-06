# ADAPTER PLAN

This document records the first event-record adapter planning pass. It should
stay focused on generator-specific mapping questions; stable core boundaries
belong in `docs/ADAPTER_CONTRACT.md`.

## First Adapter Choice

The first real event-record adapter should target HepMC3.

Reasons:

- HepMC3 is an event-record format, not a generator-specific runtime interface.
- It keeps LowMC usable for stored samples and for events produced by different
  generators.
- It matches the roadmap order: HepMC3 adapter first, then a thin Pythia8
  interface around the same core.
- It is the cleanest test of the adapter contract because the core can remain
  independent of generator execution.

Pythia8 should remain the second integration target. Its adapter can either
translate directly from the Pythia event record or reuse a HepMC3 output path,
but that choice should not affect the core API.

## Required Event-Record Concepts

Any real adapter must map generator records onto the core input:

- four-momentum as `(E, px, py, pz)` in GeV
- charge in units of the positron charge
- `inOutSign = Incoming` for incoming particles
- `inOutSign = Outgoing` for outgoing particles

For HepMC3, the adapter planning questions are:

- how to identify the two physical incoming beam particles
- how to identify final stable particles
- how to reject documentation, intermediate, beam-remnant, and other
  bookkeeping records
- how to obtain or compute particle charge from the available record fields
- how to verify momentum units and convert to GeV if needed
- how to retain source labels for diagnostics without polluting
  `ChargedParticle`

For Pythia8, the corresponding concepts are:

- incoming beam entries in the Pythia event record
- final-state particle predicates, such as final-state status helpers
- particle charge through Pythia particle data or particle helpers
- four-momentum components in the configured unit convention
- treatment of documentation lines, intermediate states, beam remnants, and
  shower history records

## Build Strategy

Do not add a mandatory generator dependency to the core.

The first implementation now uses an optional adapter target controlled by
`LOWMC_ENABLE_HEPMC3`, which is disabled by default. If HepMC3 is not
available, the default core library, diagnostics, examples, and tests continue
to build.

The adapter should live outside the core target. It may depend on HepMC3, but
`lowmc_core` must not.

To request the optional adapter build, configure with:

```sh
cmake -S . -B build-hepmc3 -DLOWMC_ENABLE_HEPMC3=ON
```

If HepMC3 is installed outside the default CMake search path, pass its prefix:

```sh
cmake -S . -B build-hepmc3 -DLOWMC_ENABLE_HEPMC3=ON -DCMAKE_PREFIX_PATH=/path/to/hepmc3
```

The current placeholder target is `lowmc_hepmc3_adapter`. It validates the
optional package boundary and contains the first one-particle mapping helper,
but it does not yet traverse full HepMC3 events.
With HepMC3 absent, enabling `LOWMC_ENABLE_HEPMC3` should fail at configure time
with a clear message; the default build should remain unaffected.

The local macOS development installation has been verified at `$HOME/hepmc3`,
configured with:

```sh
cmake -S . -B build-hepmc3 -DLOWMC_ENABLE_HEPMC3=ON -DCMAKE_PREFIX_PATH=$HOME/hepmc3
cmake --build build-hepmc3 --target lowmc_hepmc3_adapter
```

This installation exports `HepMC3::HepMC3` and `HepMC3::HepMC3_static`; LowMC
currently links the optional adapter against the shared target
`HepMC3::HepMC3` when it is available.

## Diagnostics Needed First

The first adapter prototype should report:

- selected incoming charged count
- selected outgoing charged count
- ignored neutral candidate count
- signed charge sum `sum_i inOutSign_i Q_i`
- ordinary charge sum `sum_i Q_i`
- rejected charged candidate count
- rejection reasons, at least in a verbose diagnostic path
- source event identifier if available

These diagnostics should be tested on tiny hand-made event records before any
large event sample is used.

## Open Mapping Questions

- Which HepMC3 event samples populate the explicit beam-particle list
  consistently, and which rely only on graph/status conventions?
- Should final stable HepMC3 particles require both `status() == 1` and
  `end_vertex() == nullptr`, or should one condition be sufficient with a
  diagnostic warning for disagreement?
- Are there common HepMC3 samples with unsupported but physically relevant PDG
  ids that require an optional charge override or generated charge table?

## Verified HepMC3 API Notes

The following notes were checked against the HepMC3 upstream source mirror and
should be rechecked against the installed headers when the optional dependency
is added.

- `HepMC3::GenEvent` exposes `particles()` for all particles.
- `HepMC3::GenEvent` exposes `beams()` and `set_beam_particles(...)` for beam
  particles. The first adapter should prefer `event.beams()` for incoming
  particles when it is populated.
- `HepMC3::GenParticle` exposes `pid()`, `abs_pid()`, `status()`,
  `momentum()`, `production_vertex()`, and `end_vertex()`.
- `HepMC3::FourVector` exposes `px()`, `py()`, `pz()`, and `e()`.
- `HepMC3::GenEvent` exposes `momentum_unit()` and `length_unit()`, with
  default momentum unit `Units::GEV`.
- `HepMC3::Units` supports momentum units `MEV` and `GEV`, and
  `Units::convert(...)` converts four-vectors between them. The adapter should
  convert all momenta to GeV before constructing `lowmc::FourVector`.
- `HepMC3::GenParticleData` contains PDG id, status, mass state, mass, and
  momentum, but no charge field. The adapter therefore cannot rely on a generic
  HepMC3 particle-charge member.
- `HepMC3::GenEvent` exposes event weights, optional `GenRunInfo` weight
  names, event/run attributes, and particle/vertex attributes.
- `HepMC3::GenParticle`, `GenVertex`, and `GenEvent` can carry attributes via
  `add_attribute(...)`; attributes include scalar and vector numeric/string
  types such as `DoubleAttribute`, `IntAttribute`, and `StringAttribute`.
- `HepMC3::GenVertex` supports incoming/outgoing particle lists, status,
  position, and attributes. A particle with a production vertex and no end
  vertex can be represented as a final stable particle without making it a
  beam-like root particle.

## First HepMC3 Selection Rule

The first adapter should use this conservative rule unless real sample tests
show a problem:

- Incoming particles: use `event.beams()` and select charged beam particles
  from that list with `inOutSign = Incoming`.
- Outgoing particles: iterate over `event.particles()` and select charged
  particles with `status() == 1` and no `end_vertex()`, assigning
  `inOutSign = Outgoing`.
- Ambiguous particles, such as charged particles with `status() == 1` but an
  end vertex, or particles with no end vertex but non-stable status, should be
  rejected or reported in diagnostics before being considered for inclusion.
- Charge should initially be obtained from a small PDG-id charge helper or a
  user-supplied charge service. Particles whose charge cannot be determined
  should be rejected with an explicit reason.
- Neutral particles should be ignored for the Low-current input, but the first
  adapter diagnostics should count neutral totals separately.

This rule intentionally uses both status and event-graph finality for outgoing
particles. It is stricter than using either alone, but it should catch event
record inconsistencies early.

## Charge Lookup Strategy

The default HepMC3 adapter charge lookup should be based on PDG ids, not on a
Pythia-generated table or mandatory external particle-data dependency.

The lookup should return charge as an integer `charge3 = 3 Q`, where `Q` is the
charge in positron-charge units. The adapter can then pass
`Q = charge3 / 3.0` to `lowmc::ChargedParticle`.

The first strategy is:

- use an explicit table for elementary and special particles:
  - quarks `1..6`
  - charged leptons and neutrinos
  - gluon, photon, `Z`, `W`, Higgs
  - other common special cases when they are intentionally supported
- use PDG digit decoding for ordinary mesons and baryons, based on the standard
  Monte Carlo particle-numbering quark-content fields
- support nuclei only if the PDG ion code `10LZZZAAAI` is implemented
  deliberately; otherwise reject ion codes explicitly at first
- reject private, generator-specific, exotic, or unsupported ids with an
  explicit unknown-charge diagnostic
- allow a user-supplied charge override or callback later, without making it
  necessary for ordinary Standard Model particles

This avoids a large hand-written table while keeping LowMC independent of
Pythia8, ROOT, HepPDT, or another particle-data runtime. A generated table from
Pythia or PDG can still be added later as an optional extension if real HepMC3
samples expose many unsupported ids.

The generator-independent helper lives in `lowmc/PDGCharge.h` as
`chargeTimesThreeFromPDG(int pdgId)`. It returns a result with `known` and
`chargeTimesThree` fields so callers cannot accidentally treat unsupported ids
as neutral. The implementation includes elementary and special particles,
ordinary meson/baryon decoding, and PDG ion-code support for `10LZZZAAAI`.

The ordinary-hadron decoder is covered by focused tests for common charged and
neutral particles, including `pi+`, `pi-`, `pi0`, `K+`, `K-`, `K0`, `D+`,
`D0`, `B+`, `B0`, proton, antiproton, neutron, `Lambda0`, `Sigma+`, `Sigma-`,
`Xi-`, and `Omega-`.

## Dependency-Free Selection Helper

Before writing the HepMC3 adapter itself, LowMC now has a small
dependency-free adapter-selection helper in `lowmc/AdapterSupport.h`.

The helper accepts `AdapterParticleCandidate` records with four-momentum, PDG
id, incoming/final-stable flags, optional beam-remnant flag, and optional
source label. It applies `ParticleSelectionConfig` and
`chargeTimesThreeFromPDG(...)`, returning selected `ChargedParticle` objects
plus `AdapterSelectionDiagnostic`.

Neutral candidates are ignored for the Low-current input and counted as
`ignoredNeutralCandidates`. Unsupported PDG charge is rejected with
`RejectedParticleReason::UnsupportedCharge`. Ambiguous charged candidates are
rejected with `RejectedParticleReason::AmbiguousRecord`. Source labels remain
diagnostic-only.

Focused tests cover incoming beams, final stable outgoing particles, neutral
particles, unsupported PDG ids, excluded beam remnants, disabled incoming-beam
selection, explicit unstable-intermediate inclusion, and ambiguous records.

## Next Implementation Boundary

The optional adapter now provides `makeAdapterParticleCandidate(...)`, which
converts one `HepMC3::GenParticle` plus explicit classification metadata into
one `AdapterParticleCandidate`. It copies the PDG id, stores caller-provided
incoming/final-stable/beam-remnant flags and source label, converts HepMC3
momentum units from `GEV` or `MEV` to GeV, and reorders HepMC3
`(px, py, pz, E)` into LowMC `(E, px, py, pz)`.

The optional adapter also provides `collectAdapterParticleCandidates(...)` for
one `HepMC3::GenEvent`. It uses `event.momentum_unit()`, classifies incoming
beams from `event.beams()`, classifies outgoing final-stable candidates from
`event.particles()` with `status() == 1` and `end_vertex() == nullptr`, and
skips duplicate beam particles during the particle scan. Event-level source
labels currently contain HepMC3 `id`, PDG id, and status.

The hand-built event test records one HepMC3 convention that matters for toy
events: root particles without a production vertex can be returned by
`event.beams()`. Tests that intend ordinary outgoing final particles therefore
attach them to a production vertex.

The next step should combine HepMC3 event candidate collection with the
dependency-free `selectAdapterParticles(...)` helper so one optional adapter
call returns selected charged particles plus diagnostics.

That combined helper now exists as `selectChargedParticles(...)`. It takes a
`HepMC3::GenEvent` plus optional `ParticleSelectionConfig`, uses
`collectAdapterParticleCandidates(...)`, and returns the existing
`AdapterParticleSelectionResult`. The optional tests verify selected incoming
and outgoing charged counts, neutral candidate counting, disabled incoming-beam
selection, rejected charged candidates, and preservation of HepMC3 source
labels.

The next implementation boundary is a tiny HepMC3 file-reading diagnostic that
prints these adapter-selection diagnostics for one or a few events, without
yet connecting HepMC3 input to photon generation.

The diagnostic now exists as `lowmc_hepmc3_adapter_diagnostic` in the optional
HepMC3 build. It uses `HepMC3::ReaderAscii`, accepts `--max-events N`, and
prints one compact line per event with selected incoming/outgoing charged
counts, ignored neutral candidate count, rejected charged candidate count,
signed charge sum, ordinary charge sum, and event momentum unit. It is a
selection/adapter diagnostic only; it does not generate photons.

The diagnostic is now validated on a checked-in tiny HepMC3 ASCII fixture,
`adapters/hepmc3/tests/fixtures/tiny_event.hepmc3`. The fixture was
generated once with `HepMC3::WriterAscii` and is kept static so future changes
exercise file reading rather than writer/reader round-tripping. The CTest
`lowmc_hepmc3_adapter_diagnostic_tiny_fixture` asserts the expected compact
diagnostic line:

- selected incoming charged: `2`
- selected outgoing charged: `1`
- ignored neutral candidates: `2`
- rejected charged candidates: `0`
- signed charge sum: `-1`
- ordinary charge sum: `3`

The diagnostic is also validated on a charge-conserving checked-in fixture,
`adapters/hepmc3/tests/fixtures/charge_conserving_event.hepmc3`. This
fixture uses two incoming protons, two outgoing protons, and one outgoing
neutral photon in `MEV` units. The CTest
`lowmc_hepmc3_adapter_diagnostic_charge_conserving_fixture` asserts:

- selected incoming charged: `2`
- selected outgoing charged: `2`
- ignored neutral candidates: `1`
- rejected charged candidates: `0`
- signed charge sum: `0`
- ordinary charge sum: `4`

The next implementation boundary is a minimal HepMC3-to-weighted-photon
diagnostic that reads this charge-conserving fixture, selects charged
particles, and calls the existing generator-independent weighted event API
without writing photons back into the HepMC3 event record.

That diagnostic now exists as `lowmc_hepmc3_weighted_photon_diagnostic`. It
uses `HepMC3::ReaderAscii`, applies `selectChargedParticles(...)`, reports the
same compact adapter diagnostics, and calls `generateWeightedPhotons(...)` for
events with signed charge sum consistent with zero. Events with nonzero signed
charge sum are reported with `photon_status=skipped_nonconserved_charge`
rather than silently repaired. For generated events, the diagnostic reports
sample count, weight sum, mean weight, weight RMS, effective sample size,
effective sample fraction, and grid-build `k_T`.

The diagnostic is intentionally output-only text. It does not yet write Low
photons into the HepMC3 event record.

The next implementation boundary should define the HepMC3 output convention
for added Low photons before writing a new HepMC3 file.

## HepMC3 Output Writer Contract

The production event-record output should eventually be:

```text
input.hepmc3 -> LowMC afterburner -> output.hepmc3
```

The algorithmic path from weighted samples to unweighted event-record photons
is summarized in `docs/ALGORITHM.md` under "Road To Unweighted Event Output".
HepMC3 writing is allowed only after the count-plus-shape diagnostic gate has
passed for the relevant event class.

The output event should preserve the input event record and append LowMC-added
photons. It should not rewrite the generator history, charged-particle
momenta, or existing vertices in the first implementation.

The first production output convention is settled as:

- photon PDG id: `22`
- photon status: `1`, so downstream tools can treat accepted Low photons as
  stable final-state photons
- photon four-momentum: stored in the event's existing HepMC3 momentum unit;
  LowMC's internal GeV four-vector must be converted before writing
- photon generated mass: leave unset, relying on the massless four-momentum
- vertex placement: no new production vertex in the first durable convention
- topology: preserve all existing vertices and append attributed stable LowMC
  photons as root particles in the event record
- barcode/id policy: let HepMC3 assign ids/barcodes through its normal event
  machinery; do not attempt to preserve candidate indices as HepMC3 ids
- recoil: do not modify any input particle momenta in the first output writer
- event weights: preserve the input event weight vector exactly; do not add a
  new physics weight for unweighted Low photons in the first writer

The first HepMC3 convention stores attributes directly on each added photon.
Attribute names use the `LowMC:` prefix:

- `LowMC:origin = "afterburner"`
- `LowMC:model = "leading-low"`
- `LowMC:recoil = "none"`
- `LowMC:photon_mode = "unweighted"`
- `LowMC:enhancement_factor = 1`
- `LowMC:version`, if available from the build or package metadata
- `LowMC:lambda_event`
- `LowMC:validation_class = "ordinary"` or `"stress"`
- `LowMC:shape_status`
- `LowMC:max_trial_bound_ratio`
- `LowMC:kT`
- `LowMC:eta`
- `LowMC:phi`
- `LowMC:sample_index`, numbered within the LowMC-added photons for that event

The default HepMC3 output writer uses physical multiplicity, equivalent to
enhancement factor `1`. The explicit enhanced mode draws
`N_gamma ~ Poisson(enhancementFactor * lambda_event)` to increase Low-photon
statistics in the event record. Such output is marked on each added photon and
duplicated at event level for events that receive LowMC photons:

- `LowMC:photon_mode = "enhanced"`
- `LowMC:enhancement_factor = F`
- per-photon `LowMC:weight = 1 / F`

Skipped events remain unchanged and receive no LowMC event metadata. If the
enhancement factor is constant for a file, future run-level metadata is the
natural compact duplicate, but this is deferred until a stable file-level
`GenRunInfo` convention is needed. The per-photon `LowMC:weight` remains useful
because downstream histogramming can recover the physical Low yield even
without looking up event metadata.

Directly counting Low photons in enhanced mode gives the enhanced sample yield.
Filling histograms with `LowMC:weight` gives the physical Low yield.

Failure behavior is also settled for the first writer:

- invalid charged-particle selection: write no Low photons for that event and
  report the adapter/yield diagnostic, or fail if the selected user policy is
  strict
- non-conserved signed charge for the selected current: write no Low photons
  and report `skipped_nonconserved_charge`
- failed adaptive-yield validation gate: write no Low photons and report the
  failed gate
- failed shape envelope: write no Low photons for that event
- insufficient accepted shape candidates for the requested Poisson count:
  write no partial event; report the event as incomplete

These policies intentionally prefer missing added photons with diagnostics over
silently writing a biased or ambiguous event.

The weighted generator may use `oversamplingFactor` to produce more weighted
samples per event for histogramming. This does not change the physical event
contribution because each photon weight is scaled down by the same factor. This
weighted diagnostic path remains conceptually distinct from enhanced HepMC3
event samples, even if both use per-photon weights.

An explicit diagnostic-only writer may later store weighted photon samples in
HepMC3 for inspection, but such files must be marked unambiguously:

- event or run attribute `LowMC:photon_mode = "weighted-samples"`
- per-photon attribute `LowMC:weight`
- per-photon attributes for `LowMC:kT`, `LowMC:eta`, `LowMC:phi`, and
  `LowMC:proposal_density` when available
- documentation that these particles are weighted samples, not unweighted
  physical photons

This keeps the normal output format safe for downstream analyses while still
leaving room for HepMC3-based diagnostics.

The first implementation step is complete as an in-memory optional HepMC3
adapter helper, `appendUnweightedLowPhotons(...)`. It appends accepted
shape-unweighted candidates to a `HepMC3::GenEvent`, converts momenta from
LowMC GeV to the event's momentum unit, attaches the initial `LowMC:`
attributes to the added photons, preserves event weights and original
particles, and rejects invalid enhancement factors before modifying the event.

The first file-to-file diagnostic,
`lowmc_hepmc3_afterburner_diagnostic`, now reads HepMC3 ASCII input, runs the
count-plus-shape path for charge-conserving selected events, appends unweighted
Low photons through `appendUnweightedLowPhotons(...)`, and writes a new HepMC3
ASCII file. Its default failed-event policy is `skip`: failed events are
written unchanged with diagnostics. `--failure-policy fail` is available for
strict runs.

The attributed-root-photon topology is accepted as the first durable HepMC3
convention. The originally preferred no-incoming afterburner vertex is deferred
because it did not round-trip cleanly through the local HepMC3 ASCII test. The
root-photon convention is deliberately conservative: it makes added LowMC
photons visible as ordinary stable photons, stores all LowMC provenance on the
new photons, and leaves the original event graph untouched.

## HepMC3 Output Open Questions

- Do important downstream analysis chains require a custom non-`1` status code
  for afterburner-added photons, despite the first diagnostic using status `1`?
- Should a diagnostic-only weighted HepMC3 writer use `LowMC:` attributes or a
  separate namespace to make accidental analysis use even less likely?
- Should the production command default remain `skip`, or should production
  runs use `fail` once validation gates are finalized?
- Should a later optional mode expose a tested afterburner-vertex topology for
  users who prefer an explicit graph node, if a portable convention is found?

## Direct Pythia In-Memory Adapter Plan

The next adapter target is a direct Pythia8 interface that works on an
in-memory `Pythia8::Event`. This is separate from the existing Pythia-to-HepMC3
validation path. It should be useful for users who run Pythia and fill
histograms in the same executable, without writing an intermediate HepMC3 file.

The planned layering is:

```text
Pythia8::Event
-> optional lowmc_pythia8 adapter
-> AdapterParticleCandidate vector
-> selectAdapterParticles(...)
-> LowMC core
```

The Pythia adapter should live in an optional target, for example under
`adapters/pythia8`, controlled by a CMake option such as
`LOWMC_ENABLE_PYTHIA8`. The default build must continue to build the core,
diagnostics, examples, and tests without Pythia.

The first adapter boundary should not create a new core input type. It should
reuse:

- `AdapterParticleCandidate`
- `ParticleSelectionConfig`
- `AdapterParticleSelectionResult`
- `selectAdapterParticles(...)`
- `chargeTimesThreeFromPDG(...)`

The adapter-specific work is only classification and conversion:

- identify incoming beams according to Pythia's event-record convention
- identify final stable particles with Pythia final-state helpers
- exclude documentation lines, beam remnants, intermediate shower particles,
  and decayed particles by default
- copy PDG id into the candidate
- convert Pythia momentum to `(E, px, py, pz)` in GeV
- attach source labels containing Pythia particle index, PDG id, status, and
  classification

The default selected current should match the HepMC3/default adapter policy:

- include incoming beams by default
- include final stable charged particles by default
- ignore neutral candidates after charge lookup, while counting them in
  diagnostics
- reject unsupported-charge candidates explicitly
- reject ambiguous charged candidates rather than guessing

Charge lookup should use LowMC's PDG helper by default. This keeps one charge
convention across HepMC3, direct Pythia, and direct/manual core input. Pythia's
particle-data charge is still valuable as an optional diagnostic comparison,
especially for unusual ids or user-defined particles, but it should not be the
first default charge source.

The first direct Pythia diagnostic stops at selected-particle diagnostics. It
can instantiate either a deterministic toy event or a small generated SoftQCD
event, run the adapter selection, and print the same compact fields as the
HepMC3 adapter diagnostic:

- selected incoming charged count
- selected outgoing charged count
- ignored neutral candidate count
- rejected charged candidate count
- signed charge sum
- ordinary charge sum
- event number or Pythia event label if available

This diagnostic now exists as `lowmc_pythia8_adapter_diagnostic` in the
optional Pythia8 build. The optional target is controlled by
`LOWMC_ENABLE_PYTHIA8`, disabled by default, and uses `pythia8-config` to find
headers, libraries, and the XML documentation path. The default LowMC build
remains independent of Pythia8.

The first mapping uses Pythia event-record indices `1` and `2` as incoming
beam candidates and all non-beam particles with `particle.isFinal()` as final
stable outgoing candidates. This is intentionally conservative and compact for
the first diagnostic; future work should add consistency diagnostics if real
Pythia configurations require a more nuanced beam/remnant treatment.

The deterministic toy diagnostic gives:

```text
event_index=0 mode=toy selected_incoming_charged=2 selected_outgoing_charged=2 ignored_neutral_candidates=1 rejected_charged_candidates=0 signed_charge_sum=0 ordinary_charge_sum=4
```

A one-event SoftQCD diagnostic with seed `12345` has also been exercised. It
selects two incoming charged beam particles, preserves signed charge
conservation, and reports compact adapter diagnostics without writing any
event record.

Once this selection diagnostic is stable, a second diagnostic can call
`generateWeightedPhotons(...)` on the selected charged particles and report the
existing weighted photon summary. HepMC3 output writing is not part of the
direct Pythia boundary.

That second diagnostic now exists as
`lowmc_pythia8_weighted_photon_diagnostic`. It is a separate executable rather
than an overloaded mode of `lowmc_pythia8_adapter_diagnostic`, so the simple
selection diagnostic remains small and unambiguous. The weighted diagnostic
supports deterministic toy input, a deliberately non-charge-conserving toy
input, and generated SoftQCD events. For charge-conserving selected currents it
calls `generateWeightedPhotons(...)` and reports sample count, base sample
count, oversampling factor, weight sum, mean weight, weight RMS, ESS, ESS/N,
and grid-build `k_T`. For non-conserved selected currents it reports
`photon_status=skipped_nonconserved_charge` and writes no photons.

Representative checked outputs are:

```text
event_index=0 mode=toy selected_incoming_charged=2 selected_outgoing_charged=2 ignored_neutral_candidates=1 rejected_charged_candidates=0 signed_charge_sum=0 ordinary_charge_sum=4 photon_status=generated samples=8 base_samples=8 oversampling_factor=1 sum_weight=1.47073 mean_weight=0.183841 weight_rms=0.00589109 effective_sample_size=7.99179 effective_sample_fraction=0.998974 grid_build_kt=0.02
```

```text
event_index=0 mode=toy-nonconserved selected_incoming_charged=2 selected_outgoing_charged=2 ignored_neutral_candidates=1 rejected_charged_candidates=0 signed_charge_sum=-2 ordinary_charge_sum=2 photon_status=skipped_nonconserved_charge
```

```text
event_index=0 mode=softqcd selected_incoming_charged=2 selected_outgoing_charged=46 ignored_neutral_candidates=25 rejected_charged_candidates=2 signed_charge_sum=0 ordinary_charge_sum=4 photon_status=generated samples=4 base_samples=4 oversampling_factor=1 sum_weight=3.59379 mean_weight=0.898448 weight_rms=0.0785295 effective_sample_size=3.96967 effective_sample_fraction=0.992418 grid_build_kt=0.02
```

### Direct Pythia Open Questions

- Should incoming beams remain fixed to event-record indices `1` and `2`, or
  should the adapter add status/history consistency diagnostics for unusual
  Pythia configurations?
- Is `particle.isFinal()` sufficient for the default final-stable selection, or
  should a later adapter also inspect status-code ranges for diagnostics?
- How should beam remnants be recognized and reported?
- Should a later diagnostic compare `chargeTimesThreeFromPDG(...)` with
  Pythia particle-data charge and report disagreements?
- How closely do direct Pythia selection diagnostics match the Pythia-to-HepMC3
  path for the same generated event?

## Direct Pythia Versus Pythia-To-HepMC3 Cross-Check

The first direct-Pythia/HepMC3 comparison diagnostic now exists as
`lowmc_pythia8_hepmc3_crosscheck`. It is built only when both optional adapter
paths are enabled:

```sh
cmake -S . -B build-pythia8-hepmc3 \
  -DLOWMC_ENABLE_PYTHIA8=ON \
  -DLOWMC_ENABLE_HEPMC3=ON \
  -DCMAKE_PREFIX_PATH=$HOME/hepmc3
```

The diagnostic generates a deterministic Pythia8 SoftQCD event in memory,
selects charged particles with the direct Pythia adapter, converts the same
event to HepMC3 in memory through `Pythia8ToHepMC3`, and selects charged
particles with the HepMC3 adapter. It compares:

- selected incoming charged count
- selected outgoing charged count
- ignored neutral candidate count
- rejected charged candidate count
- signed charge sum
- ordinary charge sum

The first checked seed, `12345`, matches exactly:

```text
event_index=0 status=match direct_selected_incoming_charged=2 hepmc3_selected_incoming_charged=2 direct_selected_outgoing_charged=46 hepmc3_selected_outgoing_charged=46 direct_ignored_neutral_candidates=25 hepmc3_ignored_neutral_candidates=25 direct_rejected_charged_candidates=2 hepmc3_rejected_charged_candidates=2 direct_signed_charge_sum=0 hepmc3_signed_charge_sum=0 direct_ordinary_charge_sum=4 hepmc3_ordinary_charge_sum=4
```

The cross-check now also compares weighted photon summaries after verifying
that the selected-particle diagnostics match. Both paths call
`generateWeightedPhotons(...)` with the same phase space, tuned grid
configuration, sample count, oversampling factor, and photon seed. It compares
sample count, weight sum, mean weight, weight RMS, ESS, ESS/N, and grid-build
`k_T` with a tight floating-point tolerance.

For the first checked seed, `12345`, the weighted summaries also match:

```text
event_index=0 selection_status=match direct_selected_incoming_charged=2 hepmc3_selected_incoming_charged=2 direct_selected_outgoing_charged=46 hepmc3_selected_outgoing_charged=46 direct_ignored_neutral_candidates=25 hepmc3_ignored_neutral_candidates=25 direct_rejected_charged_candidates=2 hepmc3_rejected_charged_candidates=2 direct_signed_charge_sum=0 hepmc3_signed_charge_sum=0 direct_ordinary_charge_sum=4 hepmc3_ordinary_charge_sum=4 weighted_status=match direct_samples=8 hepmc3_samples=8 direct_sum_weight=7.15094 hepmc3_sum_weight=7.15094 direct_mean_weight=0.893867 hepmc3_mean_weight=0.893867 direct_weight_rms=0.0934473 hepmc3_weight_rms=0.0934473 direct_effective_sample_size=7.91351 hepmc3_effective_sample_size=7.91351 direct_effective_sample_fraction=0.989189 hepmc3_effective_sample_fraction=0.989189 direct_grid_build_kt=0.02 hepmc3_grid_build_kt=0.02
```

No expected differences are currently documented for this event. If later
Pythia configurations differ after HepMC3 conversion, the cross-check should
print the mismatch rather than silently accepting it.
