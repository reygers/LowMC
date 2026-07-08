# STATUS

## Current State

LowMC is in early development.

The current project definition is a generator-independent C++ afterburner for
coherent soft-photon emission based on the leading Low theorem.

The initial generator-independent core kernel skeleton exists under `src/`.
It implements four-vector arithmetic, charged particles with `inOutSign`,
photons with four-momentum and weight, Low-current evaluation, intensity
evaluation, a `lowPhotonDensity` helper for
`alpha / (2 pi)^2 * S_event(k)`, charge-conservation diagnostics, and a
`PhotonPhaseSpace` value type with validation for positive, non-empty `k_T`,
`eta`, and `phi` ranges. The first weighted-photon interface can convert an
explicit proposal point and proposal density into a `Photon` with
`weight = lowPhotonDensity / q_mu`.
The documentation now defines the intended weighted-photon convention in terms
of the physical measure `k_T dk_T d eta d phi`, including a logarithmic default
proposal for the one-dimensional `k_T` factor and the planned fixed-`k_T`
`eta-phi` proposal interface. The core includes log-uniform `k_T` proposal
utilities for mapping unit random values to `k_T` and evaluating the proposal
density in both `dk_T` and `k_T dk_T` conventions. The core also includes a
uniform fixed-`k_T` `eta-phi` proposal baseline and a fixed-`k_T` weighted toy
sampler for diagnostics. A first adaptive fixed-`k_T` `eta-phi` proposal now
exists as a mixture of a uniform floor and local boxes around charged-particle
directions; it is available in the weighted toy diagnostics through
`--proposal adaptive`. A minimal full phase-space weighted path now combines
log-uniform `k_T` sampling with uniform `eta-phi` sampling and weights photons
with the full proposal density `q_mu = q_kT,mu q_eta_phi`. Basic sampler
diagnostic summaries report weight sums, ranges, variance, and effective sample
size for weighted samples.

A two-particle `eta-phi` diagnostic also exists. It evaluates `S(k)` at fixed
`k_T` for a transverse back-to-back charged-particle toy event and writes CSV
output. The default toy momentum is deliberately modest so the finite-mass
dead-cone scale is visible. The toy-particle momentum, direction, and grid
acceptance can also be configured for boosted and forward stress-test cases.
The diagnostic grid resolution is configurable through `--eta-bins` and
`--phi-bins`; defaults remain `101 x 121`.
The diagnostic now routes its acceptance bounds through `PhotonPhaseSpace`
validation before writing output. An optional Python/matplotlib script can plot
that CSV as a PNG heatmap. A second diagnostic executable writes weighted
fixed-`k_T` photon samples from the same back-to-back toy event and reports
basic weight-summary diagnostics. A toy sampler-study executable writes
weight-summary rows for several charged-pion momenta so the uniform proposal
baseline can be inspected as boosts increase. A comparison script integrates
the fixed-grid diagnostic and compares it with the weighted-sample mean. For
the default transverse two-pion toy, a 100000-sample weighted run agrees with
the fixed-grid integral within one standard error. An optional plotting script
turns sampler-study summary CSV output into a compact PNG showing mean weight,
weight RMS, and effective sample fraction versus charged-pion momentum. A
grid-vs-weighted stress-study helper now runs boosted and forward two-particle
toy cases and summarizes grid integrals, weighted means, pulls, and effective
sample fractions. The first default stress run agrees for the `p = 1 GeV`
midrapidity case, but boosted midrapidity cases already show poor effective
sample fractions and fixed-grid disagreement; a `p = 10 GeV` high-stat check
indicates the current coarse grid is under-resolved as a reference for narrow
structures. Increasing the `p = 10 GeV` midrapidity grid to `401 x 481` moves
the grid integral into agreement with the weighted estimate, while the
effective sample fraction remains below one percent.
An initial `100000`-sample sampler-study comparison shows that the adaptive
proposal improves effective sample fraction relative to the uniform proposal in
the transverse two-pion toy, but high-momentum cases remain noisy and require
further tuning. The adaptive proposal parameters are now exposed in the
weighted diagnostics, and `scripts/scan_adaptive_tuning.py` writes a compact
uniform-vs-adaptive scan CSV with weight RMS, effective sample fraction, and
fine-grid pull columns. In a first `50000`-sample midrapidity scan, the best
`p = 10 GeV` tested variant used `widthScale = 20` and `minHalfWidth = 0.02`,
reaching `ESS/N ~= 0.107`; `p = 100` and `1000 GeV` remained much less stable.
The tuning scan can now be plotted and summarized with
`scripts/plot_adaptive_tuning.py`. The generated summary favors
`adaptive_w20_min0.02` for `p = 10 GeV`, a narrower `adaptive_w2_min0.01` for
`p = 100 GeV`, and the uniform baseline for `p = 1000 GeV` under the current
`|pull| <= 3` sanity cut, reinforcing that the local-box proposal is only a
diagnostic baseline.

A focused `p = 100 GeV` check is now available through
`scripts/run_p100_focus.py`. In the default run, grid integrals continue to
increase from `401 x 481` through `1601 x 1921`; the finest grid gives
`148.78`, with a remaining `~3%` change from the previous grid. Wider local
boxes substantially improve the sampler relative to uniform sampling. The best
tested stable variant is `adaptive_w40_min0.1`, with `ESS/N ~= 0.029` and
`|pull| ~= 2.7` against the finest available grid. This supports using
`p = 100 GeV` as the main high-energy tuning benchmark and keeping
`p = 1000 GeV` as an extreme stress case.

A first grid-adaptive fixed-`k_T` `eta-phi` proposal prototype now exists. It
builds a cell map from center-sampled `S(k)`, refines high-mass cells, includes
a uniform floor, and still computes exact Low weights at sampled points. In the
focused `p = 100 GeV` study, a `24 x 24` grid with three refinement levels
reaches `ESS/N ~= 0.095`, substantially better than local boxes, but its
current pull against the finest available grid is about `4.7`; this should be
treated as a diagnostic warning rather than a production default.

`scripts/scan_grid_tuning.py` now scans grid-adaptive map parameters. In the
first `100000`-sample `p = 100 GeV` scan, aggressive `32 x 32` maps with three
refinement levels reach `ESS/N ~= 0.16`. When interpreted against the finest
available fixed grid alone, these maps appeared to have large pulls. After
retabulating the same scan against the high-stat weighted reference with
`scripts/retabulate_grid_tuning.py`, the best tested map,
`grid_32x32_l3_r0.1_u0.1`, has `ESS/N ~= 0.164`, weight RMS `~= 353`, and
pull `~= -1.1` against `157.54 +/- 0.32`. The aggressive refined grid maps now
look consistent with the updated reference; the `1601 x 1921` fixed grid should
be treated as a deterministic lower cross-check for this benchmark, not as the
sole tuning reference.

`scripts/run_p100_reference.py` now compares fixed global grids with high-stat
weighted grid-adaptive references. With `1000000` samples, two independent
grid-adaptive weighted estimates give `157.72 +/- 0.70` and `157.49 +/- 0.36`,
well above the `1601 x 1921` fixed-grid integral `148.78`. This strongly
indicates that the fixed global grid is still systematically low for the
`p = 100 GeV` benchmark. For proposal tuning, the high-stat weighted
grid-adaptive reference should be used alongside the finest fixed grid, not
replaced by the fixed grid alone.

A first smooth local-kernel fixed-`k_T` proposal now exists as a diagnostic
alternative to local boxes. It uses compact triangular kernels in `eta` and
`phi` around charged-particle directions, with analytic truncation to the
acceptance and a uniform floor. The `lowmc_weighted_sampler_study` and
`lowmc_fixed_kt_weighted_samples` diagnostics accept `--proposal smooth`.
`scripts/compare_p100_proposals.py` compares the tuned grid map, the broad
local-box baseline, and the first smooth triangular kernel against the
high-stat weighted reference. With `100000` samples at `p = 100 GeV`, the tuned
grid map remains best (`ESS/N ~= 0.160`, weight RMS `~= 362`, pull `~= 0.39`).
The smooth triangular kernel improves over the broad local-box baseline
(`ESS/N ~= 0.040` vs `0.029`, RMS `~= 749` vs `914`) but is still well below
the grid-adaptive map. The first smooth result is statistically compatible
with the weighted reference within about `1.6` combined standard errors.

`scripts/scan_smooth_tuning.py` now scans the smooth triangular-kernel
parameters for the same `p = 100 GeV` benchmark while including the tuned grid
map and broad local-box baseline as fixed comparison rows. In the first
`100000`-sample scan, the best smooth setting is
`smooth_triangular_w40_min0.1_u0.4`, with `ESS/N ~= 0.066`, weight RMS
`~= 583`, and pull `~= -1.6` against the weighted reference. This is about
`2.3` times the broad local-box baseline in effective sample fraction, but only
about `0.41` times the tuned grid-adaptive benchmark. The smooth local-kernel
idea is therefore useful as a better local diagnostic baseline, but it does not
currently challenge the grid-adaptive map for the high-energy benchmark.

A first variance-aware grid-refinement mode now exists for the grid-adaptive
diagnostic proposal. The mode probes multiple points per candidate cell,
estimates `Var[S(k)]`, and refines cells according to `area^2 Var[S(k)]`
instead of only center-sampled cell mass. The sampler-study diagnostic exposes
this through `--grid-refinement-score variance` and
`--grid-variance-probe-bins`. In the first `100000`-sample `p = 100 GeV`
comparison, the three-probe-bin variance grid reaches `ESS/N ~= 0.167`, weight
RMS `~= 347`, and pull `~= -1.7`, compared with `ESS/N ~= 0.160`, RMS `~= 362`,
and pull `~= 0.39` for the current mass-refined benchmark. This is competitive
but not yet a decisive improvement. A five-seed check with `100000` samples per
seed gives `ESS/N = 0.1628 +/- 0.0046` for the mass-refined grid and
`ESS/N = 0.1624 +/- 0.0045` for the three-probe-bin variance-refined grid. The
two modes are statistically tied in this benchmark. Since mass refinement is
cheaper and simpler, it remains the preferred diagnostic default for now; the
variance-aware mode is kept as an exploratory option.

The tuned mass-refined grid proposal is now connected to a full phase-space
weighted diagnostic path. `lowmc_full_weighted_sampler_study` samples `k_T`
log-uniformly, samples `eta-phi` either uniformly or from the tuned
grid-adaptive map, combines proposal densities as
`q_mu = q_kT,mu q_eta_phi`, and evaluates exact Low weights at every sampled
point. For the first diagnostic, the grid map is built once at the geometric
mean of the configured `k_T` range. In a `p = 100 GeV`, `k_T = 0.002..0.2 GeV`,
`100000`-sample comparison, the full grid path gives mean weight `~= 0.287`,
RMS `~= 0.643`, and `ESS/N ~= 0.167`; the full uniform path gives mean weight
`~= 0.365`, RMS `~= 31.6`, and `ESS/N ~= 1.3e-4`. The grid result is compatible
with the fixed-`k_T` angular reference scaled by `k_T^2 ln(kTMax/kTMin)`, which
predicts about `0.290`.

`scripts/multirange_full_weighted.py` now checks the full grid-adaptive path
over multiple seeds and `k_T` ranges. With five seeds and `100000` samples per
seed, the ranges `0.001..0.1 GeV` and `0.002..0.2 GeV` both give mean weight
`0.2896 +/- 0.0026`, as expected because they have the same logarithmic width.
The one-decade range `0.02..0.2 GeV` gives `0.1448 +/- 0.0013`, half the
two-decade yield, again matching the expected `ln(kTMax/kTMin)` scaling. The
mean/reference ratio is `0.998 +/- 0.009` for all three ranges, and `ESS/N`
stays near `0.160`. This supports building the grid at the geometric-mean
`k_T` for the first diagnostic path while evaluating exact Low weights at the
sampled `k_T`.

A minimal generator-independent event-level weighted API now exists. The core
function `generateWeightedPhotons` accepts already selected `ChargedParticle`
objects plus `WeightedEventGeneratorConfig`, builds the coherent current,
samples log-uniform `k_T` and tuned mass-refined grid-adaptive `eta-phi`, and
returns weighted photon sample records. Each record retains the `Photon`, the
sampled `k_T`, `eta`, `phi`, and the full proposal density `q_mu`; the result
also contains a convenience `photons` vector, a weight summary,
charge-conservation diagnostic, and the grid-build `k_T`. The number of
samples, integer `oversamplingFactor`, and random seed are explicit inputs.
Oversampling returns more weighted samples and scales each photon weight down
by the same factor, preserving the expected event contribution for histograms.
This deliberately does not include HepMC3, Pythia, recoil, unweighting, or
particle-selection policy. The
`lowmc_weighted_event_example` diagnostic exercises this API on the two-pion
toy event; with `1000` samples at `p = 100 GeV` and `k_T = 0.002..0.2 GeV`, it
gives mean weight `~= 0.280` and effective sample size `~= 176`. A compact
usage example is in `docs/API_EXAMPLE.md`.

A fixed-`k_T` visual sampler-validation diagnostic now compares the
two-particle grid shape, tuned grid-adaptive proposal occupancy, and weighted
sample shape in one figure. It overlays normalized one-dimensional `eta` and
`phi` projections and adds zoomed projection panels around the strongest grid
peak so sharp boosted structures are easier to compare by eye. First runs with
`401 x 481` truth grids and `100000` tuned-grid samples show good visual
agreement for `p = 1 GeV` and `p = 10 GeV`. At `p = 100 GeV`, the diagnostic
exposes the known resolution problem: the fixed grid gives a substantially
lower integral than the weighted reference and the projections are sensitive to
binning. The plotter can now label an external reference integral and print an
explicit warning, so `p = 100 GeV` is presented as a stress visualization, not
as a fixed-grid truth integral at that resolution.
A focused `p = 10 GeV` follow-up separates weighted-sample statistics from
grid resolution. Increasing the tuned-grid sample from `100000` to `1000000`
does not remove the visible local zoom differences when the `401 x 481` grid is
kept fixed; the coarse-grid integral is low by about one percent and gives a
large pull. Replacing the reference with an `801 x 961` grid gives
`103.657` from the grid and `103.694 +/- 0.034` from the weighted sample, a
pull of about `1.1`. Thus the remaining `p = 10 GeV` blue-green differences in
the zoomed projections are primarily grid/binning effects, not a problem with
the weighted normalization.

`docs/ADAPTER_CONTRACT.md` now defines the boundary between generator-specific
event records and the LowMC core. The core input remains a cleaned list of
charged particles with `(E, px, py, pz)`, charge in positron-charge units, and
`inOutSign`. Adapters own HepMC3/Pythia8 status interpretation, default to
incoming beams plus final stable charged particles, and should report
selection and charge diagnostics. Generator-specific status codes and event
history traversal remain out of scope for the core.
The first generator-independent adapter-support types now exist in
`lowmc/AdapterSupport.h`. They define `ParticleSelectionConfig` with incoming
beams and final stable charged particles enabled by default, plus
`AdapterSelectionDiagnostic` for selected incoming/outgoing charged counts,
signed charge sum, ordinary charge sum, and rejected charged candidates. Focused
core tests cover the default selection options, option toggling, and diagnostic
summary calculations.
`docs/API_EXAMPLE.md` now presents the direct `ChargedParticle` path as the
primary manual interface for hand-made toy events, including a two-pion example,
photon phase-space configuration, weighted result diagnostics, and the existing
CSV-writing diagnostic executable.
The same pattern is now also available as the compiled
`lowmc_direct_charged_particles_example`, which prints a compact summary and the
first weighted photon without command-line or file-I/O complexity.
`docs/ADAPTER_PLAN.md` now records the first real event-record adapter planning
pass. HepMC3 is chosen as the first adapter target because it is a portable
event-record interface, while Pythia8 remains a later thin generator interface.
The plan keeps generator dependencies optional and outside `lowmc_core`, and it
lists the exact HepMC3 mapping questions that must be answered before
implementation: incoming beams, final stable particles, charge lookup, momentum
units, ambiguous records, and diagnostics.
The HepMC3 mapping has now been checked against upstream source headers. The
first adapter should prefer `GenEvent::beams()` for incoming particles, use
`GenParticle::status() == 1` together with `end_vertex() == nullptr` for final
stable outgoing particles, convert `GenParticle::momentum()` from the event
`momentum_unit()` to GeV, and obtain charge from a PDG-id lookup or
user-supplied charge service because bare HepMC3 particle data does not carry a
generic charge member.
The charge lookup strategy is now settled for the first adapter: use a
self-contained PDG-id lookup returning integer `3Q`, with an explicit
elementary/special-particle table, algorithmic ordinary meson/baryon decoding,
deliberate ion handling or explicit ion rejection, and diagnostics for
unsupported ids. Pythia-derived or generated charge tables remain possible
future extensions, not first-adapter requirements.
The generator-independent implementation now exists in `lowmc/PDGCharge.h`.
`chargeTimesThreeFromPDG(int pdgId)` returns a result with `known` and
`chargeTimesThree`, supports elementary/special particles, ordinary
meson/baryon decoding, and PDG ion codes, and returns unsupported for private,
generator-specific, invalid, or otherwise unknown ids. Focused tests cover
common mesons, baryons, antiparticles, ions, and unsupported ids.
The PDG charge helper is now connected to a dependency-free adapter-selection
exercise in `lowmc/AdapterSupport.h`. `AdapterParticleCandidate` represents a
small generator-independent record with four-momentum, PDG id,
incoming/final-stable flags, optional beam-remnant flag, and optional source
label. `selectAdapterParticles(...)` applies `ParticleSelectionConfig` and
`chargeTimesThreeFromPDG(...)`, returns cleaned `ChargedParticle` objects, and
fills `AdapterSelectionDiagnostic`. Neutral candidates are ignored for the
Low-current input but counted as `ignoredNeutralCandidates`; unsupported PDG
charge is rejected as `UnsupportedCharge`; ambiguous charged records are
rejected as `AmbiguousRecord`. Non-incoming, non-final unstable intermediates
remain excluded by default, but the dependency-free helper can include them as
outgoing candidates when explicitly requested. Focused tests cover the default
policy, disabled incoming-beam selection, explicit unstable-intermediate
inclusion, neutral counting, unsupported charge, beam remnants, unstable
intermediates, and ambiguous records.
The optional HepMC3 build boundary now exists without making HepMC3 mandatory.
`LOWMC_ENABLE_HEPMC3` defaults to `OFF`. When enabled, CMake looks for a HepMC3
package configuration and builds the placeholder `lowmc_hepmc3_adapter` target
under `adapters/hepmc3`; when HepMC3 is absent, the enabled configure path
fails with a clear message. The default build remains independent of HepMC3 and
still builds `lowmc_core`, diagnostics, examples, and tests. The local macOS
source installation at `$HOME/hepmc3` has been verified. It exports
`HepMC3::HepMC3` and `HepMC3::HepMC3_static`, and the placeholder
`lowmc_hepmc3_adapter` target builds successfully against `HepMC3::HepMC3`.
The optional adapter now includes the first tiny HepMC3 mapping function:
`makeAdapterParticleCandidate(...)` converts one `HepMC3::GenParticle` plus an
explicit `HepMC3ParticleClassification` into `AdapterParticleCandidate`. It
copies the PDG id, classification flags, and source label; converts `GEV` and
`MEV` HepMC3 momenta to GeV; and maps HepMC3 `(px, py, pz, E)` ordering to
LowMC `(E, px, py, pz)`. The HepMC3-enabled build includes
`lowmc_hepmc3_adapter_tests` for this conversion.
The optional adapter now also collects candidates from one `HepMC3::GenEvent`
through `collectAdapterParticleCandidates(...)`. It uses the event momentum
unit, emits incoming beam candidates from `event.beams()`, emits outgoing
final-stable candidates with `status() == 1` and no end vertex, and avoids
duplicating beams during the particle scan. Source labels currently include
HepMC3 `id`, PDG id, and status. The optional test builds a small HepMC3 graph
with beams, final stable particles, a neutral final particle, and a decaying
status-1 particle whose final daughter is retained.
The optional adapter now connects event candidate collection to the
dependency-free selection helper through `selectChargedParticles(...)`, which
returns `AdapterParticleSelectionResult` for one `HepMC3::GenEvent`. The
optional tests cover the default selection diagnostics, neutral counting,
disabled incoming-beam selection, rejected charged beam diagnostics, and
preservation of HepMC3 source labels.
The HepMC3-enabled build now also includes
`lowmc_hepmc3_adapter_diagnostic`. It reads HepMC3 ASCII files with
`HepMC3::ReaderAscii`, uses `selectChargedParticles(...)`, and prints compact
adapter-selection diagnostics for one event by default or multiple events with
`--max-events`. The diagnostic is now validated on the checked-in
`tiny_event.hepmc3` fixture, with CTest asserting selected incoming/outgoing
charged counts, ignored neutral candidates, rejected charged candidates,
signed charge sum, and ordinary charge sum. The fixture deliberately has
nonzero signed charge sum, so it is useful for adapter diagnostics but not yet
ideal as the first photon-generation input.
The optional diagnostic is now also validated on the checked-in
`charge_conserving_event.hepmc3` fixture. It contains two incoming protons, two
outgoing protons, and one outgoing neutral photon in `MEV` units. The
corresponding CTest asserts two selected incoming charged particles, two
selected outgoing charged particles, one ignored neutral candidate, no rejected
charged candidates, signed charge sum `0`, and ordinary charge sum `4`. This
fixture is the first clean input for a HepMC3-to-weighted-photon diagnostic.
The HepMC3-enabled build now also includes
`lowmc_hepmc3_weighted_photon_diagnostic`. It reads HepMC3 ASCII files,
selects charged particles through the adapter, skips photon generation for
events with nonzero signed charge sum, and otherwise calls the
generator-independent `generateWeightedPhotons(...)` path. The output is one
compact text line per event with adapter diagnostics plus the weighted-photon
summary. It does not write photons back into the HepMC3 event record. A smoke
CTest validates the generated path on `charge_conserving_event.hepmc3`.
The first HepMC3 output-writer contract is now documented in
`docs/ADAPTER_CONTRACT.md` and `docs/ADAPTER_PLAN.md`. The current
implementation writes physical output with enhancement factor `1` by default
and supports explicit enhanced output through `--enhancement-factor F`. Added
Low photons should be appended as stable PDG `22`, status-`1` particles. The writer should preserve
the input event record, existing event weights, charged-particle momenta, and
vertices, and it should make the no-recoil policy visible through `LowMC:`
provenance attributes. Failure modes such as invalid selection,
non-conserved signed charge, failed yield gates, shape-envelope failures, or
insufficient accepted candidates should write no Low photons for the event or
fail according to explicit user policy; they should not produce partial
event-record output. The originally planned
dedicated no-incoming afterburner vertex did not round-trip cleanly through the
local HepMC3 ASCII reader/writer, so the first durable HepMC3 topology is
attributed stable root photons: provenance is stored on the added photons and
existing vertices remain unchanged.

The first optional HepMC3 output-writer prototype now exists as an in-memory
adapter helper, `appendUnweightedLowPhotons(...)`. It appends already accepted
shape-unweighted Low photon candidates to a `HepMC3::GenEvent` as stable PDG
`22`, status-`1` particles. It converts LowMC's internal GeV four-momenta to
the event's HepMC3 momentum unit, preserves the original event record and event
weights, applies no recoil, and attaches the first `LowMC:` provenance
attributes directly to the added photons. The helper validates that the
enhancement factor is positive and finite before modifying the event. Focused
optional tests cover MeV unit conversion, no-recoil preservation, event-weight
preservation, photon attributes, final-state photon representation, invalid
enhancement rejection, and enhanced-mode photon weights.

The first file-to-file HepMC3 afterburner diagnostic now exists as
`lowmc_hepmc3_afterburner_diagnostic`. It reads HepMC3 ASCII input, selects
charged particles through the adapter, skips non-charge-conserving or otherwise
failed events by default, runs the current adaptive-yield plus shape-unweighting
path for accepted events, appends ordinary unweighted Low photons with
`appendUnweightedLowPhotons(...)`, and writes a new HepMC3 ASCII file. The
strict alternative is available through `--failure-policy fail`. Optional tests
verify invalid-enhancement rejection, a charge-conserving physical fixture
round-trip in which the output contains the original particles plus one
attributed LowMC photon, and an enhanced fixture with four attributed LowMC
photons carrying `LowMC:weight = 0.1`.

The HepMC3 topology decision is now settled for the first production-oriented
path. Attributed stable root photons are accepted because they round-trip
cleanly through HepMC3 ASCII, are visible as ordinary stable photons to
downstream analyses, and avoid implying a recoil or production-vertex model.
A dedicated afterburner-vertex topology may be revisited only as a later
optional mode if a portable convention is found and tested.

`docs/USAGE.md` now documents the first HepMC3 afterburner workflow for users:
how to build the optional HepMC3 target, how to run
`lowmc_hepmc3_afterburner_diagnostic`, what the output file contains, the
attributed-root-photon topology, the default `skip` failed-event policy, the
strict `--failure-policy fail` alternative, and the default physical
enhancement factor `1`.

Enhanced HepMC3 output is now implemented in the afterburner diagnostic through
`--enhancement-factor F`. Default physical output remains enhancement factor
`1` and samples `N_gamma ~ Poisson(lambda_event)`. Enhanced output samples
`N_gamma ~ Poisson(F * lambda_event)`. The first implementation stores
metadata directly on each added Low photon:
`LowMC:photon_mode = "enhanced"`, `LowMC:enhancement_factor = F`, and
`LowMC:weight = 1 / F`. Direct counting gives the enhanced sample yield;
histograms filled with `LowMC:weight` recover the physical Low yield. Event-level
metadata is now also written as a compact duplicate on events that receive LowMC
photons; run-level metadata remains deferred.

The enhanced HepMC3 path is now also covered by a two-event mixed fixture. The
first event is charge conserving and receives LowMC photons; the second event
has non-conserved selected signed charge and is written unchanged with no
LowMC photon provenance. Focused CTest coverage exercises both physical
`--enhancement-factor 1` output and enhanced `--enhancement-factor 10` output,
checking per-photon metadata, skipped-event behavior, preservation of the
original topology, and preservation of the event weight vector. The afterburner
and checker also ignore the empty terminal record that the local HepMC3 ASCII
reader can expose after the final listed event.

A more realistic but still inspectable three-event HepMC3 validation fixture
now exists as `realistic_mixed_three_event.hepmc3`. It uses 7 TeV proton beams,
multiple final-state charged particles and neutrals, two charge-conserving
events that receive enhanced LowMC photons, and one non-conserved event that is
skipped unchanged. With the current cheap CTest afterburner settings and
`--enhancement-factor 10`, the first two events produce `5` and `2` LowMC
photons, while the third reports `skipped_nonconserved_charge`. A dedicated
output checker verifies the three-event topology, per-photon and event-level
metadata, neutral counts, and preserved charge diagnostics; the summary
diagnostic verifies that the skipped event has no LowMC metadata.

The same realistic three-event fixture now also has strict
`--failure-policy fail` coverage. In strict mode the afterburner returns
nonzero when it reaches the third event with non-conserved selected charge.
The current diagnostic writer is not atomic: the output path can contain the
two already-written successful events before the failing event aborts the run.
CTest documents this as partial diagnostic output by checking that the summary
diagnostic sees `5` and `2` enhanced LowMC photons for the first two events.
This is not yet a production output-file guarantee.

The HepMC3 atomic-output policy is now documented. The diagnostic afterburner
may keep its current streaming writer behavior, including readable partial
diagnostic output after strict-mode failure. A future stable production
afterburner command or alias should instead use an atomic output contract,
publishing the requested output path only after the full input file succeeds.
The diagnostic now has an explicit `--atomic-output` prototype for that
contract. It writes to `<output>.tmp`, publishes the requested output only after
success, and removes the temporary file on failure. CTest covers successful
atomic output for the charge-conserving fixture and strict atomic failure for
the realistic three-event fixture, verifying in a clean-output case that no
requested output is created and no temporary path remains after failure.

The first focused HepMC3 afterburner failure-mode fixture is now in place. The
neutral-only input event contains no selected charged particles, so the
afterburner reports `skipped_invalid_selection`. CTest coverage verifies that
the default `--failure-policy skip` writes the event back unchanged with no
LowMC photons or metadata, that the output-summary diagnostic reports
`photon_mode=none`, and that `--failure-policy fail` returns a nonzero exit
status for the same input.

The HepMC3 afterburner now also has controlled insufficient-candidate failure
coverage. Using the charge-conserving fixture with
`--enhancement-factor 100`, multiplicity seed `2`, and one shape trial per
requested photon deterministically requests `12` photons but retains only `5`
accepted shape candidates. The event reports
`skipped_insufficient_shape_candidates`, writes no LowMC photons or metadata in
default `skip` mode, round-trips unchanged through the output checker and
summary diagnostic, and returns nonzero with `--failure-policy fail`.

Compact enhancement metadata placement is now settled for the first HepMC3
output path. Per-photon metadata remains the durable fallback. Events that
receive LowMC photons also get duplicate event-level attributes for
`LowMC:photon_mode` and `LowMC:enhancement_factor`; enhanced events also get
event-level `LowMC:weight = 1 / F`. These attributes round-trip through the
local HepMC3 ASCII writer/reader and are covered by the single-event and
multi-event output checks. Skipped events remain unchanged and receive no
LowMC event metadata. Run-level metadata is deferred until a stable file-level
`GenRunInfo` convention is needed.

A read-only HepMC3 LowMC-output summary diagnostic now exists as
`lowmc_hepmc3_output_summary_diagnostic`. It reads HepMC3 ASCII files, reports
one compact line per non-empty event, and does not generate photons or modify
the input file. The summary includes the number of LowMC photons, photon mode,
enhancement factor, summed photon weight, per-photon weight-attribute count,
and whether event-level LowMC metadata is present. CTest coverage checks the
single-event enhanced output and the two-event physical/enhanced outputs,
including the skipped event with zero LowMC photons and no LowMC event
metadata.

A small HepMC3 afterburner count/yield smoke scan now exists as
`scripts/hepmc3_count_yield_smoke.py`. It runs the afterburner on one or more
fixtures for deterministic multiplicity seeds, covers both physical output and
enhanced output, checks successful events for
`requested_photons == photons_added == lowmc_photons`, checks skipped events
for zero added LowMC photons and `photon_mode=none`, and compares the aggregate
requested count with the sum of reported `enhancementFactor * lambda_event`
values using a deliberately loose Poisson band. The CTest run now covers both
`charge_conserving_event.hepmc3` and `realistic_mixed_three_event.hepmc3`,
using four seeds for enhancement factors `1` and `10`; it is a wiring and
normalization smoke test, not a precision production validation.
A second default CTest count/yield scan uses the generated ten-event sample but
limits the run to the first two generated events, two multiplicity seeds, and
enhancement factors `10` and `25`. It uses a tighter two-sigma aggregate
Poisson band and checks that the requested photon counts, afterburner-added
counts, and output-summary counts remain consistent at higher expected yield
without adding a large runtime burden. This still does not replace a dedicated
high-statistics validation job.

An opt-in high-statistics wrapper now exists as
`scripts/hepmc3_count_yield_highstat.py`. It is intentionally not run as part
of the default optional HepMC3 CTest suite, except for a cheap `--help` test.
By default it generates the ten-event HepMC3 sample in its output directory,
runs `scripts/hepmc3_count_yield_smoke.py` over five events, five deterministic
seeds, and enhancement factors `10`, `25`, and `50`, then prints a compact
`highstat` line with cases, events, expected mean, requested photons,
difference, tolerance, and runtime. This is a statistical validation of
Poisson multiplicity and HepMC3 bookkeeping, not a validation of the underlying
event-generator physics.

The first manual high-statistics HepMC3 count/yield validation log has been
recorded with:

```sh
python3 scripts/hepmc3_count_yield_highstat.py \
  --tool-dir build-hepmc3-fresh/adapters/hepmc3 \
  --output-dir build-hepmc3-fresh/adapters/hepmc3/count_yield_highstat_default
```

The final summary was
`highstat input=generated_realistic_sample_highstat.hepmc3 cases=15 events=75 expected_mean=586.509775 requested_photons=608 difference=21.490225 tolerance=72.6538916714 runtime_seconds=165.362`.
This confirms that the aggregate requested photon count is within the
three-sigma Poisson tolerance for the reproducibly generated sample. It is a
manual validation log for multiplicity and HepMC3 bookkeeping, not evidence
about an external generator's physical event distribution.

The manual Pythia-produced HepMC3 validation path is now in place. The
optional helper `scripts/generate_pythia_hepmc3_sample.py` compiles a tiny
temporary Pythia8 + HepMC3 generator using the local `$HOME/pythia8317` and
`$HOME/hepmc3` installations; it does not add a Pythia dependency to LowMC. A
five-event `SoftQCD:nonDiffractive` pp sample at `sqrt(s) = 13 TeV` was
generated with seed `24680`:

```sh
python3 scripts/generate_pythia_hepmc3_sample.py \
  --output build-pythia8-hepmc3/pythia_afterburner_validation/pythia_softqcd_5.hepmc3 \
  --work-dir build-pythia8-hepmc3/pythia_afterburner_validation/generator \
  --events 5 \
  --seed 24680 \
  --ecm 13000
```

The bounded multi-event afterburner validation then processed all five events
with explicit phase space, seeds, enhancement factor `10`, and atomic output:

```sh
lowmc_hepmc3_afterburner_diagnostic \
  --max-events 5 \
  --enhancement-factor 10 \
  --multiplicity-seed 13579 \
  --shape-seed 97531 \
  --kt-min 0.002 --kt-max 0.2 \
  --eta-min -5 --eta-max 5 \
  --phi-min 0 --phi-max 6.283185307179586 \
  --atomic-output \
  pythia_softqcd_5.hepmc3 \
  pythia_softqcd_5_lowmc.hepmc3
```

The run completed successfully after a few minutes, so it remains manual rather
than a default CTest path. It read and wrote five events, added `111` enhanced
LowMC photons, skipped no events, and reported `shape_status=ok` for every
event. The selected outgoing charged multiplicities were `270`, `80`, `40`,
`112`, and `52`, all with signed charge sum `0` and ordinary charge sum `4`.
The requested and retained enhanced photon counts were `58`, `18`, `5`, `24`,
and `6`, with maximum trial/bound ratios between about `0.35` and `0.39`. The
read-only output-summary diagnostic independently reported matching LowMC
photon counts, enhanced mode, enhancement factor `10`, per-photon weight
attributes on every added photon, event-level metadata on every generated
event, and photon weight sums `5.8`, `1.8`, `0.5`, `2.4`, and `0.6`. This
validates the real-generator multi-event HepMC3 record path and metadata
round-trip for a bounded manual sample. It does not make a production
throughput claim.

The first boosted-particle HepMC3 afterburner stress fixture now exists as
`boosted_transverse_p100_event.hepmc3`. It uses 7 TeV incoming proton beams and
two final-state transverse protons with `p = 100 GeV`, preserving signed charge
conservation while exercising a much narrower angular structure than the
earlier small fixtures. The optional HepMC3 CTest path runs the afterburner
with enhancement factor `10`, multiplicity seed `2`, and the standard cheap
shape settings. The event produces `5` enhanced LowMC photons with
`shape_status=ok`, `lambda_event ~= 0.362493`, and
`max_trial_bound_ratio ~= 0.339`; output-file and summary checks verify the
written photon metadata and event-level enhanced metadata.

Controlled HepMC3 shape-envelope failure coverage is now implemented. The
afterburner diagnostic exposes `--shape-safety-factor` and
`--boosted-shape-safety-factor` as diagnostic envelope knobs; their defaults
remain the previous production-like values `3` and `5`. CTest deliberately
sets both to `0.1` on the boosted `p = 100 GeV` fixture, which robustly reports
`skipped_failed_shape_envelope` with `shape_status=envelope_failed` and
`max_trial_bound_ratio > 10`. The test verifies that default `skip` mode writes
the event unchanged with no LowMC photons or metadata, the output-summary
diagnostic reports zero LowMC photons, and `--failure-policy fail` returns
nonzero. This validates failure plumbing; it is not a normal operating point.

The first forward high-rapidity HepMC3 stress fixture now exists as
`forward_pion_pt1_eta5_event.hepmc3`. It uses incoming `pi+ pi-` beams and
final-state `pi+ pi-` with `p_T = 1 GeV`, `eta = +5` and `eta = -5`, preserving
signed charge conservation while matching the earlier LHC-edge discussion. The
optional HepMC3 CTest path runs the afterburner with enhancement factor `10`,
multiplicity seed `2`, and the standard cheap shape settings. The event
produces `2` enhanced LowMC photons with `shape_status=ok`,
`lambda_event ~= 0.0351479`, and `max_trial_bound_ratio ~= 0.339`;
output-file and summary checks verify the written photon metadata and
event-level enhanced metadata.

The explicit extreme `p = 1000 GeV` HepMC3 edge fixture now exists as
`extreme_pion_p1000_event.hepmc3`. It uses incoming `pi+ pi-` beams and
final-state transverse `pi+ pi-` with `p = 1000 GeV`, treating this as an
edge-case stress test rather than a typical operating point. Under the current
cheap afterburner settings the event still produces successful enhanced output:
`7` LowMC photons, `shape_status=ok`, `lambda_event ~= 0.579429`, and
`max_trial_bound_ratio ~= 0.352`. Output-file and summary checks verify the
written photon metadata and event-level enhanced metadata.

The first install/discovery path for optional HepMC3 tools is implemented.
When `LOWMC_ENABLE_HEPMC3=ON`, CMake installs the five non-test HepMC3
executables to `${CMAKE_INSTALL_BINDIR}`:
`lowmc_hepmc3_adapter_diagnostic`,
`lowmc_hepmc3_selected_particle_text`,
`lowmc_hepmc3_weighted_photon_diagnostic`,
`lowmc_hepmc3_afterburner_diagnostic`, and
`lowmc_hepmc3_output_summary_diagnostic`. The installed tools use
`INSTALL_RPATH_USE_LINK_PATH` so the local HepMC3 shared library path is
preserved in the installed binaries when needed. Stable shorter aliases are
deferred; the installed names remain explicitly diagnostic for now.

A minimal end-to-end HepMC3 quickstart helper now exists as
`scripts/hepmc3_quickstart.py`. It reads the checked-in charge-conserving
fixture by default, writes the afterburner output to a selected output
directory, then runs `lowmc_hepmc3_output_summary_diagnostic` on the result.
The script uses build-tree tools by default and also supports installed tools
through `--tool-dir /path/to/bin`. The quickstart defaults to the
afterburner's `--atomic-output` mode, with `--no-atomic-output` available when
the lower-level streaming diagnostic behavior should be exercised explicitly.
A CTest smoke check validates the default atomic enhanced quickstart path and
summary output.

The HepMC3 quickstart now has an installed examples/tutorial layout. When
`LOWMC_ENABLE_HEPMC3=ON`, CMake installs
`share/lowmc/examples/hepmc3/README.md`,
`share/lowmc/examples/hepmc3/hepmc3_quickstart.py`, and the small
`charge_conserving_event.hepmc3` fixture. The installed quickstart is a
separate install-aware wrapper: it derives the installed `bin/` directory from
the install prefix by default and writes output only to a user-selected output
directory. The repository script remains the build-tree quickstart.

The current HepMC3 afterburner readiness status is summarized in
`docs/USAGE.md`. The recommended wording is: first controlled usable HepMC3
prototype, not production afterburner. It is safe for users to try on small
HepMC3 files, installed examples, physical or explicitly enhanced output,
explicit atomic output, and bounded Pythia-produced HepMC3 validation runs,
provided they understand the diagnostic command name and current provisional
boundaries. The implemented checks now cover the controlled single-event
physical output, enhanced output, multi-event skipped-event behavior, the more
realistic three-event fixture, strict failure on that fixture's non-conserved
event, explicit atomic-output success/failure behavior, invalid-selection
skip/fail behavior, insufficient-candidate skip/fail behavior, forced
shape-envelope skip/fail behavior, a generated ten-event pp-like HepMC3 sample
with nine generated events and one deliberate skipped event, multi-fixture and
generated-sample count/yield scans, the boosted transverse
`p = 100 GeV` stress fixture, the forward `p_T = 1 GeV`, `|eta| = 5` stress
fixture, the extreme transverse `p = 1000 GeV` stress fixture, per-photon and
event-level metadata round-tripping through HepMC3 ASCII, the read-only
output-summary diagnostic, installed optional HepMC3 tools, and the repository
quickstart helper. The generated sample is not a substitute for external
generator output; it mainly checks multi-event file handling, metadata,
skipped-event behavior, and summary diagnostics without committing a bulky
sample. A bounded Pythia-produced five-event afterburner validation now
confirms the real-generator multi-event record path, but it took a few minutes
and remains manual validation evidence. `scripts/validate_pythia_afterburner.py`
now wraps the Pythia-to-HepMC3 sample helper, the HepMC3 afterburner
diagnostic, and the read-only output-summary diagnostic. It writes generated
input, afterburner output, afterburner and summary logs, and
`validation_summary.csv` and `validation_events.csv` under a selected output
directory. The summary CSV records aggregate maxima such as selected outgoing
charged multiplicity, requested photons, retained photons, trial/bound ratio,
and runtime. The per-event CSV records the same event-level complexity fields
next to the output-summary comparison. The C++ afterburner diagnostic now
prints `event_seconds` on each event line, so new helper runs can connect
runtime directly to selected charged multiplicity, requested photons, retained
photons, and shape diagnostics. The helper keeps its older estimated timing
column for compatibility with logs produced before `event_seconds` existed.
A two-event helper check with enhancement factor `5` produced `8` enhanced
LowMC photons, skipped no events, reported `shape_status=ok` for both events,
matched the output summary diagnostic, and took about `4.3 s` for sample
generation plus `22.2 s` for the afterburner step.

A first small bounded multi-seed scan with the helper used three two-event
Pythia samples, enhancement factor `5`, and Pythia seeds `24681`, `24682`, and
`24683`. All six events passed: no events were skipped, all shape statuses
were `ok`, output-summary diagnostics matched the afterburner logs, and
`mismatch_count=0` in every helper row. The three runs added `8`, `11`, and
`10` enhanced LowMC photons, with photon weight sums `1.6`, `2.2`, and `2.0`.
The maximum trial/bound ratios were `0.345`, `0.372`, and `0.352`. Generator
times were about `2.2..3.0 s`, while afterburner times varied more strongly:
about `15.2 s`, `42.8 s`, and `30.1 s` for the same two-event count. The
helper's compact text plus CSV summaries are sufficient for manual validation
logs at this stage. The runtime spread is now the main practical observation
from this scan.

A one-event complexity check with Pythia seed `24684`, enhancement factor `5`,
and the same bounded helper path selected `216` outgoing charged particles,
requested and retained `16` LowMC photons, reported `shape_status=ok`, had
`max_trial_bound_ratio ~= 0.371`, matched the output-summary diagnostic, and
measured `event_seconds ~= 65.2` against a total afterburner time of about
`65.2 s`. This confirms that the helper now captures the main event-complexity
proxies needed to interpret slow bounded Pythia-produced events; future scans
should use the new true `event_seconds` field instead of the older evenly
distributed time estimate.

The full five-event default remains manual and can take minutes. Remaining
gaps before production-afterburner claims include broader runtime/complexity
characterization for bounded Pythia afterburner scans using true per-event
timing, naturally occurring shape-envelope failure cases in realistic files,
and a decision on stable production command aliases.

A first true per-event timing scan now combines four one-event Pythia-produced
helper runs with seeds `24681`, `24682`, `24683`, and `24684`, enhancement
factor `5`, multiplicity seed `13583`, and shape seed `97535`. All four
events passed with no skips, all shape statuses were `ok`, all output-summary
diagnostics matched, and the total measured event time was about `90.6 s`.
Selected outgoing charged multiplicities ranged from `34` to `216`; requested
and retained photons summed to `22`, with a maximum of `16` in the slowest
event. The slowest event was the `24684` event with `216` outgoing charged
particles, `16` requested photons, and `event_seconds ~= 65.2`. In this small
sample, event time is strongly correlated with selected outgoing charged
multiplicity and requested/retained photon count; the sample is too small to
quote a general scaling law. The helper-level CSVs were sufficient once
aggregated with the new optional
`scripts/summarize_pythia_afterburner_timing.py` reader; no runner script or
default test-suite timing benchmark is needed at this stage.

The command name decision is settled for now: keep
`lowmc_hepmc3_afterburner_diagnostic` as the only executable name. A shorter
production alias should wait until validation gates and output-mode conventions
are more mature; atomic output is now available as an explicit diagnostic
option.

The accepted direct Pythia in-memory adapter boundary is an optional adapter
target outside `lowmc_core`, translating `Pythia8::Event` records into the
existing generator-independent `AdapterParticleCandidate` vector and then
reusing `selectAdapterParticles(...)`. The default selection matches the
general adapter policy: incoming beams plus final stable charged particles,
momenta in `(E, px, py, pz)` order in GeV, charge from the LowMC PDG helper,
and `inOutSign` assigned only after selection. Pythia particle-data charge may
be useful as a later diagnostic cross-check, but it is not the default charge
source for the first boundary.

The first optional direct Pythia selected-particle diagnostic is now
implemented. `LOWMC_ENABLE_PYTHIA8` defaults to `OFF`; when enabled, CMake
finds Pythia8 through `pythia8-config` and builds `lowmc_pythia8_adapter`,
`lowmc_pythia8_adapter_tests`, and `lowmc_pythia8_adapter_diagnostic` under
`adapters/pythia8`. The adapter maps conventional Pythia event-record
indices `1` and `2` to incoming beam candidates, maps non-beam
`particle.isFinal()` records to outgoing final-stable candidates, converts
momenta to LowMC `(E, px, py, pz)` in GeV, copies PDG ids, and reuses
`selectAdapterParticles(...)` plus the LowMC PDG charge helper. The diagnostic
prints compact selected-particle counts, ignored neutrals, rejected charged
candidates, signed charge sum, and ordinary charge sum. The deterministic toy
mode reports two incoming charged beams, two outgoing charged particles, one
ignored neutral, no rejected charged candidates, and signed charge conservation.
A one-event SoftQCD mode with seed `12345` selects two incoming charged beams,
preserves signed charge conservation, and reports compact diagnostics without
writing HepMC3 output or generating photons.

The direct Pythia selected-particle path is now connected to the existing
weighted photon API through `lowmc_pythia8_weighted_photon_diagnostic`. This is
a separate optional executable, not a mode of the selection diagnostic. It
supports deterministic toy input, a deliberately non-conserved toy input, and
generated SoftQCD events. Charge-conserving events call
`generateWeightedPhotons(...)` and report selected-particle diagnostics,
sample count, weight sum, mean weight, weight RMS, effective sample size,
effective sample fraction, and grid-build `k_T`. Non-conserved selected
currents are reported as `photon_status=skipped_nonconserved_charge`, with no
silent repair. In the current smoke checks, the conserved toy mode with
`8` samples gives `sum_weight ~= 1.47073` and `ESS/N ~= 0.999`; the
non-conserved toy mode is skipped with signed charge sum `-2`; and a one-event
SoftQCD run with seed `12345` and `4` samples selects `46` outgoing charged
particles, has signed charge conservation, and gives `sum_weight ~= 3.59379`.

The direct-Pythia versus Pythia-to-HepMC3 selected-particle cross-check is now
implemented as `lowmc_pythia8_hepmc3_crosscheck`. It is built only when both
`LOWMC_ENABLE_PYTHIA8=ON` and `LOWMC_ENABLE_HEPMC3=ON`. The diagnostic
generates deterministic Pythia8 SoftQCD events in memory, runs the direct
Pythia selection, converts the same events to HepMC3 in memory with
`Pythia8ToHepMC3`, runs the HepMC3 selection, and compares selected incoming
charged count, selected outgoing charged count, ignored neutral candidates,
rejected charged candidates, signed charge sum, and ordinary charge sum. It
also compares weighted photon summaries after the selection diagnostics match.
Both paths call `generateWeightedPhotons(...)` with the same phase space, tuned
grid configuration, sample count, oversampling factor, and photon seed. The
diagnostic intentionally does not compare particle ordering, because
Pythia-to-HepMC3 conversion may reorder event-record entries; the checked
contract is selected-particle diagnostics plus weighted-summary agreement.

The optional CTest smoke check now runs three accepted SoftQCD events with
seed `12345` and four weighted samples per event. The three events select
`46`, `200`, and `90` outgoing charged particles in both paths, with matching
neutral and rejected-charged diagnostics, signed charge sum `0`, and ordinary
charge sum `4`. All weighted summaries match within the diagnostic tolerance;
the aggregate line reports `events_compared=3`, `mismatch_count=0`, and maximum
relative weighted-summary difference about `1.4e-11`. The check remains a
bounded optional diagnostic, not a production dependency.

`docs/ROADMAP.md` now contains a concise HepMC3 production-afterburner
readiness checklist. The checklist separates the current controlled usable
prototype from a future production afterburner, and names the remaining
promotion blockers: broader bounded Pythia-produced afterburner scans, natural
shape-envelope failure evidence, and final command naming.

A first unweighting diagnostic now exists as `lowmc_unweighting_diagnostic`.
It uses the direct two-pion toy input, estimates empirical maximum weights in a
calibration pass, applies a configurable safety factor, and tests
accept/reject unweighting on fresh weighted samples. It supports a global
envelope mode, a sample-calibrated cell-local mode based on the refined
grid-adaptive `eta-phi` cells, and a first `probed-cell` mode that evaluates
exact `w = f/q` on a regular lattice inside every refined cell. The probed
mode now uses a suspect-cell hybrid fallback: if random calibration finds a
larger local value than the deterministic probe lattice, that larger local
value is used for the cell envelope. It reports envelope mode, probe-bin
count, probe rule, accepted photon count, acceptance efficiency,
failed-envelope count, trial/calibration or trial/probed maximum ratios, empty
calibration cells, suspect cells, and fallback uses of the global bound.
`scripts/scan_unweighting_envelope.py` now runs that diagnostic over toy
momenta, seeds, safety factors, and envelope modes and writes CSV plus Markdown
summaries. In the first default global-envelope scan (`p = 1`, `10`,
`100 GeV`; three seeds; safety factors `1.5`, `3`, `10`; `1000 + 1000`
samples per point), all `p = 1` and `10 GeV` envelopes held. For
`p = 100 GeV`, safety factor `1.5` failed in one of three runs with two failed
trial samples, while safety factors `3` and `10` held. Mean observed
acceptance at `p = 100 GeV` was only about `1.8%`, `1.0%`, and `0.27%` for
safety factors `1.5`, `3`, and `10`, respectively.

A focused `p = 100 GeV` per-cell scan shows the expected promise and the
expected danger. With `1000 + 1000` samples, the cell envelope is much more
efficient but unreliable: safety factor `10` still fails in two of three seeds,
and many cells are empty after calibration. With `10000` calibration samples
and `1000` trial samples, safety factors `3` and `10` hold in the three-seed
scan. At safety factor `3`, the mean observed acceptance is about `22%`,
compared with about `0.5%` for the global envelope under the same scan. This
makes local cell envelopes promising enough to pursue, but not yet production
quality because the bound is still empirical and sparse cells still fall back
to the global envelope.

The first explicit per-cell probing pass is more encouraging. In the same
`p = 100 GeV`, three-seed, `1000 + 1000` sample scan, a `3 x 3` probe lattice
inside each refined cell remains too coarse: safety factor `3` still fails,
and the largest local trial/probed ratio is about `8.1`. A `7 x 7` probe
lattice plus the suspect-cell hybrid fallback reduces the largest local ratio
to about `1.94`; with safety factor `3`, all three seeds hold and the mean
observed acceptance is about `21%`. This is far better than the global
envelope at about `1.0%` and also better than the sample-calibrated cell
envelope at about `13.7%` in the same scan. An edge-inclusive `7 x 7` probe
rule reduces the suspect-cell count but gives a slightly lower acceptance and
a larger worst local ratio in this toy. A larger `10000`-trial-sample check at
safety factor `3` holds for both interior and edge-inclusive probing across
three seeds. The result supports explicit per-cell probing as the leading
shape-unweighting candidate, with broader seed, rapidity, and momentum stress
tests still needed before treating it as a production majorant.

The first broader stress scan for the hardened bound now exists as
`scripts/stress_unweighting_bound.py`. It runs the `7 x 7` interior
`probed-cell` mode with the suspect-cell hybrid fallback over five seeds and
three toy cases: transverse `p = 100 GeV`, transverse `p = 1000 GeV`, and
forward `eta = 4`, `p = 100 GeV`. With `1000` calibration samples,
`2000` trial samples, and safety factor `3`, the transverse `p = 100 GeV` and
forward `eta = 4`, `p = 100 GeV` cases have no envelope failures. The forward
case is efficient in this scan, with mean observed acceptance near `30%`. The
first failure mode is the deliberately extreme transverse `p = 1000 GeV` case:
three of five runs fail, with six failed trial samples and a largest
trial/bound ratio of about `7.1`. Repeating that case with safety factor `10`
removes the failures but reduces observed acceptance to about `6%`. This keeps
the hardened probed-cell bound as the leading candidate for realistic
diagnostic cases, while identifying the extreme boosted tail as the next
algorithmic target.

The first boosted-tail follow-up indicates that finer probing is better than a
large safety-factor penalty. For the transverse `p = 1000 GeV` toy, increasing
the interior `probed-cell` lattice from `7 x 7` to `11 x 11` removes the
failures in the same five-seed, `1000` calibration and `2000` trial sample scan
at safety factor `3`. The mean observed acceptance remains about `20%`.
By comparison, keeping `7 x 7` probes and raising the safety factor to `10`
also removes the failures but lowers acceptance to about `6%`. The same
`11 x 11` setting keeps the transverse `p = 100 GeV` and forward `eta = 4`,
`p = 100 GeV` regression cases failure-free, with observed acceptances near
`20%` and `30%`, respectively. `scripts/stress_unweighting_bound.py` now uses
`11 x 11` as its default probe lattice for this stress study.

The first probe-density policy is now documented. The diagnostic default for
`lowmc_unweighting_diagnostic` is `11 x 11` interior probing, matching the
current robust stress-study setting. The scan scripts also report
`probe_evaluations`, so the local-envelope build cost is visible. For the
current tuned two-pion grid with `2257` cells, `7 x 7` probing costs `110593`
exact probe evaluations, while `11 x 11` probing costs `273097`. The
recommended first production-like policy is fixed `11 x 11` probing with
safety factor `3` and the suspect-cell hybrid fallback. Adaptive or
boost-dependent probe density remains a future optimization, not the next API
default.

The first event-level shape-unweighting API is now implemented in the
generator-independent core as `generateShapeUnweightedPhotonCandidates(...)`.
It remains separate from both photon multiplicity and HepMC3 output. The API
returns accepted shape-unweighted photon candidates plus envelope diagnostics.
The configuration includes phase space, trial count, seed, grid config,
calibration-trial count, probe bins, probe rule, safety factor, and
envelope-failure policy. Diagnostics include accepted count, acceptance
efficiency, envelope failures, maximum trial/bound ratio, probe evaluations,
grid cell count, suspect cell count, and fallback count. `Strict` failure
handling clears accepted candidates if an envelope violation is observed;
`Diagnostic` keeps accepted candidates together with failure evidence for
sampler studies. Focused core tests cover ordinary successful generation,
strict failure behavior, diagnostic failure behavior, and invalid
configuration rejection.

A first diagnostic executable for the shape-unweighting API now exists as
`lowmc_shape_unweighting_diagnostic`. It uses the direct two-pion toy input,
calls `generateShapeUnweightedPhotonCandidates(...)`, and prints one
machine-readable key-value summary line. The diagnostic exposes trial count,
calibration-trial count, seed, momentum, particle `eta/phi`, safety factor,
probe bins, probe rule, and failure policy. It reports status, accepted
candidate count, envelope failures, maximum trial/bound ratio, observed and
expected acceptance efficiency, probe evaluations, grid cells, suspect cells,
fallback count, grid-build `k_T`, and signed charge sum. It does not write
HepMC3 output and does not introduce a multiplicity model.

`scripts/stress_unweighting_bound.py` now uses
`lowmc_shape_unweighting_diagnostic` as the preferred stress executable. The
older `lowmc_unweighting_diagnostic` remains available as an internal
envelope-development comparison tool, but future stress scans should run
through the public shape-unweighting API. The current default scan uses five
seeds, `1000` calibration trials, `2000` trial photons, fixed `11 x 11`
interior cell probes, safety factor `3`, and `Diagnostic` failure policy. In
the transverse `p = 100 GeV`, transverse `p = 1000 GeV`, and forward
`eta = 4`, `p = 100 GeV` two-pion toy cases, all `15` runs return `status=ok`
with zero envelope failures and zero fallback uses. The observed acceptance
efficiencies are about `20%`, `20%`, and `30%`, respectively. The largest
reported API `trial/bound` ratio is about `0.97` in the deliberately extreme
transverse `p = 1000 GeV` case, making that case still useful as a boosted-tail
stress test.

The public shape-unweighting API now uses the corrected cell-local envelope
construction. After per-cell bounds `M_c` are built, the sampler draws cells
with probability proportional to the original proposal mass in the cell times
`M_c`, samples uniformly inside the selected cell, and accepts with probability
`w / M_c`. This fixes the normalization problem of sampling from the original
proposal when `M_c` varies by cell. The diagnostic now also reports the
envelope integral. The older `lowmc_unweighting_diagnostic` remains available
only as an internal bound/efficiency development comparison.

With the corrected API sampler, the current `11 x 11` interior, safety-factor
`3`, five-seed stress scan still holds for transverse `p = 100 GeV` and
forward `eta = 4`, `p = 100 GeV`, with mean observed acceptances about `2.9%`
and `29%`, respectively. The deliberately extreme transverse `p = 1000 GeV`
case now fails in all five seeds, with eight failed trial samples and a maximum
trial/bound ratio about `17`. This is a stronger and more production-like
boosted-tail stress result than the older diagnostic; the next work should tune
the corrected envelope for this case and validate the accepted shape.

The first physical multiplicity convention is now chosen: for an event and a
configured `PhotonPhaseSpace`, compute the integrated Low yield
`lambda_event` and draw `N_gamma ~ Poisson(lambda_event)` for ordinary
unweighted event-record photons. This keeps shape unweighting separate from
normalization: photon momenta are drawn from the normalized event-level Low
shape conditional on `N_gamma`. The convention is a first independent-emission
model based on the leading-soft density, not a recoil prescription or a full
multi-photon exponentiation treatment. `oversamplingFactor` remains strictly
diagnostic: an oversampled histogram mode may draw
`N_sample ~ Poisson(oversamplingFactor * lambda_event)` and assign each sample
weight `1 / oversamplingFactor`, preserving the expected weighted yield
`lambda_event`. This mode must not be confused with ordinary event-record
photons.

The remaining normalization work before physical unweighted photons are written
is to estimate `lambda_event` robustly and validate that the corrected
shape-unweighted samples agree with weighted or grid references in focused toy
cases.

`scripts/validate_shape_unweighting.py` now provides the first machine-readable
shape validation for the corrected API sampler. It writes accepted candidates
from `lowmc_shape_unweighting_diagnostic --samples-csv`, compares normalized
eta-phi histograms with the fixed-grid two-particle Low-intensity diagnostic,
and reports L1 shape distance plus eta/phi projection RMS values. In a first
`50000`-trial run with `11 x 11` interior probes and safety factor `3`, the
transverse `p = 100 GeV` case has `1647` accepted candidates, no envelope
failures, L1 `~= 1.23`, eta RMS `~= 0.032`, and phi RMS `~= 0.027`; the forward
`eta = 4`, `p = 100 GeV` case has `14695` accepted candidates, no failures,
L1 `~= 0.64`, and projection RMS near `0.001`. The transverse `p = 1000 GeV`
case remains envelope-failed in this configuration, so its shape metrics are
diagnostic only.

`scripts/tune_shape_unweighting_bound.py` now scans simple fixed probe/safety
settings for the corrected transverse `p = 1000 GeV` sampler. In a first
five-seed scan with `1000` calibration trials and `2000` corrected-shape trials
per seed, `11 x 11` safety factor `3` fails in all seeds. `21 x 21` interior
probing with safety factor `10` still fails in one of five seeds. `21 x 21`
interior probing with safety factor `20` is the first tested failure-free
setting, with maximum trial/bound ratio about `0.72`, but its observed
acceptance is only about `0.14%`. This supports adaptive or targeted boosted
cell probing as the next algorithmic step rather than simply raising the safety
factor globally.

The first adaptive boosted-cell probe prototype is now implemented in the
shape-unweighting API and exposed through
`lowmc_shape_unweighting_diagnostic --probe-policy boosted-cells`. The policy
uses base probes for all cells, promotes cells near charged-particle directions
or cells where random calibration exceeds the base probe maximum, then re-probes
promoted cells with a denser lattice and a larger local safety factor. The
corrected envelope sampler itself is unchanged: cells are still sampled
proportional to `q_c M_c` and accepted with `w / M_c`.

In the expanded transverse `p = 1000 GeV`, five-seed tuning scan,
`adaptive_p11_to_p21_s3_local10` still fails in one seed. The first tested
failure-free adaptive setting is `adaptive_p11_to_p21_s3_local15`: base
`11 x 11` probes, promoted `21 x 21` probes, base safety factor `3`, promoted
local safety factor `15`. It promotes about `41` of `2257` cells on average,
uses about `291000` exact probe evaluations instead of the `995000` required by
full `21 x 21` probing, and gives observed acceptance about `0.17%`. This is
slightly better than the full `21 x 21`, safety-factor-`20` fallback
acceptance of about `0.14%`, with much lower probing cost, but it is still too
inefficient to call production-ready.

The boosted-cell shape-validation run with `50000` trials and local safety
factor `15` returns `status=ok` for the `p = 100`, `p = 1000`, and forward
`eta = 4`, `p = 100 GeV` toy cases. The `p = 1000 GeV` case has only `79`
accepted candidates, so its L1 shape metric is statistics-limited and should be
treated as a sanity check rather than a precise shape validation.

The boosted-cell promotion rule is now tunable through
`--boosted-neighbor-scale`. Setting the scale to `0` promotes only cells that
contain charged-particle directions, plus calibration-suspect cells. In the
transverse `p = 1000 GeV` scan this reduces promoted cells from about `41` to
about `18` and probe evaluations from about `291000` to about `281000`, but it
does not materially improve the envelope integral or acceptance. Lowering the
local safety factor below `15` still fails in the same five-seed scan.

The shape-unweighting diagnostic now also exposes grid-adaptive proposal
parameters. The best new candidate in the first scan is
`adaptive_48_l4_p11_to_p21_s3_local15_n0`: a `48 x 48` initial grid, four
refinement levels, base `11 x 11` probes, promoted `21 x 21` probes, promoted
local safety factor `15`, and zero charged-neighbor halo. In the transverse
`p = 1000 GeV`, five-seed scan it has no envelope failures, mean envelope
integral about `248`, maximum trial/bound ratio about `0.33`, about `20`
promoted cells, and about `806000` probe evaluations. This is a modest
acceptance improvement over the previous boosted-cell prototype and remains
less expensive than the full `21 x 21` all-cell fallback, but the cost increase
is substantial.

The corresponding `50000`-trial shape-validation run returns `status=ok` for
the `p = 100`, `p = 1000`, and forward `eta = 4`, `p = 100 GeV` toy cases. The
`p = 1000 GeV` accepted count rises from `79` to `106`, so this remains
statistics-limited but moves in the right direction.

A targeted charged-direction refinement prototype now splits only the cells
that contain selected charged-particle directions before constructing the
same corrected cell-local envelope. This is a much better cost/performance
tradeoff than globally refining the whole grid. In the transverse
`p = 1000 GeV`, five-seed scan, the best tested setting so far is
`targeted_l4_p11_to_p21_s3_local5_n0`: default `32 x 32` grid with three
ordinary refinement levels, four extra targeted split levels, base `11 x 11`
probes, promoted `21 x 21` probes, zero charged-neighbor halo, and promoted
local safety factor `5`. It has no envelope failures, mean envelope integral
about `8.0`, maximum trial/bound ratio about `0.36`, observed acceptance about
`4.6%`, expected acceptance about `4.8%`, about `17` promoted cells, about
`2353` final grid cells, and about `292000` probe evaluations.

The targeted setting also improves the first shape-validation statistics. In
the `50000`-trial validation run, the `p = 1000 GeV` accepted count rises to
`2536` with no envelope failures, compared with `106` for the global `48 x 48`,
level-4 candidate and `79` for the first boosted-cell prototype. The L1 shape
distance is still a diagnostic statistic rather than a release threshold, but
the accepted statistics are now large enough for the comparison to be more
informative.

`scripts/validate_shape_unweighting_candidates.py` now runs a reproducible
candidate comparison across configurations, toy cases, and seeds. In the first
three-seed, `30000`-trial-per-seed comparison, the targeted candidate has zero
envelope failures for transverse `p = 1 GeV`, transverse `p = 1000 GeV`, and
forward `eta = 4`, `p = 100 GeV`. For the `p = 1000 GeV` stress case it
accepts `4394` candidates, with mean observed acceptance about `4.9%`, mean
envelope integral about `8.0`, about `2353` grid cells, and about `293000`
probe evaluations. The previous `adaptive_p11_to_p21_s3_local15_n0` baseline
has one envelope failure and only `127` accepted candidates in the same stress
comparison, while the global `48 x 48`, level-4 baseline has no failures but
only `164` accepted candidates and about `806000` probe evaluations.

A focused higher-statistics transverse `p = 1000 GeV` check of the targeted
candidate uses five seeds with `100000` trials per seed. It has zero envelope
failures, `24270` accepted candidates, mean observed acceptance about `4.85%`,
mean expected acceptance about `4.83%`, mean envelope integral about `8.0`,
about `17` promoted cells, zero fallback uses, about `2353` grid cells, about
`292000` probe evaluations, and maximum trial/bound ratio about `0.365`.
This is enough to promote the targeted setting to the default diagnostic
shape-unweighting candidate. It is not yet a production majorant guarantee.

The candidate-validation workflow now has named candidate sets. The default
`targeted` set runs only `targeted_l4_p11_to_p21_s3_local5_n0`, while
`--candidate-set comparison` also runs the previous
`adaptive_p11_to_p21_s3_local15_n0` and
`adaptive_48_l4_p11_to_p21_s3_local15_n0` baselines. The lower-level
`lowmc_shape_unweighting_diagnostic` options remain explicit.

Both the fixed-grid truth diagnostic and the shape-unweighting diagnostic now
accept explicit second-pion kinematics through `--second-momentum`,
`--second-particle-eta`, and `--second-particle-phi`. Without those options,
the old back-to-back topology is unchanged. The default candidate validation
now includes an asymmetric two-pion case with a `100 GeV` positive pion at
`eta = 0.35`, `phi = 1.1` and a `35 GeV` negative pion at `eta = -0.75`,
`phi = 4.0`. In the first three-seed, `30000`-trial-per-seed targeted run,
this case has zero envelope failures, `8771` accepted candidates, observed
acceptance about `9.7%`, mean envelope integral about `2.8`, about `2281`
grid cells, about `281000` probe evaluations, and maximum trial/bound ratio
about `0.65`.

A short comparison-preset smoke run confirms that the older baselines remain
selectable. It also hints that the less-symmetric topology may not always be
best served by the targeted split alone: in that very low-stat smoke run, the
global `48 x 48`, level-4 baseline had higher asymmetric-case acceptance but
used almost three times the probe evaluations. This is a useful warning for
future topology tuning, not a change to the current diagnostic default.

The first focused topology-dependence scan compares the old targeted default,
a small targeted split halo, one additional targeted split level, a modest
`40 x 40` global grid, and the global `48 x 48`, level-4 baseline on two
cases: transverse `p = 1000 GeV` and the asymmetric two-pion topology. The
small halo improves the asymmetric case from about `9.7%` to about `20.8%`
observed acceptance without changing the symmetric `p = 1000 GeV` result. The
extra targeted split level improves the symmetric `p = 1000 GeV` case from
about `4.9%` to about `11.4%` observed acceptance without helping the
asymmetric case. Combining both changes in
`targeted_l5_halo0.5_p11_to_p21_s3_local5` improves both cases with zero
envelope failures, about `2377` grid cells, and about `293000-295000` probe
evaluations. The global `48 x 48`, level-4 baseline remains much more
expensive at about `804000-806000` probe evaluations and is still poor for the
`p = 1000 GeV` stress case.

A focused five-seed, `100000`-trial-per-seed check of
`targeted_l5_halo0.5_p11_to_p21_s3_local5` has zero envelope failures in both
tested cases. For transverse `p = 1000 GeV`, it accepts `56916` candidates,
with observed acceptance about `11.4%`, mean envelope integral about `3.4`,
and maximum trial/bound ratio about `0.365`. For the asymmetric two-pion case,
it accepts `104257` candidates, with observed acceptance about `20.9%`, mean
envelope integral about `1.28`, and maximum trial/bound ratio about `0.386`.
This combined setting is now the preferred diagnostic shape-unweighting
candidate for the next validation round.

The default `targeted` validation preset now points to this combined setting.
In the refreshed three-seed, `30000`-trial-per-seed table over all four toy
cases, it has zero envelope failures throughout. It keeps the transverse
`p = 1 GeV` and forward `eta = 4`, `p = 100 GeV` acceptances near their
previous levels, while improving the tuned `p = 1000 GeV` and asymmetric
cases.

The broader high-statistics validation pass now covers all four current toy
cases with five seeds and `100000` trials per seed. The combined candidate has
zero envelope failures and zero fallback uses in every case. The transverse
`p = 1 GeV` case accepts `141745` candidates with observed acceptance about
`28.3%` and maximum trial/bound ratio about `0.346`. The forward
`eta = 4`, `p = 100 GeV` case accepts `146248` candidates with observed
acceptance about `29.2%` and maximum trial/bound ratio about `0.354`. Together
with the previous high-statistics `p = 1000 GeV` and asymmetric checks, the
current toy envelope diagnostics support provisional gates requiring zero
envelope failures, zero fallback uses, and maximum trial/bound ratio below one
with margin. Shape L1 and projection RMS values remain diagnostic only because
some fixed-grid truth references can be resolution limited.

`scripts/validate_shape_unweighting_candidates.py` now exposes these envelope
checks as an optional gate. By default the script remains exploratory. With
`--gate-policy diagnostic`, it adds `gate_policy`,
`gate_max_trial_bound_ratio`, `gate_status`, and `gate_reasons` columns to the
CSV output and a gate column to the Markdown summary. With
`--gate-policy strict`, the same checks return nonzero if any candidate-case
summary fails. A small smoke run verified both a passing diagnostic gate and a
deliberately failing strict gate.

The first full gated toy-suite run used the default `targeted` candidate set,
three seeds, `30000` trials per seed, and
`--gate-max-trial-bound-ratio 0.8`. All `12` candidate-case-seed rows passed.
The aggregate maximum trial/bound ratios were about `0.345` for `mid_p1`,
`0.364` for `mid_p1000`, `0.352` for `forward_eta4_p100`, and `0.379` for
`asymmetric_mid_p100`, with zero envelope failures and zero fallback uses.
Because the worst observed ratio is far below `0.8`, the current provisional
choice is `0.8` rather than a near-one threshold such as `0.95`. Strict mode is
available, but not yet wired into routine CTest runs.

The validation workflow now includes the first event-like multi-particle direct
topology, `event4_inspectable`. It is a fixed four-pion final state with charges
`+ - + -`, momenta `80`, `55`, `25`, and `12 GeV`, and non-collinear directions
`(eta, phi) = (0.20, 0.70), (-0.80, 2.40), (1.60, 4.30), (-2.20, 5.40)`. Both
diagnostics accept it through `--topology event4`. In a three-seed,
`30000`-trial-per-seed gated run with the provisional `0.8` threshold, all five
cases passed. The new `event4_inspectable` case had zero envelope failures,
zero fallback uses, about `22.2%` observed acceptance, and maximum
trial/bound ratio about `0.389`.

A higher-statistics expanded-suite gate check now covers all five current cases
with five seeds and `100000` trials per seed. All `25` candidate-case-seed rows
passed the provisional `0.8` diagnostic gate, with zero envelope failures and
zero fallback uses. The aggregate maximum trial/bound ratios are about `0.346`
for `mid_p1`, `0.365` for `mid_p1000`, `0.354` for `forward_eta4_p100`,
`0.386` for `asymmetric_mid_p100`, and `0.392` for `event4_inspectable`.
The `event4_inspectable` case accepted `110977` candidates, with observed
acceptance about `22.2%` and mean envelope integral about `2.26`.

The suite now also includes `forward4_inspectable`, a fixed four-pion final
state chosen to stress one very forward charged particle. Its charged particles
have charges `+ - + -`, momenta `350`, `120`, `45`, and `30 GeV`, and
directions `(eta, phi) = (4.50, 0.40), (2.80, 1.35), (-0.60, 3.30),
(-2.00, 5.10)`. In the first three-seed, `30000`-trial-per-seed gated run over
all six current cases, all `18` candidate-case-seed rows passed. The
`forward4_inspectable` case had zero envelope failures, zero fallback uses,
accepted `20568` candidates, had observed acceptance about `22.9%`, mean
envelope integral about `1.80`, and maximum trial/bound ratio about `0.372`.

The higher-statistics six-case gate check used five seeds and `100000` trials
per seed. All `30` rows passed the provisional `0.8` diagnostic gate, with zero
envelope failures and zero fallback uses. The aggregate maximum trial/bound
ratios were about `0.346` for `mid_p1`, `0.365` for `mid_p1000`, `0.354` for
`forward_eta4_p100`, `0.386` for `asymmetric_mid_p100`, `0.392` for
`event4_inspectable`, and `0.396` for `forward4_inspectable`. The
`forward4_inspectable` case accepted `114037` candidates, with observed
acceptance about `22.8%` and mean envelope integral about `1.80`.

The suite now includes `multi8_inspectable`, the first higher-multiplicity
fixed topology. It has eight outgoing charged pions with charges
`+ - + - + - + -`, momenta `140`, `95`, `70`, `50`, `35`, `25`, `18`, and
`12 GeV`, and spread-out directions in eta and phi. In the first three-seed,
`30000`-trial-per-seed gated run over all seven current cases, all `21` rows
passed. The `multi8_inspectable` case had zero envelope failures, zero fallback
uses, accepted `21595` candidates, had observed acceptance about `24.0%`, mean
envelope integral about `5.89`, and maximum trial/bound ratio about `0.382`.
The suite-wide maximum in this run remained the `event4_inspectable` value of
about `0.389`.

The higher-statistics seven-case gate check used five seeds and `100000` trials
per seed. All `35` rows passed the provisional `0.8` diagnostic gate, with zero
envelope failures and zero fallback uses. The aggregate maximum trial/bound
ratios were about `0.346` for `mid_p1`, `0.365` for `mid_p1000`, `0.354` for
`forward_eta4_p100`, `0.386` for `asymmetric_mid_p100`, `0.392` for
`event4_inspectable`, `0.396` for `forward4_inspectable`, and `0.384` for
`multi8_inspectable`. The `multi8_inspectable` case accepted `120378`
candidates, with observed acceptance about `24.1%` and mean envelope integral
about `5.89`.

The diagnostics now also support a minimal direct-particle text input format
with columns `charge inOutSign mass p eta phi`. The first checked-in fixture,
`fixtures/event4_particles.lowmc`, duplicates the built-in `event4`
topology and is used as `event4_text_inspectable` in the candidate-validation
workflow. In the first three-seed, `30000`-trial-per-seed gated run over all
eight current cases, all `24` rows passed the provisional `0.8` diagnostic
gate, with zero envelope failures and zero fallback uses. The text fixture
matched the built-in `event4_inspectable` aggregate exactly: `19991` accepted
candidates, observed acceptance about `22.2%`, mean envelope integral about
`2.26`, and maximum trial/bound ratio about `0.389`. The suite-wide maximum in
this run was also about `0.389`.

The first clustered text fixture, `fixtures/clustered8_particles.lowmc`,
is now included as `clustered8_text_inspectable`. It contains four nearby
outgoing charged pions plus four wider-angle outgoing pions, with total charge
zero. In the first three-seed, `30000`-trial-per-seed gated run over all nine
current cases, all `27` rows passed the provisional `0.8` diagnostic gate,
with zero envelope failures and zero fallback uses. The clustered case accepted
`18779` candidates, had observed acceptance about `20.9%`, mean envelope
integral about `3.87`, and maximum trial/bound ratio about `0.404`. The
suite-wide maximum in this run was the clustered case, so it is the current
most demanding inspectable fixture by this envelope diagnostic.

The first incoming/outgoing sign-convention fixture,
`fixtures/incoming_outgoing4_particles.lowmc`, is now included as
`incoming_outgoing4_text_inspectable`. It contains two incoming `+1` charged
particles and two outgoing `+1` charged particles. The ordinary charge sum is
therefore nonzero, but the signed Low-current sum is zero when incoming
particles carry `inOutSign = incoming`. In the first three-seed,
`30000`-trial-per-seed gated run over all ten current cases, all `30` rows
passed the provisional `0.8` diagnostic gate, with zero envelope failures and
zero fallback uses. The incoming/outgoing case accepted `20976` candidates, had
observed acceptance about `23.3%`, mean envelope integral about `2.11`, and
maximum trial/bound ratio about `0.379`. The suite-wide maximum remained the
clustered case at about `0.404`.

The first generator-derived sampler-validation fixture is now in place. The
optional HepMC3 diagnostic
`lowmc_hepmc3_selected_particle_text` reads
`adapters/hepmc3/tests/fixtures/validation_event.hepmc3`, runs the HepMC3
adapter selection, and writes the selected charged particles to
`fixtures/hepmc3_validation_event_particles.lowmc`. The HepMC3 fixture is
in MEV and contains two slightly off-axis incoming charged beams, two outgoing
stable charged protons, and one neutral photon. Adapter diagnostics report two
selected incoming charged particles, two selected outgoing charged particles,
one ignored neutral candidate, zero rejected charged candidates, signed charge
sum zero, and ordinary charge sum four. The optional HepMC3 CTest path
regenerates the text fixture and compares it with the checked-in copy.

The derived sampler case is included as `hepmc3_validation_text_inspectable`.
In the first three-seed, `30000`-trial-per-seed gated run over all eleven
current cases, all `33` rows passed the provisional `0.8` diagnostic gate,
with zero envelope failures and zero fallback uses. The HepMC3-derived case
accepted `26265` candidates, had observed acceptance about `29.2%`, mean
envelope integral about `0.119`, and maximum trial/bound ratio about `0.355`.
The suite-wide maximum remained the clustered case at about `0.404`.

The validation-suite review concludes that the current eleven-case suite is
broad enough for the next first-prototype step. It covers boosted two-particle
stress cases, forward kinematics, asymmetric two-particle kinematics, fixed
four-particle and eight-particle event-like topologies, a clustered
charged-particle fixture, explicit incoming/outgoing sign convention, and a
HepMC3-derived adapter-selection input. No obvious additional toy topology is
needed before moving to event-level yield and multiplicity diagnostics.

The current first-prototype diagnostic gate is provisional: `status=ok`, zero
envelope failures, zero fallback uses, and `max_trial_bound_ratio < 0.8` for
the explicit validation run. Shape L1 and eta/phi projection RMS values should
continue to be reported but remain informative only, because fixed-grid
references can be resolution limited. This gate supports implementing
`lambda_event` estimation and Poisson multiplicity diagnostics next; it does
not yet justify writing ordinary unweighted photons into HepMC3 events.

The first event-level yield prototype is implemented in the core as
`estimateEventYield(...)`. It computes `lambda_event` by deterministic midpoint
integration in `(log k_T, eta, phi)`, including the physical-measure Jacobian
`k_T^2 dlog(k_T) deta dphi`, and reports method name, phase-space grid,
evaluation count, fine-grid estimate, half-resolution coarse estimate, and the
absolute/relative coarse-fine difference. This is a transparent baseline
integrator, not yet the final production uncertainty model.

`samplePhotonMultiplicity(...)` adds the first count-only multiplicity
diagnostic. It samples `N_gamma ~ Poisson(lambda_event)` for
`oversamplingFactor = 1`; for diagnostic oversampling it samples
`N_sample ~ Poisson(oversamplingFactor * lambda_event)` and reports per-sample
analysis weight `1 / oversamplingFactor`. The new
`lowmc_event_yield_diagnostic` writes one CSV summary row and does not write
event-record photons. Focused tests cover deterministic yield estimation,
log-`k_T` scaling, invalid inputs, and the oversampled Poisson convention.

`scripts/validate_event_yield.py` now scans the deterministic yield estimator
against the existing full weighted sampler for two-pion toy cases
`p = 1`, `10`, `100`, and `1000 GeV`, plus one text-fixture event. It writes a
summary CSV and Markdown file and clearly separates the deterministic grid
sequence from the weighted-sampler statistical uncertainty. The first default
run with `20000` weighted samples per case shows that the midpoint estimator is
reasonable only for the inspectable `p = 1 GeV` case at the current fine grid:
`lambda_event = 0.09283` versus weighted mean `0.09255`, about `1.7` weighted
standard errors. For boosted cases the deterministic grid is under-resolved:
at `p = 10 GeV`, the fine-grid value is `0.168` versus weighted mean `0.191`;
at `p = 100 GeV`, `0.173` versus `0.289`; and at `p = 1000 GeV`, `0.173`
versus `0.386`, with large grid-to-grid movement. The text fixture also shows
large grid sensitivity. Therefore the current midpoint yield estimator is a
baseline diagnostic and must not be used yet as the production event
multiplicity integral for boosted or multi-particle events.

The first adaptive event-yield prototype is implemented as
`estimateAdaptiveEventYield(...)` and exposed through
`lowmc_event_yield_diagnostic --method adaptive`. It builds a grid-adaptive
eta-phi cell partition at each log-`k_T` midpoint and then integrates
deterministically over the final cells with configurable midpoint sub-bins.
The validation script now reports both midpoint and adaptive estimates against
the weighted-sampler reference. In the default `20000` weighted-sample run,
the adaptive L7 configuration (`32 x 32` initial cells, seven refinement
levels, refine fraction `0.1`, `2 x 2` sub-bins) gives:
`p = 1 GeV`, `lambda_event = 0.09269` versus weighted mean `0.09255`;
`p = 10 GeV`, `0.19095` versus `0.19069`; `p = 100 GeV`, `0.28951` versus
`0.28939`; and the text fixture, `0.50031` versus `0.49960`.

For the extreme `p = 1000 GeV` stress case, the adaptive estimator improves
substantially over the midpoint grid: the L7 value is `0.36562` versus weighted
mean `0.38632`, compared with midpoint-fine `0.17254`. However, L5-to-L7
movement is still about `15%`, and the weighted reference has low
`ESS/N` (`0.00487`) in this run. This supports the adaptive-cell strategy but
does not yet define a final production validation gate for very high boosts.

`scripts/tune_event_yield.py` now scans adaptive refinement levels,
refinement fractions, and subcell quadrature settings. In the first tuning run
with `50000` weighted samples per case, adaptive L7 differs from the L9 tuning
reference by `0.014%` for `p = 1 GeV`, `0.014%` for `p = 10 GeV`, `0.026%`
for `p = 100 GeV`, and `0.028%` for the text fixture. The corresponding
weighted-reference ESS/N values are `0.937`, `0.903`, `0.161`, and `0.624`.
These pass the provisional ordinary-case gate: adaptive L7 or better, less
than `1%` deterministic difference from the L9 tuning reference, and
weighted-reference `ESS/N >= 0.05`.

The `p = 1000 GeV` stress case remains separate. In the same tuning run,
adaptive L7 is `6.3%` below L9, L8 is `1.3%` below L9, and L9 gives
`lambda_event = 0.39010`. The weighted reference is `0.40122` but has
`ESS/N = 0.00373`, so it is useful as a sanity check but not as the primary
gate. The current stress-case policy is to require L8/L9-level deterministic
checks and to avoid treating the weighted comparison alone as decisive.

`lowmc_count_shape_diagnostic` now composes adaptive event yield, Poisson
count sampling, and the corrected shape-unweighting API. It writes a summary
CSV row and can optionally write retained candidate photons to a second CSV.
For ordinary cases it uses the adaptive L7 yield settings; for stress cases it
uses L9-level yield refinement. It keeps the midpoint yield estimator out of
the count path and does not write HepMC3 event records. A representative
oversampled `p = 10 GeV` two-pion run with `oversamplingFactor = 100` produced
`lambda_event = 0.19095`, Poisson mean `19.09`, `11` requested diagnostic
photons, `11` retained candidates, analysis weight `0.01`, zero envelope
failures, and `max_trial_bound_ratio = 0.344`.

`scripts/validate_count_shape.py` now runs the count-plus-shape diagnostic over
controlled fixtures. The first three-seed validation uses
`oversamplingFactor = 100` for two-pion `p = 1`, `10`, `100`, and
`1000 GeV`, `event4_particles.lowmc`, and
`hepmc3_validation_event_particles.lowmc`. All `18` rows passed:
`status=ok`, `envelope_ok=1`, `enough_candidates=1`, and zero envelope
failures. Aggregated requested/retained diagnostic counts were `25/25`,
`44/44`, `71/71`, `99/99`, `130/130`, and `10/10`, respectively. The largest
observed `max_trial_bound_ratio` was `0.390` in the event4 text fixture. The
HepMC3-derived selected-particle text fixture is now included in this
diagnostic validation, but no HepMC3 event output is written.

`docs/STATUS.md` is a living file. Update it when the project state, active
milestone, or known open questions change.

## In Scope Now

- C++17-or-later core
- generator-independent event representation
- charged external particles with physical four-momenta
- coherent Low current `J^mu(k)`
- photon intensity `S(k) = -J^2`
- photon phase-space bounds for `k_T`, `eta`, and `phi`
- unit tests for the basic physics and scaling properties
- CSV diagnostics for simple controlled toy events
- optional Python plotting for diagnostic CSV output

## Out of Scope For The Current Skeleton

- production HepMC3 event writing or full HepMC3 workflow integration
- Pythia8 interface
- recoil prescription
- multiple-photon exponentiation
- virtual photons and dielectrons
- integration with existing QED radiation tools
- ROOT-dependent plotting or analysis
- unweighted event generation

## Immediate Next Milestone

Follow the single task in `docs/NEXT_TASK.md`.

## Open Questions

- Is a no-incoming-particle LowMC afterburner vertex the best HepMC3 topology
  for event-wide soft photons, or should a different diagnostic topology be
  used?
- Which exact HepMC3 attribute naming style is most convenient for downstream
  tools?
- Can adaptive/per-cell unweighting bounds improve boosted-case acceptance
  beyond the first boosted-cell prototype without introducing envelope failures?
- How stable are empirical unweighting envelopes for forward charged-particle
  topologies?
- How should the adaptive `eta-phi` proposal be tuned beyond the first local-box
  diagnostic baseline?
- Which proposal diagnostics should be release-blocking?
