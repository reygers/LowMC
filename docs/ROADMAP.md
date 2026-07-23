# ROADMAP

## Publication And Method Framing

LowMC is intended as a reusable software tool in the spirit of a Computer
Physics Communications methods paper.

The goal is not to invent a new generic Monte Carlo sampling philosophy. The
goal is to apply established adaptive Monte Carlo and event-generator methods
to the specific hard sampling problem posed by coherent leading-Low
soft-photon emission:

- formulate the event-level Low photon density cleanly in
  `k_T, eta, phi`
- expose the sharp, event-dependent structures near charged-particle
  directions as a central algorithmic problem
- use weighted adaptive importance sampling before relying on unweighting
  envelopes
- validate proposal quality, weight variance, envelope stability, and
  performance with explicit diagnostics
- keep the physics core independent of event-record adapters
- provide HepMC3/Pythia workflows once the weighted and unweighted modes are
  understood

The intended user workflows are summarized in `docs/USAGE.md`.

### Position Relative To VinciaQED

VinciaQED is a close existing comparison point, not only an adjacent decay
tool. Its multipole mode provides a fully coherent QED shower whose scalar
eikonal limit overlaps directly with LowMC's leading-Low current. LowMC should
therefore not be framed merely as a new generic coherent QED shower.

LowMC's intended niche is a generator-independent leading-soft reference and
afterburner for completed event records. It should provide explicit control of
the charged particles included in the coherent current, detector-oriented
soft-photon phase-space limits, transparent weighted densities and diagnostics,
and workflows that do not require access to an event generator's internal
shower history. VinciaQED instead supplies an ordered shower integrated with
Pythia's parton systems, recoil maps, collinear kernels, repeated emissions,
and other QED evolution.

The different definitions of the coherent system are a physics question that
must be tested. In particular, LowMC's default use of incoming beams plus final
stable charged particles is not automatically equivalent to Vincia's
system-, scale-, and decay-organized coherence. Establishing where those
descriptions agree or differ is part of the scientific validation of LowMC,
not just an implementation cross-check.

## User-Facing Examples

LowMC should eventually provide a small set of example programs in the style of
Pythia `mainXX.cc` examples. These should demonstrate the intended workflows
rather than only test internals. In particular, one early example should show
enhanced Low photon production with an explicit enhancement factor, because
physical Low photon counts can be too sparse for practical histograms in a
limited event sample. The example should make clear that enhanced HepMC output
samples

```text
N_sample ~ Poisson(enhancementFactor * lambda_event)
```

and that photon histograms must be filled with `LowMC:weight = 1 /
enhancementFactor` to recover the physical Low yield.

Candidate examples:

- minimal two-charged-particle direct-core example
- HepMC3 afterburner with physical `enhancementFactor = 1`
- HepMC3 afterburner with enhanced photons for smooth histograms
- direct Pythia8 in-memory workflow filling histograms without HepMC output
- optional Pythia8-to-HepMC3 round-trip workflow

VinciaQED and PHOTOS are important related tools and comparison points.
VinciaQED overlaps directly in the coherent soft kernel but is an integrated
ordered Pythia shower, while PHOTOS primarily targets QED radiation in decays.
LowMC's different target is a generator-independent leading-soft calculation
over a configurable set of charged external particles.

## Milestone 1: Core Physics Skeleton

Goal: implement the minimal generator-independent Low-current calculation.

- define lightweight core data structures
- implement four-vector arithmetic with metric `(+,-,-,-)`
- compute the coherent current `J^mu(k)` in `O(N)`
- compute `S(k) = -J^2`
- add focused unit tests for charge conservation, finite-mass behavior, and
  current-evaluation scaling

## Milestone 2: Weighted Photon Generation

Goal: sample soft photons from the leading Low spectrum in user-defined
phase-space limits.

- define photon phase-space configuration
- evaluate the density
  `dN_gamma / (k_T dk_T d eta d phi) = alpha / (2 pi)^2 S(k)`
- return weighted photons
- document the meaning of weights and normalization
- keep phase-space acceptance separate from internal proposal strategy
- treat efficient sampling near charged-particle directions as a design topic,
  not as an incidental implementation detail
- add sampler diagnostics such as weight spread, acceptance rate, and refined
  cell counts

### Physics-Informed Analytic Proposals

PHOTOS and VinciaQED both generate trial photons from analytically invertible
soft/collinear approximations and apply the complete target kernel only as a
final rejection or weight correction. LowMC already follows this principle in
photon `k_T` by sampling logarithmically. A diagnostic workstream should test
whether extending it to photon angles can reduce grid construction and
shape-unweighting cost without changing the leading-Low density.

The first diagnostic prototype of a PHOTOS-inspired analytic leg-centered
proposal now exists. It:

- uses the finite-mass denominator
  `p_i . k = E_i E_gamma (1 - beta_i cos(theta_i_gamma))` to define a
  normalized angular component around each selected charged particle
- mixes the charged-particle components with an explicit uniform floor so the
  proposal has support over the complete configured acceptance
- continues to evaluate the exact coherent `S(k) = -J(k)^2` for every returned
  weight or accept/reject decision, so interference is never approximated in
  the physical density
- begins as a weighted, diagnostic-only option and has unit coverage plus a
  first three-seed comparison on four controlled multi-particle fixtures
- can be compared against the uniform, local-box, smooth-triangular, and grid-adaptive
  proposals using effective sample fraction, weight RMS, unweighting
  acceptance, maximum target/proposal ratio, proposal-build time, total
  runtime, and Low-current evaluation count
- retains the current validated grid path as the default unless the analytic
  proposal passes both efficiency and correctness checks on boosted toys and
  realistic Pythia-produced events

The initial controlled comparison finds large effective-sample improvements
for boosted and clustered fixtures but a clear loss relative to the grid for a
low-boost HepMC3 validation fixture. The exact contract and first results are
recorded in `docs/LEG_CENTERED_PROPOSAL_PLAN.md`. The next gate is an
independent Pythia-event tail study. That first gate is now complete on sixty
`e+e-` hadronic and sixty SoftQCD events. The analytic proposal has lower
per-sample ESS than the grid in most realistic events, but substantially higher
effective-sample throughput and a far more stable observed residual maximum.
This justifies a separate diagnostic analytic-envelope prototype with strict
failure detection; it does not justify changing the production default. That
prototype is now complete for the same `120` Pythia events. A bound calibrated
from `20000` proposal draws and tested on `100000` independent draws per event
has zero overshoots at all three predeclared safety factors. At factor `1.2`,
the largest trial/bound ratio is `0.844` for hadronic Z events and `0.915` for
SoftQCD, with mean expected acceptance `7.97%` and `6.58%`. This clears a
finite-sample diagnostic gate only: the bound is empirical, failure detection
must remain mandatory, and the current cell-local grid stays the production
default. Further work should be driven by a concrete need for an independent
reference or by broader stress coverage, rather than making this PHOTOS-
inspired branch the main LowMC workstream.

A second, more exploratory prototype may use Vincia-inspired charged-pair
antennae. With signed effective charges `c_i = eta_i Q_i`, candidate pairs can
be proposed using positive strengths such as `max(0, -c_i c_j)` and angular
shapes based on

```text
(p_i . p_j) / ((p_i . k) (p_j . k)).
```

The exact Low current must again remain the target. This prototype should only
follow the leg-centered study because a naive pair mixture has `O(N^2)`
components and may reproduce Vincia's known conservative-overestimate cost at
high charged multiplicity. Its proposal-density evaluation and memory scaling
must be included in the comparison.

Sudakov ordering, hardest-emission competition, recoil maps, virtual/no-photon
normalization, spin-dependent collinear kernels, and iterated shower evolution
solve different physics problems and should not enter this sampling workstream.
They require separate explicit physics-model decisions if considered later.

## Milestone 3: Event Interfaces

Goal: connect the core to external event records without introducing generator
dependencies into the core.

- add a HepMC3 adapter
- add a thin Pythia8 interface outside the core
- define configurable charged-particle selection policies
- keep adapter behavior separate from the physics kernel

### HepMC3 Production-Afterburner Readiness

The current `lowmc_hepmc3_afterburner_diagnostic` is a first controlled usable
prototype. It is appropriate for small-file user trials and validation studies,
but not yet a production afterburner. Before promoting its behavior to a
production-style command or alias, the minimum readiness checklist is:

- satisfied: optional HepMC3 build keeps `lowmc_core` generator-independent
- satisfied: default selection includes incoming beams and final stable charged
  particles, with diagnostics for charge conservation and skipped events
- satisfied: output topology uses attributed stable root photons and preserves
  original vertices, momenta, event weights, and no-recoil behavior
- satisfied: physical and enhanced output modes are explicitly labeled, with
  per-photon metadata and event-level duplicate metadata for generated events
- satisfied: skip/fail behavior is covered for invalid selection, non-conserved
  selected charge, insufficient candidates, and forced envelope failure
- satisfied: read-only output-summary diagnostics and repository quickstart
  checks exist for small fixtures
- satisfied: explicit `--atomic-output` mode publishes the requested output
  only after success and leaves no temporary output after strict failure in the
  checked fixture
- satisfied: validation includes a reproducibly generated ten-event HepMC3
  sample with varied pp-like final states, generated LowMC photons, and a
  deliberate skipped event
- partial: a Pythia-produced HepMC3 sample can be generated manually; a
  bounded five-event afterburner run has passed with matching output-summary
  diagnostics, and an opt-in helper now records logs, a CSV summary, and
  runtime; a first three-seed, two-event-per-seed scan passed with no skips,
  no summary mismatches, all shape statuses `ok`, and visible afterburner
  runtime spread; the helper now also writes per-event complexity rows and
  aggregate maxima for selected charged multiplicity, requested/retained
  photons, shape diagnostics, and true per-event afterburner timing; a first
  four-event timing scan passed and showed the expected association between
  runtime, selected charged multiplicity, and requested photons; a first
  three-event stage-timing scan shows shape unweighting accounts for about
  `90%` of measured event time
- open: naturally occurring shape-envelope failures should be observed or
  bounded in realistic event files, not only forced by unsafe diagnostic knobs
- satisfied: default count/yield checks include both a loose multi-fixture
  smoke scan and a tighter generated-sample aggregate scan
- satisfied: an opt-in high-statistics count/yield validation script exists and
  is kept out of the default CTest runtime except for `--help`
- satisfied: a manual high-statistics count/yield validation log has been
  recorded for the reproducibly generated HepMC3 sample
- deferred: broader boosted HepMC3 stress coverage is lower priority because
  boosted sampling is already covered in the core and single-event HepMC3
  fixtures
- satisfied: installed examples use a self-contained
  `share/lowmc/examples/hepmc3/` layout with README, quickstart script,
  fixture, and install-tree smoke coverage
- satisfied: first optional direct Pythia in-memory selected-particle
  diagnostic translates `Pythia8::Event` into `AdapterParticleCandidate`
  records and reuses the dependency-free selection helper
- satisfied: first optional direct Pythia in-memory weighted-photon diagnostic
  calls the generator-independent weighted API and reports weighted summaries
  without HepMC3 output writing
- satisfied: direct Pythia versus Pythia-to-HepMC3 selected-particle and
  weighted-summary cross-check matches for a bounded three-event SoftQCD
  sample
- open: final command naming should be settled before a production alias such
  as `lowmc_hepmc3_afterburner` is added

The next speed milestone should prototype a validated spatial index for
shape-unweighting cell lookup/sampling, because stage timing now shows shape
unweighting dominates bounded Pythia-produced afterburner runtime.

## Milestone 4: Validation

Goal: compare LowMC output against analytic and event-level expectations.

- validate simple charge-neutral toy systems
- add simple diagnostics that write CSV or JSON for plotting
- provide optional Python/matplotlib plotting scripts for diagnostics
- check finite-mass collinear regularization
- validate sampler behavior for highly boosted charged particles
- test coherent interference patterns
- compare selected event samples against independent Low-theorem calculations
- add a controlled VinciaQED multipole comparison using the same
  charge-conserving particle configuration in both programs
- first compare the one-photon, very-soft limit with Vincia spin-dependent
  collinear terms disabled, where the scalar eikonal shapes should agree up to
  documented normalization and phase-space conventions
- then isolate differences from charged-system selection, spin-collinear
  terms, recoil, ordering, and repeated emissions
- determine whether the global current built from incoming beams and final
  stable charged particles is an adequate approximation for the intended
  experimental soft-photon observables

### DELPHI Hadronic Bremsstrahlung Benchmark

The next central physics-validation target is the efficiency-corrected
inner-bremsstrahlung spectrum in DELPHI hadronic `Z` decays. The publication
provides ten `8 MeV` bins for photon transverse momentum relative to the
closest accepted jet and an integrated prediction of
`17.10e-3 photons / accepted jet` for `0.2 < E_gamma < 1 GeV` and
`pT_jet < 80 MeV`.

The first benchmark must compare the paper's practical per-jet
forward-hemisphere construction with a complete charge-conserving current made
from incoming beams plus final charged hadrons, using LUCLUS-like jets with
`djoin = 3 GeV`. LowMC's general beams-plus-all-stable-charged selection should
be reported as a further diagnostic alternative.
The existing exploratory analysis source supplies useful Pythia and jet
definitions, but its active primary-quark-only configuration and arbitrary
plot rescaling are not normalization references.

The implementation and predeclared gates are recorded in
`docs/DELPHI_HADRONIC_BENCHMARK_PLAN.md`. This work should first prove that
LowMC reproduces the calculated bremsstrahlung baseline. Comparison with the
measured excess comes only afterward and must not be described as an
explanation of that excess.

### PHOTOS Soft-Decay Benchmark

Add a direct `Z -> mu+ mu-` comparison with PHOTOS as the first external
afterburner benchmark. For a neutral parent and the same two outgoing muons,
PHOTOS and LowMC contain the same charge-conserving leading-soft current, so
their real-photon predictions should converge in the strict soft limit. This
is a clean test of LowMC's charge signs, normalization, finite-mass behavior,
interference pattern, sampling, and unweighting in a domain where PHOTOS is a
mature independently validated reference.

The detailed benchmark contract, staging, and success criteria are recorded in
`docs/PHOTOS_BENCHMARK_PLAN.md`.

The first comparison should:

- use identical fixed Born events with a `Z` at rest and back-to-back physical
  muon four-momenta
- configure PHOTOS for interference-enabled single-photon generation with its
  universal kernel before comparing process-specific hard-photon corrections
- align the fine-structure constant, muon mass, infrared cutoff, photon
  acceptance, and normalization convention explicitly
- compare the accepted real-photon yield, photon-energy spectrum, angle to the
  nearest muon, and wide-angle interference-sensitive distributions
- repeat the comparison for progressively smaller maximum photon energy and
  demonstrate that the LowMC/PHOTOS ratio approaches one as
  `E_gamma / m_Z -> 0`
- define angles relative to the common pre-radiation Born muons, or otherwise
  quantify PHOTOS recoil effects, because LowMC currently preserves the input
  charged-particle momenta
- build on the existing `lowmc_delphi_muon_fsr_diagnostic`, while extending it
  from a fixed-topology shape/data comparison to a reproducible rate-and-shape
  comparison with PHOTOS

Agreement would validate LowMC where the two physical models overlap. It would
not imply that PHOTOS's decay-by-decay event traversal reproduces LowMC's
single configurable current over a general scattering event with distinct
incoming legs. PHOTOS remains the natural tool for ordinary decay radiation;
LowMC's broader target remains a generator-independent leading-soft
calculation for a user-selected scattering configuration.

## Milestone 5: Extensions

Potential future work:

- unweighted generation
- multiple-photon exponentiation
- recoil handling
- virtual photons
- dielectron generation
- double-counting treatment with VinciaQED, other QED showers, PHOTOS, and
  radiative decays

## Open Questions

- Which validation observables should be considered release-blocking?
- Should exponentiation be part of LowMC's core model or a separate layer?
- How much event-record policy should live in adapters versus user configuration?
- Which charged particles form a physically coherent system for each intended
  experimental use case, and when are lifetime or formation-length effects
  needed?
- Which existing-radiation metadata or configuration is required to prevent
  LowMC from being applied to overlapping QED photon sources?
