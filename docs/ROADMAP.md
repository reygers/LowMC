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

PHOTOS is an important related tool and comparison point, but LowMC has a
different target: event-wide coherent radiation from a configurable set of
charged external particles, rather than QED radiation in decays.

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
- partial: a Pythia-produced HepMC3 sample can be generated manually and the
  first event has passed the afterburner; broader multi-event Pythia runtime
  validation remains open
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
- satisfied: first direct Pythia versus Pythia-to-HepMC3 selected-particle and
  weighted-summary cross-check matches for a deterministic one-event SoftQCD
  sample
- open: final command naming should be settled before a production alias such
  as `lowmc_hepmc3_afterburner` is added

The next adapter milestone should broaden the direct-vs-HepMC3 Pythia
cross-check to a small multi-event or multi-seed sample while keeping it a
bounded diagnostic, not a production dependency.

## Milestone 4: Validation

Goal: compare LowMC output against analytic and event-level expectations.

- validate simple charge-neutral toy systems
- add simple diagnostics that write CSV or JSON for plotting
- provide optional Python/matplotlib plotting scripts for diagnostics
- check finite-mass collinear regularization
- validate sampler behavior for highly boosted charged particles
- test coherent interference patterns
- compare selected event samples against independent Low-theorem calculations
- define comparison points with PHOTOS where decay radiation overlaps
  technically, while keeping LowMC's event-wide current distinct

## Milestone 5: Extensions

Potential future work:

- unweighted generation
- multiple-photon exponentiation
- recoil handling
- virtual photons
- dielectron generation
- double-counting treatment with showers or PHOTOS

## Open Questions

- Which validation observables should be considered release-blocking?
- Should exponentiation be part of LowMC's core model or a separate layer?
- How much event-record policy should live in adapters versus user configuration?
