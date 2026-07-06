# ADAPTER CONTRACT

This document defines the boundary between generator-specific event records and
the generator-independent LowMC core.

## Purpose

Adapters translate external event records, such as HepMC3 or Pythia8, into the
minimal charged-particle input required by the LowMC core.

The core should not know generator-specific status codes, barcode conventions,
mother/daughter record structure, or particle-selection policy. It should only
receive charged external particles that are already classified as incoming or
outgoing.

## Core Input

For each selected charged particle, an adapter must provide:

- four-momentum in GeV, ordered as `(E, px, py, pz)`
- charge `Q` in units of the positron charge, so a positron has `Q = +1`
- `inOutSign`, with `-1` for incoming particles and `+1` for outgoing
  particles

The formula notation may call this sign `eta_i`, but code and public APIs
should use `inOutSign` to avoid confusion with pseudorapidity.

Neutral particles must not be passed to the Low-current core. Generator
metadata may be retained for diagnostics, but it is not part of the physics
kernel input.

## Default Selection Policy

The default adapter policy is:

- include physical incoming beam particles
- include final stable charged particles
- exclude neutral particles
- exclude unstable intermediate charged particles
- keep `ParticleSelection::IncludeIncomingBeams` enabled by default

Other selections may be exposed through adapter or user policy, but they should
remain explicit options. The default should be conservative and close to the
external charged particles entering the Low-current sum.

## Adapter Responsibilities

Adapters are responsible for:

- identifying incoming particles according to the generator record
- identifying final stable charged particles according to the generator record
- assigning the LowMC `inOutSign`
- converting momenta and charges into LowMC conventions
- applying optional user selection policy
- rejecting unsupported or ambiguous records before they reach the core
- reporting diagnostics about selected and rejected particles

Generator-specific concepts must terminate at the adapter boundary. Examples
include HepMC3 status codes, Pythia8 status codes, event-record history links,
beam-remnant conventions, and documentation-specific particle-state labels.

## Diagnostics

An adapter should be able to report, at least in diagnostic mode:

- number of selected incoming charged particles
- number of selected outgoing charged particles
- number of ignored neutral candidates
- total signed charge sum `sum_i inOutSign_i Q_i`
- total selected ordinary charge sum `sum_i Q_i`
- number of rejected charged candidates
- rejection reasons for charged candidates when practical
- source event identifier if the generator record provides one

The signed charge sum is a useful Low-current sanity check. It should be
reported, not silently repaired, because a nonzero value may reflect either a
real selection choice or an adapter mistake.

The generator-independent support types live in `lowmc/AdapterSupport.h`.
They currently provide `ParticleSelectionConfig`,
`AdapterParticleCandidate`, `AdapterParticleSelectionResult`,
`RejectedChargedCandidate`, and `AdapterSelectionDiagnostic`. The dependency-free
selection helper applies `ParticleSelectionConfig` and the PDG charge helper to
produce cleaned `ChargedParticle` input plus diagnostics. These types are
deliberately independent of HepMC3 and Pythia8.

## Beam Remnants And Special Records

Physical incoming beams are included by default.

Beam remnants, documentation lines, intermediate shower particles, and other
generator-specific bookkeeping records should not be included by default unless
the adapter can justify that they represent external charged particles for the
chosen Low-current calculation.

If such records are exposed, they should be behind explicit selection options
and documented per adapter.

## Out Of Scope For The Core

The core must not depend on:

- HepMC3
- Pythia8
- generator-specific status codes
- generator-specific particle-history traversal
- adapter diagnostics formatting
- event-record I/O

The core may accept the cleaned charged-particle list, photon phase space,
sampling configuration, and random seed.

## Direct Pythia In-Memory Boundary

A direct Pythia adapter is an optional adapter layer around the same core
contract. It should translate an in-memory `Pythia8::Event` into
`AdapterParticleCandidate` records and then call the dependency-free
`selectAdapterParticles(...)` helper.

The Pythia adapter must not add a Pythia dependency to `lowmc_core`.

For each Pythia event, the adapter boundary should provide:

- incoming beam candidates identified from Pythia's event-record convention,
  with `isIncomingBeam = true`
- outgoing final-stable candidates identified from Pythia's final-state
  predicates, with `isFinalStable = true`
- four-momenta converted to LowMC order `(E, px, py, pz)` in GeV
- PDG ids for charge lookup through `chargeTimesThreeFromPDG(...)`
- source labels containing enough Pythia information for diagnostics, such as
  event index, particle index, id, status, and classification

The cleaned core input remains unchanged:

- charge `Q` in units of the positron charge
- `inOutSign = -1` for incoming particles
- `inOutSign = +1` for outgoing particles

The default direct Pythia selection should match the general adapter default:
incoming beams plus final stable charged particles. Beam remnants,
documentation entries, intermediate shower particles, and decayed particles
should be excluded by default and exposed only through explicit diagnostic or
selection options.

The first adapter should use the LowMC PDG charge helper by default. Pythia's
particle-data charge can be reported as an optional cross-check later, but the
first boundary should not make Pythia particle data the charge source of truth.

## Event-Record Output Contract

Event-record writing is adapter-owned. The core may return generated photon
candidates and diagnostics, but it must not know how a generator-specific
record stores particles, vertices, units, status codes, attributes, or barcodes.

Before any adapter writes photons to an event record, the event must have passed
the current generator-independent count-plus-shape validation path:

- valid selected charged-particle input
- accepted `lambda_event` estimate for the event class
- successful Poisson multiplicity draw
- successful shape unweighting with no envelope failure
- enough accepted shape candidates for the requested photon count

If any of these conditions fail, the output adapter must not silently write
ordinary final-state photons. It should either skip writing photons for that
event with diagnostics or fail according to an explicit user policy.

The current diagnostic HepMC3 afterburner may stream directly to the requested
output path. In strict failure mode this means a later failed event can leave a
readable partial diagnostic output containing earlier successful events. Such a
file is not a completed LowMC event sample.

The diagnostic afterburner also exposes an explicit atomic mode. In atomic
mode it writes to `<output>.tmp` and publishes the requested output path only
after the whole input file has completed successfully. Future production-style
HepMC3 output should use this atomic file contract by default.

For physical event-record photons, the default writer contract is:

- write only unweighted physical Low photons
- use enhancement factor `1`
- preserve the original event record
- do not modify input particle momenta, vertices, or event weights
- do not apply recoil
- write photon four-momenta using the event record's momentum-unit convention
- attach provenance metadata identifying LowMC, the leading-Low model, and the
  no-recoil policy

For HepMC3 output, the first durable topology is attributed stable root
photons. Provenance metadata is attached directly to each added photon. A
dedicated no-incoming afterburner vertex did not round-trip cleanly through the
local HepMC3 ASCII writer/reader, so the first convention preserves all
existing vertices and appends attributed stable photons instead. A future
writer may revisit the vertex topology only if a portable HepMC3 convention is
identified and tested.

Enhanced HepMC3 output is an explicit mode, not the physical default. It draws
`N_gamma ~ Poisson(enhancementFactor * lambda_event)`. The durable fallback
metadata location is per added photon, with
`LowMC:photon_mode = "enhanced"`, `LowMC:enhancement_factor = F`, and
`LowMC:weight = 1 / F`. The first compact duplicate convention is event-level
metadata with the same mode, enhancement factor, and, in enhanced mode, weight
on events that receive LowMC photons. Skipped events remain unchanged and do
not receive LowMC event metadata. Weighted histograms recover the physical Low
yield, while direct counting of Low photons in such a file gives the enhanced
yield. Run-level metadata may be added later if a stable file-level
`GenRunInfo` convention is needed.

Diagnostic weighted samples remain separate from enhanced event-record output.
A future diagnostic writer may store weighted samples in an event record only
if the file and particles are explicitly marked as weighted diagnostic output.

## Open Questions

- Which HepMC3 status codes and attributes should identify incoming beams in
  the first adapter?
- Which HepMC3 records should count as final stable charged particles?
- Which Pythia8 record entries should be exposed as incoming particles?
- Which Pythia8 record entries should be exposed as final stable charged
  particles?
- Which Pythia8 status/helper combination should identify beam remnants
  robustly enough for diagnostics?
- Should a future Pythia diagnostic compare LowMC's PDG charge result with
  Pythia particle-data charge and warn on disagreement?
- How should beam remnants be represented in adapter diagnostics?
- Should rejected charged particles be reported by default or only in debug
  diagnostics?
- Should the adapter attach source metadata to selected particles for later
  tracing, and if so, what minimal metadata is generator-independent?
