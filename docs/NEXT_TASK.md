# NEXT_TASK

This file should contain exactly one immediate task for Codex or another coding
agent. Broader plans belong in `docs/ROADMAP.md`; project state belongs in
`docs/STATUS.md`.

Run a bounded multi-event HepMC3 afterburner validation on a Pythia-produced
sample.

Requirements:

- review `docs/ADAPTER_CONTRACT.md`, `docs/ADAPTER_PLAN.md`,
  `docs/DECISIONS.md`, `docs/STATUS.md`, and `docs/USAGE.md`
- use the existing Pythia-to-HepMC3 sample-generation helper or an equivalent
  bounded build-tree path
- run `lowmc_hepmc3_afterburner_diagnostic` on more than one generated event
  with explicit photon phase space, seed, and enhancement settings
- keep the validation bounded enough for manual use first; do not add a slow
  default CTest path
- summarize event count, skipped events, added LowMC photons, photon weights,
  and runtime-relevant observations
- inspect whether any naturally occurring selection, yield, or shape-envelope
  failures appear in the generated sample
- compare output-summary diagnostics against the afterburner text diagnostics
  rather than relying on visual inspection
- keep this as validation evidence for the diagnostic afterburner, not as a
  production-readiness claim
- update `docs/STATUS.md` and `docs/NEXT_TASK.md` afterward

Open questions:

- How many generated events and which enhancement factor give useful evidence
  without making the manual validation slow?
- Should this become a helper script after the first manual run, or remain a
  documented command sequence?
