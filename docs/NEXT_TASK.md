# NEXT_TASK

This file should contain exactly one immediate task for Codex or another coding
agent. Broader plans belong in `docs/ROADMAP.md`; project state belongs in
`docs/STATUS.md`.

Prototype a faster base-probe stage for shape-unweighting envelope
construction.

Requirements:

- review `docs/ADAPTER_CONTRACT.md`, `docs/ADAPTER_PLAN.md`,
  `docs/DECISIONS.md`, `docs/STATUS.md`, `docs/ALGORITHM.md`, and
  `docs/USAGE.md`
- keep the factorized event-yield estimator and the spatial-indexed
  shape-unweighting path as the current reference defaults; do not change
  physics densities, `max_trial_bound_ratio` accounting, or failure policy
  semantics
- use the new shape substage timings as the motivation: base probing is about
  `75%` of shape time in the current four-seed Pythia-produced scan, while
  final rejection trials are small
- inspect why the base probe pass evaluates every cell on the full base
  `probeBins x probeBins` lattice, and identify which parts are unavoidable
  envelope-safety work versus bookkeeping overhead
- prototype one low-risk faster base-probe path behind an explicit option or
  diagnostic/reference comparison, not as an unvalidated silent replacement
- compare the prototype against the current full base-probe reference on toy
  cases and the bounded Pythia-produced seeds `24681`, `24682`, `24683`, and
  `24684`
- require unchanged or more conservative envelope diagnostics before considering
  the prototype as a future default: no new envelope failures, no larger
  accepted `max_trial_bound_ratio`, and unchanged accepted/retained photon
  counts for fixed seeds
- report speed impact using the shape substage timings, especially
  `shape_base_probe_seconds`, total `shape_seconds`, and total `event_seconds`
- update `docs/STATUS.md`, `docs/ALGORITHM.md`, and `docs/NEXT_TASK.md`
  afterward

Open questions:

- Can any base-probe bookkeeping be removed without changing the set of probe
  coordinates or the exact evaluated weights?
- If reducing the number of base probes is necessary for a meaningful speedup,
  what validation would make that acceptable, given envelope safety?
- Should a reduced-probe prototype be combined with calibration/promoted-cell
  safeguards, or should it first be tested as a pure reference-vs-prototype
  diagnostic?
