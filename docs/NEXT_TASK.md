# NEXT_TASK

This file should contain exactly one immediate task for Codex or another coding
agent. Broader plans belong in `docs/ROADMAP.md`; project state belongs in
`docs/STATUS.md`.

Prototype a reduced or targeted base-probe lattice behind a full-reference
comparison.

Requirements:

- review `docs/ADAPTER_CONTRACT.md`, `docs/ADAPTER_PLAN.md`,
  `docs/DECISIONS.md`, `docs/STATUS.md`, `docs/ALGORITHM.md`, and
  `docs/USAGE.md`
- keep the factorized event-yield estimator and the spatial-indexed
  shape-unweighting path as the current reference defaults; keep the
  known-cell interior probe shortcut enabled because it is semantics-preserving
- do not change physics densities, `max_trial_bound_ratio` accounting, or
  failure policy semantics
- treat the current full `probeBins x probeBins` base lattice as the reference
  envelope construction
- prototype a reduced or targeted base-probe lattice only behind an explicit
  diagnostic/reference comparison; do not silently replace the full lattice
- compare reference and prototype bounds per cell where practical, and at least
  compare event-level envelope integral, accepted/retained photon counts,
  envelope failures, `max_trial_bound_ratio`, suspect/promoted-cell counts, and
  substage timings
- use controlled toy cases plus the bounded Pythia-produced seeds `24681`,
  `24682`, `24683`, and `24684`
- require no new envelope failures, unchanged accepted/retained photon counts
  for fixed seeds, and no larger accepted `max_trial_bound_ratio` before
  considering the prototype as a candidate for broader validation
- report speed impact using `shape_base_probe_seconds`, total `shape_seconds`,
  and total `event_seconds`
- update `docs/STATUS.md`, `docs/ALGORITHM.md`, and `docs/NEXT_TASK.md`
  afterward

Open questions:

- Is a reduced lattice with calibration/promoted-cell safeguards conservative
  enough on realistic events, or does it need an explicit compensation factor?
- Should the first prototype reduce probes uniformly, or target full probing to
  high-risk cells near charged-particle directions and calibration suspects?
- How much per-cell bound agreement is necessary before the prototype can move
  from diagnostic comparison to an optional generation mode?
