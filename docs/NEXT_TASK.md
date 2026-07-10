# NEXT_TASK

This file should contain exactly one immediate task for Codex or another coding
agent. Broader plans belong in `docs/ROADMAP.md`; project state belongs in
`docs/STATUS.md`.

Prototype a neighbor-aware evidence base-probe policy.

Requirements:

- review `docs/STATUS.md`, `docs/ALGORITHM.md`, `docs/DECISIONS.md`, and
  `docs/USAGE.md`
- keep the factorized event-yield estimator, spatial-indexed shape-unweighting
  path, known-cell interior probe shortcut, and full `11 x 11` base-probe
  reference as current defaults
- do not change physics densities, `max_trial_bound_ratio` accounting,
  multiplicity conventions, or failure policy semantics
- keep all reduced or targeted base-probe modes behind explicit diagnostic
  options; do not make them defaults
- start from the current evidence: charged-particle proximity alone is not
  sufficient, and calibration-suspect-cell promotion alone is too sparse
- extend the evidence-targeted diagnostic so suspicious cells can promote
  nearby cells to the full `11 x 11` base lattice
- define the neighborhood rule explicitly, for example by cell adjacency,
  eta-phi distance in units of local cell size, or reuse of the existing
  `cellNearCoordinate` style with suspicious-cell centers
- compare at least two neighborhood strengths against the full-probe reference
  on the bounded Pythia-produced seeds `24681`, `24682`, `24683`, and `24684`
- record accepted and retained photon counts, envelope failures,
  `max_trial_bound_ratio`, suspect/promoted/full-probed-cell counts, probe
  evaluations, and substage timings
- require no new envelope failures, unchanged accepted/retained photon counts
  for fixed seeds, and no larger accepted `max_trial_bound_ratio` before
  treating the prototype as a candidate for broader validation
- if neighbor expansion cannot satisfy the criterion without losing the speed
  advantage, propose a direct reduced-versus-full probe audit as the next
  diagnostic instead of adding more tuning knobs
- update `docs/STATUS.md`, `docs/ALGORITHM.md`, and `docs/NEXT_TASK.md`
  afterward

Open questions:

- Which neighborhood definition best catches under-bounded nearby cells without
  collapsing back to full probing?
- Should suspicious-cell expansion depend on reduced-probe severity, such as
  `calibrationMax / reducedProbeMax`, rather than a fixed radius?
- Would a small direct audit of reduced-versus-full probe maxima be more
  informative than another heuristic promotion rule?
