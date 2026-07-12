# NEXT_TASK

This file should contain exactly one immediate task for Codex or another coding
agent. Broader plans belong in `docs/ROADMAP.md`; project state belongs in
`docs/STATUS.md`.

Characterize severe reduced/full base-probe audit cells.

Requirements:

- review `docs/STATUS.md`, `docs/ALGORITHM.md`, `docs/DECISIONS.md`, and
  `docs/USAGE.md`
- keep the factorized event-yield estimator, spatial-indexed shape-unweighting
  path, known-cell interior probe shortcut, and full `11 x 11` base-probe
  reference as current defaults
- do not change physics densities, `max_trial_bound_ratio` accounting,
  multiplicity conventions, or failure policy semantics
- keep the reduced/full audit diagnostic-only; do not make reduced probing a
  default
- start from the current audit result: strict `full > reduced` occurs in nearly
  every cell, while full/reduced ratios above `1.1` are sparse, especially for
  `7 x 7` base probes
- add a compact per-event or optional sparse per-cell diagnostic for severe
  audit cells, for example cells with full/reduced ratio above `1.1`
- report cheap candidate predictors for severe cells, such as cell eta/phi
  center and size, cell area, proposal mass, reduced-probe maximum,
  calibration maximum, calibration-suspect status, and charged-particle
  proximity if already available without a costly new search
- run the diagnostic on the bounded Pythia-produced seeds `24681`, `24682`,
  `24683`, and `24684`, comparing at least `5 x 5` and `7 x 7` reduced base
  probes with boosted-cell probes disabled
- decide whether severe audit cells are localized or predictable enough for a
  targeted full-probe promotion rule, or whether reduced base probing should
  remain only a diagnostic until a different envelope construction exists
- update `docs/STATUS.md`, `docs/ALGORITHM.md`, and `docs/NEXT_TASK.md`
  afterward

Open questions:

- Which severity threshold is most useful for policy design: `>1.1`, `>1.5`,
  or a threshold tied to the observed `max_trial_bound_ratio` budget?
- Are severe cells concentrated in large cells, high proposal-mass cells, cells
  near charged-particle directions, or cells where calibration already hints
  that the reduced probe missed structure?
- Is a sparse severe-cell report enough, or is a small standalone analysis
  script needed to compare predictor distributions against all audited cells?
