# NEXT_TASK

This file should contain exactly one immediate task for Codex or another coding
agent. Broader plans belong in `docs/ROADMAP.md`; project state belongs in
`docs/STATUS.md`.

Prototype a reduced-versus-full base-probe audit.

Requirements:

- review `docs/STATUS.md`, `docs/ALGORITHM.md`, `docs/DECISIONS.md`, and
  `docs/USAGE.md`
- keep the factorized event-yield estimator, spatial-indexed shape-unweighting
  path, known-cell interior probe shortcut, and full `11 x 11` base-probe
  reference as current defaults
- do not change physics densities, `max_trial_bound_ratio` accounting,
  multiplicity conventions, or failure policy semantics
- keep the audit diagnostic-only; do not make reduced probing a default
- start from the current evidence: charged-particle proximity, calibration
  suspect cells, and suspicious-cell neighborhoods all fail to reproduce the
  full-probe reference on the bounded Pythia high-multiplicity seed
- add an explicit diagnostic mode that evaluates both a reduced base lattice
  and the full `11 x 11` lattice for audit purposes, then reports where the
  reduced lattice underestimates the full probe maximum
- report per-event audit summaries such as audited cell count, number/fraction
  of under-bounded cells, maximum full/reduced probe ratio, and optional
  severity buckets
- keep the normal envelope construction separate from the audit result unless
  an explicit prototype option uses the audit to choose full-probed cells
- compare the audit on the bounded Pythia-produced seeds `24681`, `24682`,
  `24683`, and `24684`
- use the audit result to recommend the next targeted policy, for example a
  severity threshold, a cell-size or mass predictor, or abandoning reduced base
  probes if the under-bounded cells are too diffuse
- update `docs/STATUS.md`, `docs/ALGORITHM.md`, and `docs/NEXT_TASK.md`
  afterward

Open questions:

- Are under-bounded cells localized enough that a targeted rule can catch them?
- Does reduced/full disagreement correlate with cell size, proposal mass,
  calibration evidence, charged-particle proximity, or another cheap feature?
- Should the audit drive a new promotion rule, or should LowMC keep full base
  probing until a different envelope construction is available?
