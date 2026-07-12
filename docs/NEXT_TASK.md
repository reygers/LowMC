# NEXT_TASK

This file should contain exactly one immediate task for Codex or another coding
agent. Broader plans belong in `docs/ROADMAP.md`; project state belongs in
`docs/STATUS.md`.

Prototype a reduced-probe compensation factor.

Requirements:

- review `docs/STATUS.md`, `docs/ALGORITHM.md`, `docs/DECISIONS.md`, and
  `docs/USAGE.md`
- keep the factorized event-yield estimator, spatial-indexed shape-unweighting
  path, known-cell interior probe shortcut, and full `11 x 11` base-probe
  reference as current defaults
- do not change physics densities, `max_trial_bound_ratio` accounting,
  multiplicity conventions, or failure policy semantics
- keep the reduced/full audit and sparse severe-cell report diagnostic-only;
  do not make reduced probing a default
- start from the sparse audit result: severe cells are stable across `5 x 5`
  and `7 x 7`, sparse, low proposal-mass, mostly not calibration-suspect, and
  only partly explained by charged-particle proximity
- add an explicit prototype option that multiplies the probed maximum for
  reduced-probed cells by a compensation factor before envelope construction
- do not apply the compensation to cells already probed with the full
  `11 x 11` base lattice or to boosted-cell re-probes unless the code path
  makes that behavior explicit
- keep the current full-probe reference unchanged and report the compensation
  factor in diagnostics
- scan a small set of factors, for example `1.05`, `1.10`, `1.20`, and `1.30`,
  on bounded seeds `24681`, `24682`, `24683`, and `24684`, comparing `5 x 5`
  and `7 x 7` with boosted-cell probes disabled
- compare runtime, `max_trial_bound_ratio`, and accepted photon counts against
  the full `11 x 11` reference and the uncompensated reduced-probe runs
- update `docs/STATUS.md`, `docs/ALGORITHM.md`, and `docs/NEXT_TASK.md`
  afterward

Open questions:

- Is a single global reduced-probe compensation factor enough for both `5 x 5`
  and `7 x 7`, or should any factor depend on the reduced lattice size?
- Should the compensation be applied only to uniformly reduced cells, or also
  to cells reduced after targeted/evidence-based promotion?
- What criterion should be required before considering a compensated reduced
  probe mode beyond diagnostics: matching the full-reference
  `max_trial_bound_ratio`, no envelope failures, or a stronger safety margin?
