# NEXT_TASK

This file should contain exactly one immediate task for Codex or another coding
agent. Broader plans belong in `docs/ROADMAP.md`; project state belongs in
`docs/STATUS.md`.

Wait for independent review feedback on the reduced-probe recommendation.

Requirements:

- review `docs/STATUS.md`, `docs/ALGORITHM.md`, `docs/DECISIONS.md`, and
  `docs/USAGE.md`
- keep the factorized event-yield estimator, spatial-indexed shape-unweighting
  path, known-cell interior probe shortcut, and full `11 x 11` boosted-cell
  reference as the current defaults
- do not change physics densities, `max_trial_bound_ratio` accounting,
  multiplicity conventions, output metadata, or failure policy semantics
- keep reduced probing and compensation diagnostic-only; do not make
  `5 x 5` with any compensation factor a default
- send or discuss `docs/reduced-probe-review-note-stefan.md` with Stefan
- ask Stefan to focus on:
  - whether the no-per-seed-overshoot criterion relative to the full `11 x 11`
    reference is the right current safety criterion
  - whether factor `1.50`'s accepted-trial penalty is acceptable for a
    non-default diagnostic option
  - whether the next validation should use more fresh Pythia seeds, different
    Pythia processes, or targeted topology classes
  - whether additional diagnostics are needed to identify calibration,
    boosted-cell-promotion, or residual reduced-probe under-bounding failures
- do not run more validation or change defaults until review feedback is
  available
- update `docs/STATUS.md`, `docs/ALGORITHM.md`, and `docs/NEXT_TASK.md`
  afterward

Open questions:

- Does Stefan agree that factor `1.50` is the only current conservative
  reduced-probe candidate?
- What validation should be required before reduced probing becomes a
  documented non-default diagnostic option?
- Should the next validation broaden Pythia seeds, process types, or topology
  classes?
