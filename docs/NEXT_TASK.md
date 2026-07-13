# NEXT_TASK

This file should contain exactly one immediate task for Codex or another coding
agent. Broader plans belong in `docs/ROADMAP.md`; project state belongs in
`docs/STATUS.md`.

Validate compensated reduced-probe candidates.

Requirements:

- review `docs/STATUS.md`, `docs/ALGORITHM.md`, `docs/DECISIONS.md`, and
  `docs/USAGE.md`
- keep the factorized event-yield estimator, spatial-indexed shape-unweighting
  path, known-cell interior probe shortcut, and full `11 x 11` base-probe
  reference as current defaults
- do not change physics densities, `max_trial_bound_ratio` accounting,
  multiplicity conventions, or failure policy semantics
- keep compensated reduced probing diagnostic-only; do not make reduced probing
  or compensation a default
- start from the four-seed result: `5 x 5` with factor `1.30` and `7 x 7`
  with factors `1.10` or `1.20` are the first plausible candidates
- run a broader bounded Pythia seed sample for these candidates and for the
  full `11 x 11` reference, with boosted-cell probes disabled
- include the uncompensated `5 x 5` and `7 x 7` runs only if needed for a clean
  comparison table
- compare total event time, shape time, `max_trial_bound_ratio`, envelope
  failures, accepted/retained photon counts, and compensated-cell counts
- decide whether `7 x 7` with `1.10` is already sufficient or whether the
  safer candidate is `7 x 7` with `1.20`; keep `5 x 5` with `1.30` as the
  aggressive speed candidate
- if the broader sample still looks promising, update the docs with a candidate
  validation criterion for future compensated reduced-probe studies
- update `docs/STATUS.md`, `docs/ALGORITHM.md`, and `docs/NEXT_TASK.md`
  afterward

Open questions:

- How many bounded seeds are enough for this diagnostic phase: `10`, `20`, or
  more?
- Should the comparison criterion be matching the full-reference
  `max_trial_bound_ratio`, staying below it by a margin, or using a direct
  full/reduced audit on the same events?
- Should compensated reduced probing be studied with boosted-cell probes
  enabled after the no-boosted validation, or only after the no-boosted case is
  clearly stable?
