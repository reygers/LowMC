# NEXT_TASK

This file should contain exactly one immediate task for Codex or another coding
agent. Broader plans belong in `docs/ROADMAP.md`; project state belongs in
`docs/STATUS.md`.

Validate compensated candidates with boosted-cell probes.

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
- start from the ten-seed no-boosted result: all candidates stayed below the
  full-reference ratio, but `7 x 7` with factor `1.10` had only a small margin;
  `7 x 7` with `1.20` is the safer candidate and `5 x 5` with `1.30` is the
  aggressive speed candidate
- run the same bounded seeds `24681..24690` with the default boosted-cell
  probe path enabled for the full `11 x 11` reference, `7 x 7` with
  compensation factor `1.20`, and `5 x 5` with compensation factor `1.30`
- include `7 x 7` with factor `1.10` only if the first two candidates leave
  enough ambiguity to justify the extra runs
- compare total event time, shape time, `max_trial_bound_ratio`, envelope
  failures, accepted/retained photon counts, compensated-cell counts, promoted
  cell counts, and boosted-probe time
- decide whether boosted-cell probes still add useful safety on top of
  compensation, or whether they mostly erase the speed advantage
- update `docs/STATUS.md`, `docs/ALGORITHM.md`, and `docs/NEXT_TASK.md`
  afterward

Open questions:

- With boosted-cell probes enabled, does compensation still affect enough cells
  to matter?
- Does the boosted-cell path reduce `max_trial_bound_ratio` further, or only
  increase runtime for these candidates?
- Should the eventual validation criterion require both no-boosted and boosted
  scans to beat the full-reference ratio, or is the no-boosted comparison the
  cleaner criterion for the reduced-base-probe policy?
