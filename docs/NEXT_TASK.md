# NEXT_TASK

This file should contain exactly one immediate task for Codex or another coding
agent. Broader plans belong in `docs/ROADMAP.md`; project state belongs in
`docs/STATUS.md`.

Prototype an evidence-based targeted base-probe policy.

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
- start from the current evidence: uniform `5 x 5` and `7 x 7` reductions are
  fast but raise `max_trial_bound_ratio`, and charged-particle-proximity
  targeting alone is not sufficient
- design a two-pass diagnostic policy that first uses a reduced base lattice,
  then promotes cells to full probing using evidence from the reduced probe
  pass and/or calibration
- candidate high-risk evidence may include calibration exceeding the reduced
  probe maximum, large reduced-probe variation within a cell, high estimated
  envelope contribution, small cells near sharp structures, or neighboring
  cells with suspicious reduced-probe behavior
- keep the comparison against the full-probe reference explicit: accepted and
  retained photon counts, envelope failures, `max_trial_bound_ratio`,
  suspect/promoted-cell counts, probe evaluations, full-probed-cell counts,
  and substage timings
- run controlled toy checks plus the bounded Pythia-produced seeds `24681`,
  `24682`, `24683`, and `24684`
- require no new envelope failures, unchanged accepted/retained photon counts
  for fixed seeds, and no larger accepted `max_trial_bound_ratio` before
  treating the prototype as a candidate for broader validation
- update `docs/STATUS.md`, `docs/ALGORITHM.md`, and `docs/NEXT_TASK.md`
  afterward

Open questions:

- Which reduced-pass evidence best predicts cells where full probing matters?
- Is a separate compensation/safety factor still needed after evidence-based
  promotion?
- How much speedup remains once the prototype matches the full-probe reference
  on the bounded Pythia seeds?
