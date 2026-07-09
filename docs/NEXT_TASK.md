# NEXT_TASK

This file should contain exactly one immediate task for Codex or another coding
agent. Broader plans belong in `docs/ROADMAP.md`; project state belongs in
`docs/STATUS.md`.

Prototype a targeted base-probe policy for high-risk versus low-risk cells.

Requirements:

- review `docs/ADAPTER_CONTRACT.md`, `docs/ADAPTER_PLAN.md`,
  `docs/DECISIONS.md`, `docs/STATUS.md`, `docs/ALGORITHM.md`, and
  `docs/USAGE.md`
- keep the factorized event-yield estimator, spatial-indexed shape-unweighting
  path, and known-cell interior probe shortcut as current defaults
- do not change physics densities, `max_trial_bound_ratio` accounting, or
  failure policy semantics
- treat the current full `probeBins x probeBins` base lattice as the reference
  envelope construction
- start from the evidence that uniform `5 x 5` and `7 x 7` base-probe
  reductions are fast but raise the observed `max_trial_bound_ratio`
- prototype a targeted policy that keeps the full `11 x 11` base lattice for
  high-risk cells and uses a reduced lattice only for low-risk cells
- define high-risk cells using information available before the base-probe pass,
  such as proximity to charged-particle directions, refined-cell size, or
  proposal/target mass; if calibration-suspect information is needed, keep it
  in a separate two-pass diagnostic rather than silently depending on it
- keep the targeted policy behind an explicit diagnostic/prototype option; do
  not make it the default
- compare reference and prototype diagnostics, including event-level envelope
  integral, accepted/retained photon counts, envelope failures,
  `max_trial_bound_ratio`, suspect/promoted-cell counts, probe evaluations, and
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

- Which pre-base-probe high-risk marker best predicts where full probing is
  needed?
- Does a targeted reduced lattice still need an explicit compensation factor?
- How much per-cell bound agreement is necessary before the prototype can move
  from diagnostic comparison to an optional generation mode?
