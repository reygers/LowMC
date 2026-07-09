# NEXT_TASK

This file should contain exactly one immediate task for Codex or another coding
agent. Broader plans belong in `docs/ROADMAP.md`; project state belongs in
`docs/STATUS.md`.

Instrument the remaining shape-unweighting bottleneck inside the HepMC3
afterburner diagnostic.

Requirements:

- review `docs/ADAPTER_CONTRACT.md`, `docs/ADAPTER_PLAN.md`,
  `docs/DECISIONS.md`, `docs/STATUS.md`, `docs/ALGORITHM.md`, and
  `docs/USAGE.md`
- keep the factorized event-yield estimator and the spatial-indexed
  shape-unweighting path as the current defaults; do not change physics
  densities or envelope safety policy in this task
- add timing and counter diagnostics that split the current `shape_seconds`
  stage into the main pieces of `generateShapeUnweightedPhotonCandidates(...)`,
  such as grid construction, envelope probing/calibration, trial generation,
  accepted-candidate bookkeeping, and envelope-failure checks
- make the new diagnostics available to the HepMC3 afterburner timing summary
  without cluttering ordinary user-facing output more than necessary
- rerun the bounded four-seed Pythia-produced timing scan and identify the
  dominant remaining substage for seeds `24681`, `24682`, `24683`, and `24684`
- preserve existing shape-envelope diagnostics, especially
  `max_trial_bound_ratio`, accepted/retained counts, and failure status
- use the result to choose the next optimization target, but avoid implementing
  that optimization in the same step unless it is only a trivial bookkeeping
  cleanup
- update `docs/STATUS.md`, `docs/ALGORITHM.md`, and `docs/NEXT_TASK.md`
  afterward

Open questions:

- Is the current shape cost dominated by envelope setup/probing or by the
  number of rejection trials needed for accepted photons?
- Should shape diagnostics live directly in `ShapeUnweightingResult`, or should
  they be exposed through a narrower diagnostic-only structure?
- Once the dominant substage is known, is the next low-risk optimization a
  bookkeeping/cache improvement, or does it require a sampling-policy change
  that needs broader validation?
