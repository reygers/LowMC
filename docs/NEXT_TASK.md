# NEXT_TASK

This file should contain exactly one immediate task for Codex or another coding
agent. Broader plans belong in `docs/ROADMAP.md`; project state belongs in
`docs/STATUS.md`.

Audit and prototype a faster event-yield estimator using the leading-Low
`k_T` factorization.

Requirements:

- review `docs/ADAPTER_CONTRACT.md`, `docs/ADAPTER_PLAN.md`,
  `docs/DECISIONS.md`, `docs/STATUS.md`, `docs/ALGORITHM.md`, and
  `docs/USAGE.md`
- keep the spatial-indexed shape-unweighting path unchanged; this task is about
  the yield estimator now visible in stage timing
- use the leading-Low scaling in which the real-photon angular shape is
  independent of `k_T` and the integrated yield factorizes into an angular
  integral times `ln(kTMax / kTMin)`
- compare the current adaptive `log k_T` midpoint estimator with a factorized
  estimator that builds/integrates the angular grid once at a representative
  `k_T`
- preserve the existing estimator as a reference while validating numerical
  agreement on toy and bounded Pythia-produced events
- decide whether the afterburner can use the factorized estimator by default,
  or whether it should remain a diagnostic option until broader validation
- report any speed impact using `yield_seconds` and `event_seconds`
- update `docs/STATUS.md`, `docs/ALGORITHM.md`, and `docs/NEXT_TASK.md`
  afterward

Open questions:

- Does the existing adaptive-yield code intentionally keep multiple `k_T`
  slices for future non-leading effects, or can the leading-Low afterburner use
  an explicitly factorized path now?
- Which angular grid configuration gives agreement with the current reference
  without rebuilding maps at many `k_T` midpoints?
- Should the currently unused coarse-yield comparison be made optional in the
  afterburner path as part of the same change, or handled separately?
