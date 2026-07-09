# NEXT_TASK

This file should contain exactly one immediate task for Codex or another coding
agent. Broader plans belong in `docs/ROADMAP.md`; project state belongs in
`docs/STATUS.md`.

Broaden validation of the shape-unweighting spatial-index prototype and measure
its speed impact on a bounded Pythia-produced scan.

Requirements:

- review `docs/ADAPTER_CONTRACT.md`, `docs/ADAPTER_PLAN.md`,
  `docs/DECISIONS.md`, `docs/STATUS.md`, `docs/ALGORITHM.md`, and
  `docs/USAGE.md`
- treat the current spatial index as a bookkeeping optimization; do not change
  physics densities, envelope safety factors, or yield-integration settings
- rerun a bounded multi-seed Pythia-produced afterburner scan with stage timing
  and compare against the pre-index timing notes in `docs/STATUS.md`
- include at least the known seeds `24681`, `24683`, and `24684`; add more
  seeds only if runtime remains practical
- check that photon counts, `shape_status`, `max_trial_bound_ratio`, skipped
  events, and output-summary agreement remain consistent with the pre-index
  validation expectations
- report `event_seconds`, `shape_seconds`, and `yield_seconds` before/after
  where direct comparison data exists
- decide whether the spatial index is ready to keep as the default internal
  path or whether it needs a runtime switch/reference mode for more testing
- update `docs/STATUS.md`, `docs/ALGORITHM.md`, and `docs/NEXT_TASK.md`
  afterward

Open questions:

- How many realistic Pythia seeds are enough before calling this a stable
  bookkeeping speedup?
- Should the diagnostic consistency helper remain public, or should it move to
  a narrower test-support interface later?
- After the spatial-index speedup, is yield integration now the next meaningful
  bottleneck?
