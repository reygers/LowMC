# NEXT_TASK

This file should contain exactly one immediate task for Codex or another coding
agent. Broader plans belong in `docs/ROADMAP.md`; project state belongs in
`docs/STATUS.md`.

Run a bounded Pythia-produced HepMC3 afterburner timing scan using true
per-event `event_seconds` and summarize the runtime/complexity picture.

Requirements:

- review `docs/ADAPTER_CONTRACT.md`, `docs/ADAPTER_PLAN.md`,
  `docs/DECISIONS.md`, `docs/STATUS.md`, and `docs/USAGE.md`
- use `scripts/validate_pythia_afterburner.py`, `validation_summary.csv`,
  `validation_events.csv`, and the afterburner logs from bounded runs
- keep the workflow manual/opt-in; do not add slow default CTest coverage
- run a small multi-seed or multi-event scan that is large enough to exercise
  several charged multiplicities but still practical on a laptop
- use the new true `event_seconds` field rather than the helper's older
  evenly distributed timing estimate
- summarize, at minimum, runtime versus selected outgoing charged multiplicity,
  requested photons, retained photons, and `max_trial_bound_ratio`
- decide whether a small aggregation script is useful or whether the existing
  CSVs plus a documented manual summary are enough
- do not tune physics or sampling parameters in this step unless the scan
  exposes a clear failure
- update `docs/USAGE.md`, `docs/STATUS.md`, and `docs/NEXT_TASK.md` afterward

Open questions:

- What bounded scan size gives useful runtime spread without turning manual
  validation into a long benchmark?
- Should future runtime summaries report simple correlations, worst-event
  rows, or both?
