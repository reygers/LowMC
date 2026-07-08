# NEXT_TASK

This file should contain exactly one immediate task for Codex or another coding
agent. Broader plans belong in `docs/ROADMAP.md`; project state belongs in
`docs/STATUS.md`.

Decide whether true per-event afterburner timing is needed for bounded
Pythia-produced HepMC3 validation, and add it only if the helper-level estimate
is not sufficient.

Requirements:

- review `docs/ADAPTER_CONTRACT.md`, `docs/ADAPTER_PLAN.md`,
  `docs/DECISIONS.md`, `docs/STATUS.md`, and `docs/USAGE.md`
- use `scripts/validate_pythia_afterburner.py`, `validation_summary.csv`,
  `validation_events.csv`, and the afterburner logs from bounded runs
- keep the workflow manual/opt-in; do not add slow default CTest coverage
- inspect whether the current per-event complexity CSV plus aggregate runtime
  is enough to diagnose slow bounded Pythia-produced events
- if it is enough, document the decision and choose the next validation task
- if it is not enough, add narrowly scoped timing around per-event work in
  `lowmc_hepmc3_afterburner_diagnostic` and include an `event_seconds`-style
  key in the existing compact text output
- keep any new output backwards-compatible key-value text
- add only cheap parser or `--help` tests unless the C++ timing change itself
  needs a focused fixture test
- update `docs/USAGE.md`, `docs/STATUS.md`, and `docs/NEXT_TASK.md` afterward

Open questions:

- Should timing include HepMC3 reading and writing, or only LowMC event
  selection, count sampling, and shape sampling?
- Should true per-event timing always be printed by the diagnostic, or only
  behind an explicit diagnostic flag?
