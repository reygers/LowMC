# NEXT_TASK

This file should contain exactly one immediate task for Codex or another coding
agent. Broader plans belong in `docs/ROADMAP.md`; project state belongs in
`docs/STATUS.md`.

Search for naturally occurring shape-envelope failures or near-failures in
bounded realistic Pythia-produced HepMC3 afterburner validation runs.

Requirements:

- review `docs/ADAPTER_CONTRACT.md`, `docs/ADAPTER_PLAN.md`,
  `docs/DECISIONS.md`, `docs/STATUS.md`, and `docs/USAGE.md`
- use `scripts/validate_pythia_afterburner.py`, `validation_summary.csv`,
  `validation_events.csv`, afterburner logs, and output-summary logs from
  bounded runs
- keep the workflow manual/opt-in; do not add slow default CTest coverage
- run a small bounded seed scan that is likely to include varied realistic
  charged multiplicities but remains practical on a laptop
- keep the default safe envelope settings first; do not start by weakening
  safety factors
- report counts of `shape_status`, `max_trial_bound_ratio`, skipped events,
  and any mismatches with the output-summary diagnostic
- if no natural failures appear, document the observed upper bound for this
  small scan and decide whether a larger manual scan or a targeted stress case
  is the next sensible step
- if a natural failure or near-failure appears, preserve the seed, event
  number, logs, and CSV rows needed to reproduce it
- update `docs/USAGE.md`, `docs/STATUS.md`, and `docs/NEXT_TASK.md` afterward

Open questions:

- How large should the first natural-failure scan be before it becomes too slow
  for ordinary collaborator validation?
- Should near-failures be defined by `max_trial_bound_ratio` alone, or also by
  low acceptance/retained-candidate margin once those fields are exposed?
