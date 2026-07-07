# NEXT_TASK

This file should contain exactly one immediate task for Codex or another coding
agent. Broader plans belong in `docs/ROADMAP.md`; project state belongs in
`docs/STATUS.md`.

Run a small bounded multi-seed Pythia-produced HepMC3 afterburner scan with the
validation helper.

Requirements:

- review `docs/ADAPTER_CONTRACT.md`, `docs/ADAPTER_PLAN.md`,
  `docs/DECISIONS.md`, `docs/STATUS.md`, and `docs/USAGE.md`
- use `scripts/validate_pythia_afterburner.py`
- keep the scan manual/opt-in; do not add slow default CTest coverage
- choose a small set of Pythia seeds and event counts that runs in reasonable
  interactive time
- record event counts, skipped events, added LowMC photons, photon weight sums,
  shape-status counts, maximum trial/bound ratios, and runtimes
- inspect whether any natural selection, yield, or shape-envelope failures
  appear without forcing unsafe diagnostic knobs
- decide whether the helper's current text plus CSV summary is sufficient for
  manual validation logs or whether a multi-run aggregator is needed
- update `docs/USAGE.md`, `docs/STATUS.md`, and `docs/NEXT_TASK.md` afterward

Open questions:

- How many seeds/events are enough to be informative without burning too much
  time?
- Should natural shape-envelope failures in this scan become a separate focused
  failure-mode task?
