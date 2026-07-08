# NEXT_TASK

This file should contain exactly one immediate task for Codex or another coding
agent. Broader plans belong in `docs/ROADMAP.md`; project state belongs in
`docs/STATUS.md`.

Add runtime/complexity reporting for bounded Pythia-produced HepMC3
afterburner validation.

Requirements:

- review `docs/ADAPTER_CONTRACT.md`, `docs/ADAPTER_PLAN.md`,
  `docs/DECISIONS.md`, `docs/STATUS.md`, and `docs/USAGE.md`
- use `scripts/validate_pythia_afterburner.py` and its existing logs/CSV
- keep the workflow manual/opt-in; do not add slow default CTest coverage
- expose or derive enough per-event information to relate afterburner runtime
  to selected charged multiplicity, requested photon count, retained photon
  count, and shape-unweighting diagnostics
- decide whether to extend the helper's `validation_summary.csv`, add a
  per-event CSV, or add a small multi-run aggregation script
- preserve the current compact text summary for quick manual logs
- add only cheap `--help` or parser tests if new scripts are introduced
- update `docs/USAGE.md`, `docs/STATUS.md`, and `docs/NEXT_TASK.md` afterward

Open questions:

- Is total afterburner runtime enough for now, or do we need per-event timing
  inside `lowmc_hepmc3_afterburner_diagnostic`?
- Which complexity proxy is most useful first: selected outgoing charged
  multiplicity, requested photon count, retained candidate count, or trial
  count?
