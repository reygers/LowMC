# NEXT_TASK

This file should contain exactly one immediate task for Codex or another coding
agent. Broader plans belong in `docs/ROADMAP.md`; project state belongs in
`docs/STATUS.md`.

Add a helper script for bounded Pythia-produced HepMC3 afterburner validation.

Requirements:

- review `docs/ADAPTER_CONTRACT.md`, `docs/ADAPTER_PLAN.md`,
  `docs/DECISIONS.md`, `docs/STATUS.md`, and `docs/USAGE.md`
- wrap the existing `scripts/generate_pythia_hepmc3_sample.py`,
  `lowmc_hepmc3_afterburner_diagnostic`, and
  `lowmc_hepmc3_output_summary_diagnostic` workflow
- keep it opt-in/manual; do not add the full validation run to default CTest
- expose bounded defaults for event count, enhancement factor, seeds, and
  photon phase space
- write generated input, afterburner output, afterburner log, and summary log
  under a user-selected output directory
- parse the afterburner and output-summary logs enough to report event count,
  skipped events, added LowMC photons, photon weight sums, shape statuses, and
  mismatch warnings
- print runtime for the generated-sample step and afterburner step
- add only a cheap `--help` test if CTest coverage is added
- update `docs/USAGE.md`, `docs/STATUS.md`, and `docs/NEXT_TASK.md` afterward

Open questions:

- Should the helper fail on summary mismatches by default, or only report them
  until the log parser is more mature?
- Which output format is most useful for manual logs: compact text only, CSV,
  or both?
