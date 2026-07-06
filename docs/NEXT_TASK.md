# NEXT_TASK

This file should contain exactly one immediate task for Codex or another coding
agent. Broader plans belong in `docs/ROADMAP.md`; project state belongs in
`docs/STATUS.md`.

Broaden the direct Pythia versus Pythia-to-HepMC3 cross-check beyond one event.

Requirements:

- review `docs/ADAPTER_CONTRACT.md`, `docs/ADAPTER_PLAN.md`,
  `docs/DECISIONS.md`, `docs/STATUS.md`, and `docs/USAGE.md`
- reuse `lowmc_pythia8_hepmc3_crosscheck`
- support checking more than one accepted Pythia event and/or more than one
  seed in a bounded smoke configuration
- keep the default CTest path quick
- compare the same selected-particle diagnostics and weighted summary fields
  already used in the one-event cross-check
- print an aggregate summary with event count, mismatch count, and maximum
  relative weighted-summary difference
- document any expected differences caused by particle ordering or
  Pythia-to-HepMC3 conversion rather than hiding them
- keep this as a diagnostic/cross-check target, not a production dependency
- add focused CTest coverage only when both optional adapters are enabled
- update `docs/STATUS.md` and `docs/NEXT_TASK.md` afterward

Open questions:

- Should the first broader check use several events from one seed or one event
  from several seeds?
- What maximum runtime should be allowed for this optional CTest smoke check?
