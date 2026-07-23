# NEXT_TASK

This file should contain exactly one immediate task for Codex or another coding
agent. Broader plans belong in `docs/ROADMAP.md`; project state belongs in
`docs/STATUS.md`.

Define a release-blocking validation matrix for the first public LowMC
afterburner workflow.

Requirements:

- inventory the existing core, weighted-sampling, yield, count-plus-shape,
  HepMC3 round-trip, Pythia, DELPHI, PHOTOS, envelope, and performance tests
- distinguish release-blocking correctness gates from scientific benchmarks
  and optional stress diagnostics
- identify stale or contradictory status statements about already implemented
  Pythia, HepMC3, unweighting, and external-tool support
- specify exact commands, optional dependencies, expected outcomes, and
  acceptable failure criteria for a clean release candidate
- keep DELPHI, PHOTOS, analytic-bound, and eta-center conclusions unchanged
- do not tune physics cuts or sampling defaults as part of the inventory
- produce one concise validation document and update the roadmap/status links
- leave implementation of any newly discovered missing gate as a subsequent
  task rather than expanding the inventory without a predeclared contract
