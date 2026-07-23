# NEXT_TASK

This file should contain exactly one immediate task for Codex or another coding
agent. Broader plans belong in `docs/ROADMAP.md`; project state belongs in
`docs/STATUS.md`.

Implement a single release-validation driver for the first public LowMC
research-prototype release.

Requirements:

- follow `docs/RELEASE_VALIDATION_MATRIX.md` exactly; do not redefine its
  physics criteria
- add one script that creates fresh Release build directories and runs `R1`,
  `R2`, and `R3`
- provide explicit switches and dependency-prefix arguments for the applicable
  `C1` combined HepMC3/Pythia8 gate and the optional PHOTOS build
- stop on a blocking failure, preserve the failing logs, and return nonzero
- write a machine-readable manifest containing the LowMC revision, dirty-tree
  state, compiler/CMake versions, enabled components, dependency paths or
  reported versions, commands, exit statuses, durations, and discovered test
  counts
- write a concise Markdown summary suitable for attaching to a release
- test command construction and failure propagation without requiring all
  optional dependencies in the default CTest path
- do not run DELPHI or PHOTOS production-statistics studies automatically
- do not tune physics cuts, envelopes, or sampling defaults
