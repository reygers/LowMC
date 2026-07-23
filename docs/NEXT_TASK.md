# NEXT_TASK

This file should contain exactly one immediate task for Codex or another coding
agent. Broader plans belong in `docs/ROADMAP.md`; project state belongs in
`docs/STATUS.md`.

Review and, with explicit maintainer approval, tag the first public LowMC
`0.1.0` research-prototype release candidate.

Requirements:

- inspect the retained clean-tree `manifest.json`, `summary.md`, and failing-log
  absence before taking any release action
- review `docs/RELEASE_NOTES_0.1.0.md` and the public repository contents
- confirm the intended tag name, either `v0.1.0-rc1` or `v0.1.0`
- create and publish a tag or hosted release only after explicit maintainer
  approval
- attach or archive the retained validation manifest and summary with the
  release
- do not change physics defaults, add a production HepMC3 alias, or rerun
  DELPHI/PHOTOS production-statistics studies as part of tagging
