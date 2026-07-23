# NEXT_TASK

This file should contain exactly one immediate task for Codex or another coding
agent. Broader plans belong in `docs/ROADMAP.md`; project state belongs in
`docs/STATUS.md`.

Prepare the clean-tree evidence and release notes for the first public LowMC
research-prototype release candidate.

Requirements:

- start from a committed tree and use
  `scripts/validate_release_candidate.py` with the combined HepMC3/Pythia8 and
  optional PHOTOS switches
- require `dirty_tree=false` and all applicable gates to pass in the retained
  manifest
- draft concise release notes for version `0.1.0` that distinguish the stable
  generator-independent core from provisional diagnostic interfaces
- state the no-recoil convention, leading-soft scope, envelope failure policy,
  optional dependencies, and known boosted/high-multiplicity limitations
- summarize the frozen DELPHI and PHOTOS validation conclusions without
  rerunning their production-statistics studies
- include build/install/quickstart commands and links to the release evidence
- check that README, installation, license, version, and public header
  information are mutually consistent
- do not add a production HepMC3 alias, change physics defaults, or expand the
  release scope
