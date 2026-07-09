# NEXT_TASK

This file should contain exactly one immediate task for Codex or another coding
agent. Broader plans belong in `docs/ROADMAP.md`; project state belongs in
`docs/STATUS.md`.

Prototype a validated spatial index for shape-unweighting cell lookup and
envelope-cell sampling.

Requirements:

- review `docs/ADAPTER_CONTRACT.md`, `docs/ADAPTER_PLAN.md`,
  `docs/DECISIONS.md`, `docs/STATUS.md`, `docs/ALGORITHM.md`, and
  `docs/USAGE.md`
- focus on the shape-unweighting hot path identified by stage timing, not on
  changing physics, envelope safety factors, or yield-integration accuracy
- identify the current linear cell operations in `src/ShapeUnweighting.cpp`,
  especially sampled-point cell lookup and envelope-cell sampling
- implement the smallest internal index/data structure that can replace those
  linear scans without changing the public API or sampled distribution
- add focused tests proving the indexed path agrees with the existing linear
  behavior for representative refined grids and boundary points
- keep a linear/reference path or helper available for tests until confidence
  is high
- rerun focused shape-unweighting and HepMC3 afterburner checks after the
  prototype
- update `docs/STATUS.md`, `docs/ALGORITHM.md`, and `docs/NEXT_TASK.md`
  afterward

Open questions:

- Is a tree over the existing refined cells enough, or is there a simpler
  regular-bin acceleration compatible with the current split pattern?
- Can envelope-cell sampling use a prefix-sum/binary-search table without
  changing candidate probabilities?
- Which exact validation run is sufficient to claim the prototype is a
  bookkeeping speedup rather than a sampling change?
