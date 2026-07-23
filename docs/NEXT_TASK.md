# NEXT_TASK

This file should contain exactly one immediate task for Codex or another coding
agent. Broader plans belong in `docs/ROADMAP.md`; project state belongs in
`docs/STATUS.md`.

Measure the effect of true-rapidity rather than pseudorapidity proposal centers
without changing the default sampler.

Requirements:

- add an explicit diagnostic switch that changes only leg-centered proposal or
  refinement centers from pseudorapidity to true rapidity
- use the same generated event and random stream for paired A/B measurements
  wherever possible
- include controlled massive low-velocity fixtures where the two centers
  visibly differ, ordinary relativistic fixtures, hadronic Z events, and
  SoftQCD events
- compare envelope integral, expected and observed acceptance, cell/probe
  count, failure count, and runtime
- verify that target densities and exact Low-current values are identical
  between the two modes
- state separately whether the change helps controlled nonrelativistic cases
  and whether it matters for realistic selected particles
- keep pseudorapidity as the default until the paired result justifies a change
- do not alter the analytic-bound audit conclusion, reduced-probe policy,
  physics density, event yield, multiplicity, or output semantics
