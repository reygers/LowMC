# NEXT_TASK

This file should contain exactly one immediate task for Codex or another coding
agent. Broader plans belong in `docs/ROADMAP.md`; project state belongs in
`docs/STATUS.md`.

Validate Stefan's analytic global and cell-local Low-current majorants without
changing sampling behavior.

Requirements:

- review `docs/analytic-dipole-envelope-proposal-2026-07-23.md`,
  `src/LowCurrent.cpp`, `src/ShapeUnweighting.cpp`, and the physical proposal
  density convention in `docs/ALGORITHM.md`
- add pointwise tests of
  `p_i.k = kT D_i(eta,phi)` and `S(k) = g(eta,phi)/kT^2`
- formulate the rejection-weight bound for `g/q_eta_phi`, including all
  normalization constants, rather than testing a bound on `g` alone
- evaluate `M_T-p_T` stably as `m^2/(M_T+p_T)` and define explicit diagnostic
  fallback behavior for massless, off-shell, or numerically non-timelike legs
- implement global and rectangular cell-local analytic bounds as an
  audit-only diagnostic beside the existing probed envelope
- compare each analytic bound with dense observed maxima on existing toy
  fixtures and held-out Pythia events
- report bound violations, analytic/current envelope ratios, implied
  acceptance, construction cost, and fallback counts
- test true rapidity versus pseudorapidity proposal centers separately; do not
  change the current default as part of the bound audit
- defer the cosh-cos inverse-CDF proposal until the bound study establishes
  whether an analytic path is practically useful
- keep the full `11 x 11` envelope as the default and keep `5 x 5`, factor
  `1.50` diagnostic-only
- do not change physics densities, multiplicity conventions, output metadata,
  or envelope-failure semantics

The go/no-go decision for an analytic sampler should be based on both proof of
coverage and measured tightness. A formally safe but impractically loose
global bound is not sufficient.
