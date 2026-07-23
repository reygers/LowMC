# Analytic Dipole-Channel Proposal For Shape Unweighting

## Context

This note starts from a question about whether the current shape-unweighting
sampler (`src/ShapeUnweighting.cpp`) is the standard way to unweight the
Low-photon `eta`-`phi` distribution, or whether something faster/less
computation-intensive is possible.

The current sampler treats the target density as a numerical black box: it
builds a grid-adaptive proposal, probes each cell for a max weight, runs
calibration trials to catch cases the probe grid misses, and audits severe
cells where a full probe lattice disagrees with the reduced one. This is a
recognized generic technique (comparable to FOAM-style or multichannel
cell-based unweighting), but it is generic precisely because it does not use
the known closed form of the Low current.

This note works out that closed form and what it would buy us. It is a
proposal for discussion, not a plan to implement yet.

## The key algebraic fact

With the photon parameterized as in `photonMomentumFromKTEtaPhi` (kT, eta,
phi) and the mostly-minus metric in `FourVector::dot`, for any charged
particle `i`:

$$
p_i \cdot k = k_T \, D_i(\eta,\phi)
$$

Expanding `p_i . k` and using the transverse-mass decomposition
`E_i = M_i^T cosh(y_i)`, `p_i^z = M_i^T sinh(y_i)` (true rapidity
`y_i = atanh(p_i^z/E_i)`, transverse mass `M_i^T = sqrt(m_i^2 + (p_i^T)^2)`),
this reduces to a clean closed form:

$$
D_i(\eta,\phi) = M_i^T \cosh(\eta - y_i) - p_i^T \cos(\phi - \phi_i)
$$

`intensityPairSum` in `src/LowCurrent.cpp` is `sum_ij c_ij / ((p_i.k)(p_j.k))`
with `c_ij` constants (masses, momenta, charges only). Since each factor in
the denominator contributes one power of `kT`, the full intensity factorizes
exactly:

$$
S(k) = g(\eta,\phi) / k_T^2, \qquad g(\eta,\phi) = \sum_{ij} c_{ij} / (D_i D_j)
$$

This is exact, not an approximation. It also explains something already in
the code: `logUniformKTDensityPhysicalMeasure` (`include/lowmc/KTProposal.h`)
is already `∝ 1/kT^2`, so the existing `kT` sampling is already the exact
answer for that dimension, not a heuristic.

## A free, exact bound

Because `cosh(x) >= 1` and `cos(x) <= 1` for all real `x`:

$$
D_i(\eta,\phi) \;\ge\; M_i^T - p_i^T \;=:\; D_{i,\min} \;>\; 0
$$

for every `(eta, phi)`, not just near particle `i`'s own direction (equality
only in the massless limit). This is the closed-form statement of "mass
regulates the collinear singularity," and it costs nothing to compute — it's
two numbers per particle from the input four-momenta and masses.

From this, every term in `g` is boundable in closed form:

$$
|c_{ij}| / (D_i D_j) \;\le\; |c_{ij}| / (D_{i,\min} D_{j,\min})
$$

Summing over all `i,j` gives a single global number

$$
G_{\max} = \sum_{ij} |c_{ij}| / (D_{i,\min} D_{j,\min})
$$

which is a **provably valid majorant** for `g(eta,phi)` everywhere, computed
once per event in `O(N^2)` — the same cost as one `intensityPairSum`
evaluation. Tighter, still-closed-form bounds are available per grid cell by
restricting the `cosh`/`cos` minimization to that cell's `eta`/`phi` range,
which mirrors the existing grid-refinement structure but replaces "probe grid
points, take the max" with "evaluate a closed-form minimum."

## What this would replace

If this holds up, it removes the reason for most of the machinery in
`ShapeUnweighting.cpp`:

- no probe lattice (`probeSingleCellMaxWeight`, `probeBins`, `baseProbeBins`)
- no calibration trials (`calibrationCellMaxWeights`) to catch peaks the probe
  lattice missed between grid points
- no severe-cell audit (`auditSevereRatioThreshold`,
  `auditSevereCellReportMaxRows`, the reduced/full-probe comparison)
- no risk of an *underestimated* envelope from a finite probe grid, since the
  bound is proven rather than sampled

The `G_max` route trades some acceptance efficiency (a global or coarse-cell
bound is less tight than an empirically tuned per-cell one) for removing the
entire class of failure the audit machinery currently exists to catch.

## A closed-form proposal shape, as a secondary benefit

The diagonal term `-charge_i^2 m_i^2 / D_i(eta,phi)^2` is exactly a
"cosh-cos" peak centered at `(y_i, phi_i)`. This is the true shape that
`SmoothEtaPhiProposal`'s triangular kernel (`include/lowmc/EtaPhiProposal.h`)
is currently approximating by hand. The cosh-cos form has a standard closed
integral (tan-half-angle substitution) in both `eta` and `phi | eta`, so it
can be sampled by exact inverse-CDF instead of the triangular-kernel
heuristic, then combined with the other legs by ordinary multichannel
importance-sampling mixing. Interference (`i != j`) cross terms don't reduce
to a clean single-leg form, but they are already covered by the `G_max` bound
above, and are cheap to evaluate exactly at each sampled point for the
importance weight.

## A concrete finding along the way

Deriving the above surfaced a mismatch worth flagging on its own, independent
of whether the rest of this proposal goes anywhere: the correct center of the
diagonal singularity is the particle's true **rapidity**
`y_i = atanh(p_z/E)`, but `particlePseudorapidity()`
(`include/lowmc/EtaPhiProposal.h`) and `chargedParticleEtaPhi()`
(`src/ShapeUnweighting.cpp`) both compute **pseudorapidity**
`asinh(p_z/p_T)`. These agree only in the massless limit. For relativistic
pions/muons the difference is likely small, but for any charged leg with
non-negligible mass relative to its momentum, the existing adaptive kernels
are centered slightly off from where the true collinear enhancement actually
peaks.

## Effort and risk

This is a real rewrite of the proposal/envelope layer, not a local tweak:

- deriving and validating the cosh-cos inverse-CDF (both the `eta` marginal
  and `phi | eta` conditional) is bounded, standard calculus work, but needs
  to be done and checked carefully (branch cuts in the tan-half-angle
  substitution, etc.)
- `G_max` and any per-cell analytic bounds should be validated against the
  existing numeric probe grid's observed maxima on current test fixtures as
  a regression gate before trusting them (`G_max` should always be >= what
  the probe grid has ever measured)
- a new analytic-channel proposal class would need to be built and A/B
  compared against the current grid-based path using the existing
  `scripts/validate_shape_unweighting*.py` scripts for both correctness
  (shape) and efficiency

Given shape unweighting is currently diagnostic-only
(`docs/DECISIONS.md`), this changes the core sampling design rather than
tuning an existing knob, which is why it's coming to you before any code is
written.

## Questions for Klaus / Codex

1. Does the closed-form factorization and the `D_i,min` bound above look
   right, or is there a subtlety in the current's definition or the
   phase-space measure convention that breaks it?
2. Is the pseudorapidity/rapidity mismatch in the adaptive-kernel centers
   worth fixing independently of the larger proposal, given it's a small,
   self-contained change?
3. Is replacing the probe/calibrate/audit machinery with a provably valid
   (but possibly less tight) analytic bound a direction worth prototyping, or
   is the current empirical approach preferred because it's already
   validated and the efficiency loss from a looser global bound is a bigger
   concern than the removed audit complexity?
4. If this is worth prototyping, should it start as a from-scratch
   diagnostic (parallel to `lowmc_shape_unweighting_diagnostic`) rather than
   touching the current production-track code path at all?

## Proposed next step

Hold off on implementation. Circulate this note for review of the derivation
first, since the main risk is a sign or convention error somewhere in the
factorization, not the overall strategy.
