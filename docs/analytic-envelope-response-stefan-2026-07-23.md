# Response to the Analytic Dipole-Envelope Proposal

Prepared for discussion with Stefan, 23 July 2026.

## Executive summary

We followed up the analytic proposal in
`analytic-dipole-envelope-proposal-2026-07-23.md` with two separate studies:

1. an audit of the exact factorization and the proposed global and cell-local
   analytic bounds; and
2. a paired study of true-rapidity versus pseudorapidity centers in the
   current shape-unweighting algorithm.

The main conclusions are:

- The factorization
  `p_i.k = kT D_i(eta,phi)` and `S(k) = g(eta,phi)/kT^2` is correct.
- The absolute-value global and rectangular cell-local constructions are
  genuine bounds for finite massive timelike charged legs when all proposal
  normalization factors are included.
- The bounds were safe in all controlled tests and in a predeclared block of
  thirty held-out Pythia events.
- The global bound is much too loose for sampling. The cell-local analytic
  bound is also substantially looser than the current envelope on realistic
  events.
- True rapidity is the exact center of the massive-leg denominator, but using
  it only for the current charged-direction targeting does not improve the
  sampler.
- We therefore keep the present grid/cell-local sampler and pseudorapidity
  targeting as defaults. The analytic formulas remain useful as independent
  audit and possible future cell-certification tools.
- The proposed cosh-cos inverse-CDF sampler is not the next priority because
  the bound study did not establish a competitive analytic rejection path.

This work also fits into a broader comparison with PHOTOS and VinciaQED.
PHOTOS has already provided both an independent physics validation in
`Z -> mu+ mu-` and useful sampling ideas. VinciaQED is an important future
soft-kernel comparison and helps define LowMC's niche, but no direct
LowMC/Vincia event comparison has yet been performed.

## 1. Answer to the factorization and bound question

For each charged leg we verified point by point that

```text
p_i.k = kT D_i(eta,phi)

D_i = mT_i cosh(eta-y_i) - pT_i cos(phi-phi_i)
```

and that

```text
S(k) = g(eta,phi) / kT^2.
```

The tests compare the factorized expression with the existing four-vector
calculation in `LowCurrent`. They cover asymmetric momenta, incoming and
outgoing signs, finite masses, dense angular grids, and the complete photon
weight.

For numerical stability the denominator is evaluated as

```text
D_i =
    (mT_i-pT_i)
  + 2 mT_i sinh^2((eta-y_i)/2)
  + 2 pT_i sin^2((phi-phi_i)/2),

mT_i-pT_i = m_i^2 / (mT_i+pT_i).
```

This avoids subtracting two nearly equal numbers for highly relativistic
particles.

The bound must be applied to the complete rejection weight, not just to
`g`. For LowMC's log-uniform `kT` proposal it is

```text
w =
  alpha/(2 pi)^2
  x log(kTmax/kTmin)
  x g(eta,phi) / q_eta_phi(eta,phi).
```

The implementation includes all these factors and uses the local
`q_eta_phi` density in each proposal cell.

The analytic certification currently requires every selected charged leg to
be finite, massive, and timelike. Massless, off-shell, numerically
non-timelike, or non-finite inputs produce an explicit unavailable/fallback
diagnostic rather than an asserted finite bound.

One conceptual qualification is important. With LowMC's metric and current
conventions, the diagonal pair terms are negative. Their absolute cosh-cos
shapes can still be useful proposal components, but they are not a positive
decomposition of the physical coherent intensity. Interference cancellations
remain essential. This is one reason the absolute-value sum can be very
conservative.

### Held-out bound result

The predeclared test used Pythia seeds `46801`, `46802`, and `46803`, five
events per seed, for both hadronic Z and 13 TeV SoftQCD events. Every event
used `1000` calibration trials and `10000` independent rejection trials.

| Process | Events | Invalid legs | Fallback cells | Probe/calibration violations | Trial violations | Maximum observed/bound | Maximum trial/bound |
|---|---:|---:|---:|---:|---:|---:|---:|
| hadronic Z | 15 | 0 | 0 | 0 | 0 | 0.1278 | 0.1212 |
| SoftQCD | 15 | 0 | 0 | 0 | 0 | 0.03021 | 0.03010 |

Thus the tested analytic construction passed the safety gate with a
substantial margin.

It failed the predeclared practical tightness gate:

| Process | Median analytic/current cell envelope | 90th percentile | Maximum | Median current acceptance | Median analytic acceptance |
|---|---:|---:|---:|---:|---:|
| hadronic Z | 8.34 | 13.98 | 86.96 | 0.102 | 0.0123 |
| SoftQCD | 40.52 | 69.39 | 1339.38 | 0.0672 | 0.00156 |

The global-bound/current-envelope ratio is of order `10^9`. The cell-local
bound is much better than the global one and is relatively cheap to construct,
but the acceptance loss makes it uncompetitive as a standalone rejection
sampler for realistic multi-particle events.

Our answer to the first question is therefore: the algebra and proof strategy
are correct, with the finite-mass/timelike and full-proposal-normalization
qualifications above.

## 2. Answer to the rapidity-center question

The observation about the center is also correct. At fixed photon `kT`, the
minimum of the massive-leg denominator is at the charged particle's true
rapidity `y_i`, not its pseudorapidity.

We added a diagnostic-only center convention to the shape-unweighting
configuration and compared paired runs with identical events, seeds, phase
space, calibration, and trial counts. The switch affects only:

- charged-direction targeted refinement;
- targeted base probing;
- boosted-cell identification; and
- charged-near audit labels.

It does not affect `LowCurrent`, the physical density, event yield,
multiplicity, or acceptance formula. Unit tests confirm identical Low
densities on a fixed point set.

For deliberately slow pions with `p = 0.1 GeV` and pseudorapidity `+/-2`,
true rapidity moves the conservative boosted-cell treatment onto the exact
peak. This makes the envelope about `9.9%` larger and reduces expected
acceptance. With explicit local refinement instead of the boosted-cell policy,
the envelope changes by only `0.02%`. For ordinary `10 GeV` pions the two
center conventions are numerically identical in the tested configurations.

The realistic paired study used seeds `57911`, `57912`, and `57913`, again
with fifteen hadronic-Z and fifteen SoftQCD events:

| Process | Median `|eta-y|` | Largest event `|eta-y|` | Median true/pseudo envelope | Maximum ratio | Failures pseudo/true |
|---|---:|---:|---:|---:|---:|
| hadronic Z | 0.00441 | 1.516 | 1.00031 | 1.01215 | 0 / 0 |
| SoftQCD | 0.0973 | 3.815 | 1.01765 | 1.06795 | 0 / 0 |

For SoftQCD, true-rapidity targeting also adds a median of `6174` probe
evaluations and `14` promoted cells.

Our answer to the second question is therefore: true rapidity is
mathematically correct and should be used in exact analytic constructions,
but a center-only change does not improve the current practical algorithm.
Pseudorapidity remains the default. We did not retune the older
local-box/triangular weighted diagnostic proposals because they are not the
production path and the center-only production-track test showed no gain.

## 3. Should the analytic bound replace probe/calibrate/audit?

No, not in its present absolute-value form.

The result separates formal safety from practical usefulness:

- The analytic bound removes finite-probe undercoverage by construction.
- The global version is unusably conservative.
- The rectangular cell-local version is safe but still loses too much
  acceptance on realistic events.
- The current empirical cell-local envelope is not formally proven, but it is
  much tighter and has extensive failure detection and stress validation.

The analytic machinery is nevertheless worth keeping. Plausible future roles
are:

- an independent audit of empirical cell bounds;
- certification of selected cells;
- a fallback for invalid or suspicious empirical cells; or
- an ingredient in a more refined branch-and-bound construction that retains
  more cancellation information.

Any such hybrid should be motivated by a concrete failure mode or performance
need. The present results do not motivate replacing the current sampler.

## 4. Was a parallel prototype the right development path?

Yes. The analytic calculation was added as an opt-in audit beside the
unchanged shape sampler. It never supplies the active sampling bound. Tests
confirm that enabling it does not change the empirical envelope, accepted
candidates, status, or defaults.

This let us answer both the proof and efficiency questions without committing
to a sampler rewrite. Based on the result, we have not implemented the
cosh-cos inverse-CDF proposal described in the note.

## 5. Relation to PHOTOS

PHOTOS overlaps with LowMC most cleanly for radiation in a specified neutral
decay. We therefore built an independent benchmark for a fixed

```text
Z at rest -> mu+ mu-.
```

The PHOTOS 3.64 baseline uses:

- interference enabled;
- exponentiation disabled;
- double and higher bremsstrahlung disabled;
- pair emission disabled;
- the exact Z matrix-element correction disabled; and
- the same fine-structure constant and common photon analysis acceptance as
  LowMC.

This isolates the shared leading-soft current as closely as the installed
fixed-order implementation permits. Primary angular comparisons use the
common pre-radiation Born muons because PHOTOS applies recoil while LowMC
currently does not.

### Physics comparison

In the initial `10000`-decay pilot, PHOTOS gave
`0.0791 +/- 0.0028` accepted photons per decay, while the converged LowMC
integral was `0.0801377`, corresponding to `PHOTOS/LowMC = 0.987`.
An independent PHOTOS seed pair gave `0.995`. Normalized `kT`, `eta`, `phi`,
and nearest-Born-muon angular shapes showed no reproducible discrepancy.

The subsequent three-seed scan used `300000` total decays in each of six
nested photon-energy windows:

| `E_gamma,max` [GeV] | LowMC / PHOTOS |
|---:|---:|
| 0.060 | `0.995 +/- 0.019` |
| 0.075 | `1.010 +/- 0.013` |
| 0.100 | `1.013 +/- 0.010` |
| 0.150 | `1.007 +/- 0.008` |
| 0.250 | `1.009 +/- 0.007` |
| 0.500 | `1.009 +/- 0.007` |

All windows are compatible with unity, and the softest point is within `0.5%`
of unity. This is an important independent validation of LowMC's signs,
normalization, finite-mass treatment, interference, and integration in a
domain where the two models should overlap.

The comparison does not imply that PHOTOS can replace LowMC for a general
full scattering event. PHOTOS traverses and radiates individual decays,
whereas LowMC evaluates one configurable coherent current over a user-selected
set of incoming and outgoing charged particles. PHOTOS is the natural mature
tool for ordinary decay radiation; LowMC's intended use is broader and more
explicitly configurable at leading-soft level.

### Sampling lessons from PHOTOS

PHOTOS motivated a separate analytic leg-centred proposal proportional to

```text
1 / (1 - beta_i cos(theta_i,gamma)).
```

The proposal mixes normalized charged-leg components with a uniform support
floor and always evaluates the exact coherent Low current for the physical
weight.

This proposal is not universally superior to the grid:

- it performs well for boosted and clustered controlled fixtures;
- it loses to the grid for a low-boost HepMC3 validation fixture;
- on realistic Pythia events it usually has lower effective sample fraction
  per draw, but its low construction cost gives competitive or better
  effective samples per second; and
- its residual upper tail is much more stable than rare grid outliers.

An empirical envelope calibrated with `20000` proposal draws and checked on
`100000` independent draws per event had zero overshoots in sixty hadronic-Z
and sixty SoftQCD events. With safety factor `1.2`, the largest trial/bound
ratios were `0.844` and `0.915`, and mean expected acceptances were `7.97%`
and `6.58%`.

This is useful as an opt-in independent reference, but it remains an empirical
majorant. It does not provide the formal guarantee of Stefan's bound, and it
has not displaced the current grid.

The combined lesson is informative:

| Approach | Main strength | Main limitation | Current role |
|---|---|---|---|
| Stefan global/cell bound | formal coverage | very loose on realistic events | audit/certification research |
| PHOTOS-inspired leg proposal | cheap sampling and stable empirical tails | no formal event-independent majorant | opt-in diagnostic reference |
| Current grid/cell envelope | best validated practical balance | empirical probing and calibration machinery | default |

## 6. Relation to VinciaQED

VinciaQED is a closer physics overlap than PHOTOS for a general coherent
charged system. Its multipole mode contains the same scalar soft/eikonal
structure that underlies LowMC's leading current.

It is nevertheless a broader and differently organized model:

- VinciaQED is integrated into Pythia's parton systems and shower history.
- It supplies ordering, competition between emissions, recoil maps, repeated
  emissions, and spin-dependent collinear contributions.
- Its coherent systems are organized by shower systems, scales, and decays.
- LowMC instead applies a transparent leading-soft calculation to a
  configurable set of charged particles in a completed event record.

Consequently, LowMC's default current built from incoming beams plus final
stable charged particles is not automatically equivalent to Vincia's
system-by-system coherence definition. That difference is a physics question,
not merely an implementation detail.

Vincia also suggests a possible pair- or antenna-centred proposal based on

```text
(p_i.p_j) / ((p_i.k)(p_j.k)).
```

We have not implemented this. A naive pair mixture has `O(N^2)` components,
can become expensive for high charged multiplicity, and faces the same
positive-majorant versus interference-cancellation tension exposed by
Stefan's absolute-value bound. The present LowMC work therefore stops at
recording it as a later exploratory option.

A meaningful future Vincia comparison would use the same charge-conserving
particle configuration in both programs and first isolate the one-photon,
very-soft scalar-eikonal limit, with spin-dependent collinear terms disabled
where possible. Only afterwards should one turn on or interpret differences
from charged-system selection, spin terms, recoil, ordering, and repeated
emissions.

No such controlled direct comparison has yet been run, so the Vincia
discussion is presently a positioning and design conclusion rather than a
numerical validation result.

## 7. Overall project implication

The work supports a fairly clear division:

- Use LowMC as a generator-independent, configurable leading-soft reference
  and afterburner with transparent densities and diagnostics.
- Use PHOTOS as the mature reference and natural production tool for
  radiation in decays; the completed `Z -> mu+ mu-` benchmark validates the
  shared soft limit.
- Treat VinciaQED as the principal future comparison for coherent shower
  physics beyond an isolated decay, while recognizing that its system
  definition and shower evolution are not the same as LowMC's event-wide
  current.
- Keep the current grid/cell shape sampler as the LowMC default.
- Preserve Stefan's analytic formulas and the PHOTOS-inspired analytic
  proposal as complementary diagnostics rather than competing production
  replacements.

This means the analytic-envelope study has reached a useful stopping point.
Further work on analytic proposals should be driven by a specific unresolved
failure, a need for independent certification, or a planned Vincia comparison,
not pursued as a sampler rewrite for its own sake.

## Reproducibility and project records

The main implementation revisions are:

- `4cb5956`: analytic global/cell-bound audit;
- `3785b7f`: true-rapidity targeting audit.

Detailed records:

- `docs/ANALYTIC_BOUND_AUDIT_PLAN.md`
- `docs/RAPIDITY_CENTER_AUDIT.md`
- `docs/PHOTOS_BENCHMARK_PLAN.md`
- `docs/LEG_CENTERED_PROPOSAL_PLAN.md`
- `docs/ROADMAP.md`
