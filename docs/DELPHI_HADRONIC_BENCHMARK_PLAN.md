# DELPHI Hadronic Soft-Photon Benchmark

## Purpose

The key validation target is the calculated inner-bremsstrahlung spectrum in
DELPHI hadronic `Z` decays, not the anomalous excess itself. LowMC implements
the same leading-soft current used for that calculation, so it should reproduce
the published bremsstrahlung reference when the charged system, event
selection, jet definition, photon cuts, and normalization are aligned.

Only after this baseline works should a LowMC curve be shown beside the
efficiency-corrected `RD - MC` excess. Agreement with the bremsstrahlung curve
would validate LowMC; it would not explain the excess.

Primary reference:

- DELPHI Collaboration, *Evidence for an Excess of Soft Photons in Hadronic
  Decays of Z0*, Eur. Phys. J. C 47 (2006) 273-294,
  [arXiv:hep-ex/0604038](https://arxiv.org/abs/hep-ex/0604038),
  [DOI:10.1140/epjc/s2006-02568-8](https://doi.org/10.1140/epjc/s2006-02568-8)

The machine-readable Table 4 values are stored in
`fixtures/reference/delphi_hadronic_bremsstrahlung_table4.csv`.
The Figure 4b angular `Brems` markers are not tabulated in the paper. Their
vector-digitized values are stored separately in
`fixtures/reference/delphi_hadronic_bremsstrahlung_figure4b_theta_digitized.csv`
and must be labelled as digitized whenever they are shown.

## Published Target

The efficiency-corrected generator-level reference uses:

- `e+ e- -> Z0 -> hadrons` at the `Z` pole;
- photons with `0.2 < E_gamma < 1.0 GeV`;
- photon polar angle to the beam `20 deg <= theta_gamma <= 160 deg`;
- `pT_jet < 80 MeV`, where `pT_jet` is relative to the closest reconstructed
  jet;
- photon-to-parent-jet angle below `400 mrad`, as used by the analysis code
  and reflected in the paper's `0-400 mrad` angular range;
- ten `8 MeV`-wide `pT_jet` bins from `0` to `80 MeV`;
- normalization per accepted jet.

The published efficiency-corrected inner-bremsstrahlung integral is

```text
17.10e-3 photons / accepted jet
```

in this range. Table 4 provides the ten bin integrals. Their sum is
`17.079e-3 photons / accepted jet`, consistent with the quoted rounded total.
The published systematic components for the bremsstrahlung calculation are
`4%` from the soft formula, `5%` from the event generator, and `3%` from the
jet finder.

The much smaller `0.340e-3 photons / jet` number is the prediction after
convolution with DELPHI's conversion and reconstruction efficiency. LowMC
should first compare with the efficiency-corrected `17.10e-3` value; detector
efficiency convolution is a separate task.

Figure 4b provides a complementary angular projection in `25` bins of
`16 mrad` from `theta_gamma = 0` to `0.4 rad`. Its digitized `Brems` markers
sum to `31.663e-3 photons/jet`. Although the paper defines the quoted signal
rate by integrating `pT < 80 MeV`, Figure 4 is explicitly the set of
general-selection distributions and its pT panels extend to `0.4 GeV`.
The angular markers cannot be the `pT < 80 MeV` projection: that projection
must have the same integral as Table 4, about `17.1e-3 photons/jet`.

The paired diagnostic therefore records a joint histogram with `8 MeV` pT
bins out to `0.4 GeV` and `16 mrad` angular bins. Table 4 is the projection
of its first ten pT bins; Figure 4b is compared with the projection over all
fifty displayed pT bins. The summarizer also records the restricted
`pT < 80 MeV` angular projection as a consistency diagnostic, but does not
compare it with Figure 4b. It verifies that this restricted projection
reproduces the primary Table 4-phase-space total.

The angular plot has no tabulated bin uncertainties or covariance matrix.
It is consequently a shape and normalization cross-check, not an additional
chi-square validation gate.

## Published Event and Jet Definition

The DELPHI event selection requires:

- charged multiplicity `Nch >= 5`;
- visible energy `Evis >= 0.2 Ecm`;
- thrust-axis polar angle `30 deg <= theta_thrust <= 150 deg`.

Jets use LUCLUS with `djoin = 3 GeV`. The paper reconstructs them from selected
charged and neutral detector objects. Accepted jets satisfy:

- `30 deg <= theta_jet <= 150 deg`;
- `|p_jet| >= 5 GeV`;
- no identified electron or positron in the jet;
- an additional electron rejection when the jet charged multiplicity is one.

The exact detector track-quality and particle-identification selections cannot
be reproduced at generator level. The first LowMC benchmark should therefore
be labelled a generator-level approximation and should record the effect of
the truth-level electron veto separately.

## Charged System

The primary LowMC comparison uses its standard complete, charge-conserving
event current:

```text
incoming e+ and e- + all final stable charged particles
```

The paper writes the corresponding Low formula for the incoming beams and
charged outgoing particles with generated momenta. Its discussion of
forward-hemisphere particles is useful historical context, but a quark-current
or non-charge-conserving per-jet current is not part of this first benchmark.
Alternative charged-system definitions should be studied only if the standard
full-event comparison exposes a specific discrepancy.

## What the Existing Analysis Code Provides

An earlier exploratory DELPHI analysis source outside this repository provided
useful implementation details:

- Pythia's `ClusterJet("Lund")` with `analyze(event, 0., 3.)`, corresponding
  to a `3 GeV` joining scale;
- the thrust, jet-angle, jet-momentum, photon-energy, photon-angle, and
  closest-jet definitions;
- per-accepted-jet normalization;
- a direct angular integration of the coherent Low current;
- an optional `0.4 rad` photon-to-jet association cut.

It is not itself the published-rate reference in its current active
configuration:

- it enables only `Z -> u ubar`, rather than all five light and heavy quark
  flavours;
- it actively selects the incoming beams plus primary quarks, which is not
  relevant to this first full-event LowMC comparison;
- it generates only ten events;
- it omits the published visible-energy cut and detector-level object
  selections;
- its histogram uses `16 MeV` bins up to `400 MeV`, not the published ten
  `8 MeV` bins below `80 MeV`;
- the plotting macro multiplies the result by an arbitrary factor of `125`,
  so that plot can compare shapes but cannot validate normalization.

An older source version does use incoming beams plus all final charged
particles. That confirms the intended event-wide Low-current implementation
but still does not replace the published Table 4 target.

## First Implementation Contract

Add an opt-in Pythia diagnostic outside `lowmc_core`. It should:

1. generate all `Z -> u,d,s,c,b` modes at the Pythia `Z` mass;
2. disable generator QED photon radiation for the clean baseline and record
   that setting in the output;
3. reconstruct truth jets with Pythia `ClusterJet("Lund")`,
   `yScale = 0`, and `pTscale = 3 GeV`;
4. apply the generator-level event and accepted-jet cuts above;
5. construct the standard charge-conserving LowMC system from the beams plus
   all final stable charged particles;
6. integrate or sample the exact Low density with
   `0.2 < E_gamma < 1.0 GeV`,
   `20 deg <= theta_gamma <= 160 deg`, and `pT_jet < 80 MeV`;
7. assign each photon direction to the nearest accepted jet;
8. fill the ten published `8 MeV` bins and normalize by the total number of
   accepted jets;
9. write both per-seed CSV output and an aggregate comparison with Table 4;
10. record event counts, accepted-jet counts, charged multiplicity, flavour,
    statistical uncertainties, integration diagnostics, and runtime.

The `theta_gamma_jet < 0.4 rad` association condition is mandatory for the
Figure 4d comparison. It is not stated clearly in the publication's headline
selection, but the earlier reproduction study found it necessary and it
matches the upper range of the paper's `theta_gamma` projection. Record this
as an inferred two-dimensional-histogram projection limit, not as an explicit
published cut.

## Staged Validation

### Stage A: Same-event implementation cross-check

For a small fixed set of Pythia events, compare:

- direct numerical integration of the Low density into each complete
  `8 MeV` bin; and
- standard weighted LowMC photon sampling using the identical charged system
  and cuts.

The integrated yield and all populated bins should agree within combined
numerical uncertainties. This isolates implementation and Jacobian mistakes
before using the publication as an external target. The direct calculation is
a validation oracle; weighted photon sampling remains the main LowMC method.

This stage is implemented by
`lowmc_pythia8_delphi_hadronic_benchmark`. The executable writes one row per
accepted event and `pT` bin, including both estimates, the photon-sampling
uncertainty, the number of accepted jets, grid diagnostics, and timings. It
also writes a sibling `_joint.csv` containing the event-level
`pT`-versus-`theta_gamma` histogram. It does not emit photons or change any
LowMC default.

The fixed-seed pilot with `20` accepted events and `47` accepted jets gives:

```text
direct Low integral       0.01915796 +/- 0.00159 photons / jet
weighted LowMC sampling   0.01919287 photons / jet
LowMC / direct            1.001822
LowMC integration error   0.0000286 photons / jet
coarse/fine direct shift  1.53%
DELPHI Table 4 sum        0.017079 photons / jet
```

The uncertainty on the direct result is an event-level delete-one jackknife;
it is much larger than the photon-integration uncertainty. Thus the pilot's
`12.2%` higher central value than DELPHI is only about `1.3` event-statistical
standard errors and is not evidence for a discrepancy. A descriptive
bin-by-bin comparison gives diagonal `chi2/ndf = 11.5/10`; correlations are
not available, so this is not yet the Stage B validation statistic.

The paired result clears the intended Stage A implementation check. The
direct and sampled calculations use the same complete event current and agree
without an arbitrary normalization factor. The first direct-integration
coarse/fine comparison is weakest in the second bin (`4.25%`), but the fine
direct value agrees with the independent LowMC estimate there to `0.42%`.

For the angular projection, the same `20`-event pilot gives:

```text
direct Low integral       0.0355813 +/- 0.00266 photons / jet
weighted LowMC sampling   0.0358333 photons / jet
LowMC / direct            1.00708
DELPHI Figure 4b sum      0.0316626 photons / jet (digitized)
direct / DELPHI           1.12376
normalized-shape L1       0.0313
```

The angular normalization offset is therefore essentially the same as the
Table 4 offset, while the normalized shapes differ by only `3.1%` in L1
distance. The last two angular bins have large event fluctuations in this
small sample and should not be interpreted separately.

As an internal cross-check, the restricted `pT < 80 MeV` angular projection
integrates to `0.01915796 photons/jet`, exactly reproducing the primary pT
total by construction. It is not compared with Figure 4b; doing so would
give a misleading normalization ratio because Figure 4b spans the displayed
`0-0.4 GeV` pT range.

Because `16 mrad` bins are more demanding for deterministic quadrature, a
separate three-event convergence run used `64` log-energy bins, a `40 x 64`
base angular grid, four refinement levels, four sub-bins, and one million
weighted samples per event. It gives `LowMC/direct = 0.99892` for the full
angular projection and `0.99774` for the Table 4 projection. This clears the
angular implementation check without making the expensive high-resolution
quadrature the production estimator.

### Stage B: Published bremsstrahlung benchmark

Run at least three independent seeds and compare the standard LowMC current
with Table 4. The predeclared first-pass gates for claiming that LowMC
reproduces the published bremsstrahlung baseline are:

- the sum of the ten bins agrees with `17.10e-3 photons / jet` within `10%`;
- the binwise comparison has `chi2/ndf < 2` when the published
  bremsstrahlung systematic uncertainties and LowMC numerical uncertainties
  are combined;
- changing seeds does not move the integral by more than its estimated
  statistical uncertainty;
- LUCLUS-to-Durham and LUCLUS-to-JADE variations remain at the few-percent
  level reported by DELPHI.

These are validation gates, not a fit. Do not rescale the LowMC prediction to
the DELPHI curve.

The nominal production block was predeclared as three Pythia seeds
`24680`, `24681`, and `24682`, with `500` accepted events per seed. Use the
pilot's fine direct settings (`48` log-energy bins, `32 x 48` base angular
grid, three refinement levels, and two sub-bins) and `100000` independent
weighted LowMC draws per event. The pilot indicates that this reduces the
event-statistical uncertainty on the combined total to roughly one percent
while keeping the independent sampling uncertainty subdominant.

This nominal block is now complete. It contains `1500` accepted events and
`3794` accepted jets. The combined Table 4-phase-space result is:

```text
direct Low integral       0.01762249 +/- 0.000226 photons / jet
weighted LowMC sampling   0.01768899 +/- 0.000227 photons / jet
LowMC / direct            1.003774
LowMC integration error   0.00000526 photons / jet
coarse/fine direct shift  0.815%
DELPHI Table 4 sum        0.017079 photons / jet
direct / DELPHI           1.031822
diagonal chi2/ndf         13.702/10
```

Thus the total-rate and diagonal-chi-square gates pass without rescaling.
The three direct seed results are:

```text
seed 24680   0.01765605 +/- 0.00038587 photons / jet   1264 jets
seed 24681   0.01826387 +/- 0.00040736 photons / jet   1244 jets
seed 24682   0.01696909 +/- 0.00038134 photons / jet   1286 jets
```

The largest pairwise difference, between seeds `24681` and `24682`, is
`2.32` combined standard errors. A global constant-rate consistency check
gives `chi2 = 5.41` for two degrees of freedom (`p = 0.067`). This is not
compelling evidence for seed dependence, but the literal predeclared
one-standard-error stability wording is not cleanly satisfied. Record the
gate as inconclusive rather than silently weakening it after inspecting the
result.

For the new angular projection, weighted LowMC is the production estimator.
Retain the existing medium-resolution direct calculation in the nominal
files, and additionally run the high-resolution settings above on the first
ten accepted events of each seed as a bounded numerical cross-check. Do not
increase the full production cost by applying that quadrature to all `1500`
events.

The bounded high-resolution block is also complete. Across its `30` events
and `69` jets, `LowMC/direct = 0.99949` in the Table 4 phase space and
`0.99982` for the full Figure 4b angular projection. Seed by seed, the
high-resolution direct total is only `0.28-0.34%` above the medium-resolution
calculation on the same events.

For Figure 4b, the nominal full-range angular projection gives:

```text
direct Low integral       0.03058345 +/- 0.000345 photons / jet
weighted LowMC sampling   0.03072840 photons / jet
LowMC / direct            1.004739
DELPHI Figure 4b sum      0.03166265 photons / jet (digitized)
direct / DELPHI           0.965916
normalized-shape L1       0.042899
```

The integrated normalization is close, but the direct/DELPHI ratio changes
systematically from about `1.24` at small angle to about `0.80` near
`0.4 rad`. With no published bin errors or covariance, this is a diagnostic
shape difference rather than a formal failed gate. It should be tracked when
testing the bounded jet-definition variations.

After the nominal block, repeat bounded truth-level electron-veto and jet
algorithm variations. Report event and accepted-jet statistics separately
for every configuration. Do not tune the photon-to-jet angle cut or rescale
the prediction after inspecting the DELPHI comparison.

At this point in the study, Stage B was not yet closed: the
LUCLUS-to-Durham/JADE gate remained to be run, and the seed-stability result
needed a predeclared independent extension rather than an unplanned choice
made from the observed seed values.

The independent extension is fixed before running as follows:

- use new Pythia seeds `24683`, `24684`, and `24685`;
- generate `500` accepted events per seed with exactly the nominal LUCLUS,
  direct-integration, and weighted-sampling settings;
- require the aggregate rate of the new block to agree with the original
  three-seed block within two combined event-level jackknife standard errors;
- evaluate a constant-rate consistency statistic across all six seed
  estimates and require `p > 0.05`.

These criteria supersede neither the original recorded seed values nor the
other Stage B gates. They provide a statistically explicit resolution of the
otherwise ambiguous literal one-standard-error wording.

The bounded jet/electron variations are also fixed before running:

- use the original seeds `24680`, `24681`, and `24682`;
- process `200` accepted events per seed and retain the nominal photon
  integration and sampling settings;
- compare LUCLUS (`djoin = 3 GeV`) with Durham and JADE at `ycut = 0.01`,
  matching the settings in the existing DELPHI reproduction example;
- repeat LUCLUS with the truth-level electron veto disabled;
- compare unscaled per-accepted-jet totals, accepted-jet counts, pT shapes,
  and the Figure 4b angular diagnostic.

The variation block is diagnostic: no jet cut, photon cut, or normalization
may be tuned after inspecting its outcome.

Both predeclared follow-ups are now complete. The independent extension has
`1500` accepted events and `3729` accepted jets. Its direct rate is
`0.0176347 +/- 0.000222 photons/jet`, only `0.039` combined jackknife standard
errors from the original block. Across all six seeds, the constant-rate test
gives `chi2 = 8.38` for five degrees of freedom (`p = 0.136`). The extension
therefore resolves the seed-stability gate as passed.

Combining all six nominal seeds gives:

```text
accepted events / jets    3000 / 7523
direct Low integral       0.01762855 +/- 0.000158 photons / jet
weighted LowMC sampling   0.01768964 photons / jet
LowMC / direct            1.003465
direct / DELPHI           1.032177
diagonal chi2/ndf         12.082/10
```

The six-seed Figure 4b projection is also stable: its direct integral is
`0.0306705 +/- 0.000244 photons/jet`, LowMC/direct is `1.00381`,
direct/DELPHI is `0.96867`, and the normalized-shape L1 distance is `0.04139`.

The bounded `600`-event variations give:

| configuration | accepted jets | rate estimator | photons/jet | ratio to LUCLUS |
|---|---:|---|---:|---:|
| LUCLUS, electron veto | 1525 | weighted LowMC | 0.0174565 | 1.0000 |
| Durham, `ycut = 0.01` | 1265 | weighted LowMC | 0.0193359 | 1.1077 |
| JADE, `ycut = 0.01` | 1415 | weighted LowMC | 0.0179980 | 1.0310 |
| LUCLUS, no electron veto | 1659 | weighted LowMC | 0.0192899 | 1.1050 |

For the two jet-algorithm comparisons, the direct ratios to the bounded
LUCLUS baseline are `1.1082` for Durham and `1.0322` for JADE. JADE therefore
meets the intended few-percent check, while Durham does not. The normalized
shape changes are much smaller than the rate changes: the pT L1 distances
from LUCLUS are `0.00971` and `0.00553`, and the angular L1 distances are
`0.01019` and `0.00675`, for Durham and JADE respectively. The principal
difference is thus the accepted-jet count and per-jet normalization, not a
large deformation of the photon spectra.

The no-electron-veto row must use weighted LowMC as its rate estimator.
The deterministic angular grid is not sufficiently resolved for the narrow
electron-collinear structures: its aggregate coarse/fine shift is `10.3%`,
and one seed shifts by `27.2%`. Weighted LowMC shows that removing the veto
raises the rate by `10.5%`; this is a useful truth-selection diagnostic, not
evidence of a LowMC/direct inconsistency.

The total-rate, spectral, implementation, and seed-stability gates pass. The
Durham variation remains a literal non-pass of the predeclared few-percent
expectation; no `ycut` or photon cut was tuned in response.

A final accounting-only run applies all three algorithms to exactly the same
`600` LUCLUS-selected Pythia events, without repeating a Low-current
calculation:

| algorithm | accepted events | all jets | accepted jets | accepted jets / LUCLUS |
|---|---:|---:|---:|---:|
| LUCLUS-like Lund | 600 | 1821 | 1525 | 1.0000 |
| Durham | 595 | 1458 | 1255 | 0.8230 |
| JADE | 598 | 1813 | 1412 | 0.9259 |

Durham has fewer accepted jets than LUCLUS in `233` events, the same count in
`366`, and more in only one. Its median closest-axis displacement from a
LUCLUS accepted jet is only `0.033 rad` (`0.121 rad` at the 90th percentile),
so the main effect is a lower jet multiplicity rather than uniformly unrelated
jet axes. In the independent photon samples, Durham reduces the total direct
photon numerator to `0.9193` of LUCLUS but reduces the accepted-jet
denominator further, to `0.8295`; their ratio is the observed
`1.1082` per-jet shift. JADE gives the analogous ratios `0.9578/0.9279 =
1.0322`.

This resolves the source of the Durham result. Pythia's truth-level Durham
configuration produces substantially fewer jets at the fixed `ycut`, whereas
the few-percent statement in DELPHI refers to an algorithm variation inside
its reconstructed-jet analysis. The present truth-level switch is therefore
not a like-for-like reproduction of that detector-level systematic. Record it
as a documented limitation of the proxy, not as a LowMC physics discrepancy.

The DELPHI benchmark is now complete. It establishes the intended result:
unscaled LowMC reproduces the published inner-bremsstrahlung pT spectrum and
normalization. The final six-seed figure includes the DELPHI curve, direct Low
integral, weighted LowMC result, and binwise ratios. No further DELPHI
production run is planned.

### Stage C: Interpret the excess

For interpretation and presentation, one may overlay:

- the efficiency-corrected DELPHI inner-bremsstrahlung prediction;
- the validated LowMC prediction; and
- the efficiency-corrected `RD - MC` excess.

The expected outcome is that LowMC reproduces the bremsstrahlung baseline,
which lies well below the excess. That establishes the Standard Model
leading-soft reference used in the soft-photon puzzle; it does not generate or
explain the anomalous component.

The validated unscaled bremsstrahlung comparison is complete. A future overlay
with the efficiency-corrected `RD - MC` points is a presentation extension,
not another LowMC validation gate.

## Scope Boundary

This benchmark is independent of the reduced-probe and analytic-envelope
work. It does not change production envelope defaults, particle-selection
defaults, or failure policies.
