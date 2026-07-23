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

The exploratory source at
`/Users/reygers/uni/my_papers/published/Soft photon puzzle/Delphi/src_Low/delphi_soft_photons.cpp`
provides useful implementation details:

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
does not emit photons or change any LowMC default.

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

The next nominal production block is predeclared as three Pythia seeds
`24680`, `24681`, and `24682`, with `500` accepted events per seed. Use the
pilot's fine direct settings (`48` log-energy bins, `32 x 48` base angular
grid, three refinement levels, and two sub-bins) and `100000` independent
weighted LowMC draws per event. The pilot indicates that this reduces the
event-statistical uncertainty on the combined total to roughly one percent
while keeping the independent sampling uncertainty subdominant.

After the nominal block, repeat bounded truth-level electron-veto and jet
algorithm variations. Report event and accepted-jet statistics separately
for every configuration. Do not tune the photon-to-jet angle cut or rescale
the prediction after inspecting the DELPHI comparison.

### Stage C: Interpret the excess

Only after Stages A and B pass, overlay:

- the efficiency-corrected DELPHI inner-bremsstrahlung prediction;
- the validated LowMC prediction; and
- the efficiency-corrected `RD - MC` excess.

The expected outcome is that LowMC reproduces the bremsstrahlung baseline,
which lies well below the excess. That establishes the Standard Model
leading-soft reference used in the soft-photon puzzle; it does not generate or
explain the anomalous component.

## Scope Boundary While Review Is Pending

This plan is independent of the reduced-probe recommendation sent to Stefan.
Preparing and later running this benchmark must not change reduced-probe
settings, production envelope defaults, particle-selection defaults, or
failure policies. The PHOTOS and analytic-envelope studies remain closed
unless a specific discrepancy in this benchmark motivates reopening them.
