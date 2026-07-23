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

The paper first writes the complete Low formula for the two incoming beam
leptons and all charged outgoing hadrons with their generated momenta. It then
describes a practical per-jet calculation using generated charged particles in
the jet's forward hemisphere, including the corresponding beam particle, and
projecting the photon onto the reconstructed jet direction. It reports
agreement between the Low and Haissinski formulae at generator level.

This leaves an important implementation detail to expose rather than assume.
The benchmark must calculate both:

```text
A. publication-style per-jet forward-hemisphere construction
B. charge-conserving full event:
   incoming e+ and e- + all final stable charged hadrons
```

Construction A is the direct reproduction reference for Table 4. Construction
B is the primary LowMC result because it is a complete charge-conserving
scattering current. If a per-jet subset is not charge conserving, it must not
be passed through LowMC as though it were a valid standalone current; it
belongs in the explicitly labelled publication-reference calculator.

LowMC's general default of beams plus all final stable charged particles should
also be reported as a third selection. This determines whether stable charged
leptons excluded from B matter.

The comparison between A and B determines whether the paper's per-jet
construction is numerically equivalent to the full event current in the
selected region. The comparison between B and the general default isolates
non-hadronic stable charged particles.

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
- it actively selects the incoming beams plus primary quarks, while the
  published bremsstrahlung calculation uses charged outgoing hadrons;
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
5. evaluate both the publication-style per-jet construction and the
   charge-conserving LowMC system made from beams plus final stable charged
   hadrons;
6. integrate or sample the exact Low density with
   `0.2 < E_gamma < 1.0 GeV`,
   `20 deg <= theta_gamma <= 160 deg`, and `pT_jet < 80 MeV`;
7. assign each photon direction to the nearest accepted jet;
8. fill the ten published `8 MeV` bins and normalize by the total number of
   accepted jets;
9. write both per-seed CSV output and an aggregate comparison with Table 4;
10. record event counts, accepted-jet counts, charged multiplicity, flavour,
    statistical uncertainties, integration diagnostics, and runtime.

The `theta_gamma_jet < 0.4 rad` association condition from the exploratory
source should be enabled for the primary comparison and exposed as a named
switch. The effect of removing it must be reported because the headline
kinematic definition emphasizes `E_gamma` and `pT_jet`, whereas the paper's
angular studies cover the range only up to `400 mrad`.

## Staged Validation

### Stage A: Same-event implementation cross-check

For a small fixed set of Pythia events, compare:

- a direct angular integration following the exploratory source; and
- the new LowMC diagnostic using the identical charged system and cuts.

The integrated yield and all populated bins should agree within combined
numerical uncertainties. This isolates implementation and Jacobian mistakes
before using the publication as an external target. Separately compare the
publication-style per-jet construction with the charge-conserving full-event
current; do not include their physics difference in the numerical
implementation-agreement test.

### Stage B: Published bremsstrahlung benchmark

Run at least three independent seeds and compare both construction A and
construction B with Table 4. The predeclared first-pass gates for claiming
that LowMC reproduces the published bremsstrahlung baseline are:

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
