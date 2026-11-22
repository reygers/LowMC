# ALICE soft-dielectron benchmark

This workflow compares the cocktail-subtracted ALICE Figure 4 data with two
calculations:

1. the analytic Balek-Pisutova-Pisut (BPP) dipole term; and
2. the event-by-event virtual-photon prediction from LowMC, using either the
   leading Low current or Lichard Eq. (2.8).

It is a validation benchmark. The selected current must be stated with every
result. Neither option is a complete finite-energy prediction: both use the
soft external-leg framework and omit hadronic recoil.

## Balek reference

The core implements BPP Eq. (9a), including the exact lepton-mass factor from
their Eq. (3):

\[
\frac{d^3\sigma_{ee}}{dy\,dp_T\,dM}
=
\sigma_{\rm hadr}
\frac{4\alpha^2\langle\Delta Q^2\rangle}{3\pi^2}
\frac{p_T}{M(p_T^2+M^2)}
\left(1+\frac{2m_e^2}{M^2}\right)
\sqrt{1-\frac{4m_e^2}{M^2}}.
\]

The plotting script uses the ALICE inputs
`<Delta Q^2> = 1.32` and `sigma_hadr = 57.8 mb`. The latter is the V0 visible
cross section cited by the ALICE analysis. Both values are configurable.

## LowMC configuration

The benchmark executable has two comparison modes. The default
`--comparison detector` mode uses:

- pp collisions at `sqrt(s) = 13 TeV`;
- Pythia8 Monash tune (`Tune:pp = 14`);
- `SoftQCD:nonDiffractive`;
- Pythia QED showers and hadron-level QED radiation disabled;
- incoming protons plus all final stable charged particles in the Low current;
- the leading Low current by default, with Lichard Eq. (2.8) selectable through
  `--current lichard`;
- `0.15 < M < 0.6 GeV`, `0 < pT,ee < 0.8 GeV`;
- both electrons satisfying `|eta_e| < 0.8` and
  `0.075 < pT,e < 10 GeV`.

The mass projection additionally requires `pT,ee < 0.4 GeV`. Results are
divided by the sampled rapidity width `Delta y = 1.6` and normalized with
Pythia's generated non-diffractive cross section. There is no separate cut on
pair rapidity or pair pseudorapidity. The sampler uses
`-0.8 < y_ee < 0.8` only as proposal support: requiring both daughters to
satisfy `|eta_e| < 0.8` already implies that the accepted pair lies in this
rapidity interval.

For comparison with the published, single-electron-acceptance-corrected
Figure 4 spectra, use `--comparison publication`. This mode:

- applies no cuts to the daughter electron and positron;
- uses the angle-integrated Kroll-Wada conversion density;
- explicitly integrates over the provisional midrapidity interval
  `-0.8 < y_ee < 0.8` and reports the result per unit rapidity;
- retains the same pair-level mass and transverse-momentum projection cuts.

The numerical pair-rapidity interval is not stated explicitly with Figure 4
in the publication and should be confirmed with the analysis authors. The
publication mode records its rows as `decay=angle_integrated`; detector mode
continues to write both polarization-aware and isotropic rows.

For every event, the executable calculates the same binned integrals twice:

- ordinary pseudorandom weighted sampling;
- a deterministically seeded, randomly shifted six-dimensional Sobol rule.

Both polarization-aware and isotropic virtual-photon decays are retained. The
CSV contains event-level standard errors and weighted quantiles of the maximum
soft-expansion parameters

\[
\max_i \frac{M^2}{2|p_i\cdot q|},
\qquad
\max_i \frac{q^0}{E_i}.
\]

Large values warn that at least one selected external leg is outside the
strict soft regime. They are diagnostics, not uncertainty estimates.

For an exploratory two-dimensional agreement map, add
`--extended-agreement-map`. This appends `observable=mass_qT` rows over
`0.02 < M < 0.60 GeV` and `0 < pT,ee < 2 GeV`. The rows contain the mass and
pair-transverse-momentum bin edges and densities in `mb / GeV^2 / y`; the
standard one-dimensional projections remain unchanged.

## Run

```sh
cmake -S . -B build-pythia8 -DLOWMC_ENABLE_PYTHIA8=ON
cmake --build build-pythia8 --target lowmc_pythia8_alice_dielectron_benchmark

./build-pythia8/adapters/pythia8/lowmc_pythia8_alice_dielectron_benchmark \
  build-pythia8/alice_dielectron_prediction.csv \
  --events 1000 --samples 10000 --direct-samples 16384 --seed 12345 \
  --current leading --comparison publication

MPLCONFIGDIR=build-pythia8/matplotlib \
python3 scripts/plot_alice_dielectron_comparison.py \
  --prediction build-pythia8/alice_dielectron_prediction.csv \
  --output-prefix build-pythia8/alice_fig4_lowmc
```

The plotting script writes PDF and PNG versions. HEPData statistical errors
are shown as bars and systematic errors as boxes. The separate 5% luminosity
uncertainty is noted on the figure. HEPData bins available only as upper limits
are omitted because their numerical limits are absent from the published YAML
tables.

## Interpretation boundary

The non-diffractive result is not yet the inclusive inelastic prediction.
Single-, double-, and central-diffractive Pythia components should be evaluated
separately before forming an inclusive curve. Repeat the integration with
`--current lichard` to expose the finite-mass current dependence. The
difference is a model comparison, not by itself a complete theory uncertainty.

To identify where the two current prescriptions agree, run matched leading
and Lichard jobs with `--extended-agreement-map`, then plot the joint rows:

```sh
MPLCONFIGDIR=build-pythia8/matplotlib \
python3 scripts/plot_alice_dielectron_agreement_map.py \
  --leading output/alice_dielectron_agreement_leading.csv \
  --lichard output/alice_dielectron_agreement_lichard.csv \
  --output-prefix output/alice_dielectron_agreement_map
```

The first panel shows `100 * |N_Lichard / N_leading - 1|`, with 3%, 5%, and
10% agreement contours. The second panel shows the yield-weighted median of
the strict event-level energy-softness diagnostic `max_i(q0/E_i)`. By default,
the script tests the joint requirements of 5% current agreement and a median
energy fraction below 0.3. Cells whose event-level relative error exceeds 15%
are masked. Agreement of the currents tests the finite-mass prescription only;
it does not establish energy softness.
