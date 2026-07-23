# PHOTOS BENCHMARK PLAN

## Purpose

Define a reproducible comparison between LowMC and PHOTOS for
`Z -> mu+ mu-` before implementation begins.

For a neutral parent and the same two outgoing muons, both programs contain
the same charge-conserving leading-soft current. Their accepted real-photon
predictions should therefore converge as the photon-energy range is restricted
to the strict soft limit. The comparison is intended to validate LowMC's
charge signs, normalization, finite-mass behavior, interference pattern,
sampling, and unweighting against an independent mature afterburner.

This benchmark does not test whether PHOTOS's decay-by-decay event traversal
is equivalent to LowMC's configurable current for a general scattering event
with distinct incoming legs.

## Benchmark Configuration

Use a fixed Born configuration rather than a distribution of `Z` production
events:

```text
Z at rest -> mu+ along +x, mu- along -x
```

The implementation must record the exact numerical values used for:

- `Z` mass
- muon mass
- fine-structure constant
- input four-momenta
- random seeds
- PHOTOS version and build configuration
- LowMC revision and build configuration

The muon momenta must be computed from the chosen masses rather than rounded
independently. A copy of every pre-radiation Born four-vector must be retained
for common analysis.

## PHOTOS Mode

The first comparison should use the PHOTOS universal kernel in the simplest
mode that overlaps with LowMC:

- interference enabled
- single real photon
- exponentiation disabled
- process-specific exact hard-photon matrix-element corrections disabled
- fine-structure constant matched to LowMC
- PHOTOS generation cutoff below the lower edge of the common analysis
  acceptance

The exact PHOTOS API calls and their defaults must be verified against the
installed version and recorded in the benchmark output. If the selected
PHOTOS version cannot disable a feature independently, that limitation must be
reported rather than silently approximated.

Later comparisons may enable multiphoton generation, exponentiation, or exact
matrix-element corrections, but those are comparisons to physics beyond
LowMC's current leading-soft model.

## Common Photon Acceptance

Generate PHOTOS photons with a cutoff safely below the analysis threshold and
apply the same analysis-level acceptance to both samples. Do not assume that a
PHOTOS fractional energy cutoff is identical to a LowMC `k_T` cutoff.

For every accepted photon, compute in the original `Z` rest frame:

- photon energy `E_gamma`
- `k_T`, `eta`, and `phi` relative to the fixed benchmark axes
- angle to the pre-radiation `mu+`
- angle to the pre-radiation `mu-`
- angle to the nearest pre-radiation muon

The initial acceptance should fit within LowMC's existing
`PhotonPhaseSpace` representation. Any additional photon-energy windows used
for the soft-limit scan should be imposed as common analysis cuts.

## Recoil Convention

PHOTOS changes daughter momenta to conserve four-momentum, while LowMC
currently preserves the input charged-particle momenta. Primary angular
comparisons must therefore use the shared pre-radiation Born muon directions.

PHOTOS angles relative to its post-radiation muons should be recorded as a
secondary diagnostic. The difference between pre- and post-radiation
definitions should decrease with the maximum accepted photon energy and gives
an explicit measure of recoil effects.

## Normalization

The primary rate observable is the mean number of accepted real photons per
input Born decay. The comparison must distinguish:

- PHOTOS's probability for leaving an event without a resolved photon
- the PHOTOS generation cutoff
- the common analysis acceptance
- LowMC's integrated event yield `lambda_event`
- any LowMC enhancement or oversampling factor

The first physical-rate comparison should not use enhanced event output. If
weighted or oversampled LowMC photons are used for variance reduction, their
weights must be applied explicitly and the unweighted PHOTOS denominator must
remain the number of input Born decays.

## Observables

Compare both absolute rates and normalized shapes:

- accepted photons per Born decay
- `dN / dE_gamma`
- `dN / dlog(E_gamma)`
- angle to the nearest Born muon
- wide-angle photon distribution
- `eta` and `phi` projections in the common acceptance
- LowMC/PHOTOS ratio for every reported distribution

The existing `lowmc_delphi_muon_fsr_diagnostic` should be reused where
practical, but the PHOTOS benchmark must add absolute normalization and a
reproducible external-generator comparison rather than only a shape overlay.

## Soft-Limit Scan

Repeat the analysis for a nested series of decreasing upper photon-energy
cuts while keeping the lower analysis threshold and all other settings fixed.
The precise windows should be chosen after a small numerical pilot confirms
adequate PHOTOS statistics and stable cutoff separation.

The central expected result is

```text
LowMC / PHOTOS -> 1 as E_gamma_max / m_Z -> 0.
```

Rate and shape ratios should be shown as functions of the upper energy cut.
Harder windows are useful controls, but disagreement there is not a LowMC
failure: PHOTOS includes recoil, collinear terms, exact phase-space mappings,
and possibly other corrections outside the leading-soft approximation.

## Validation Stages

1. Verify the fixed Born kinematics and charge-conservation diagnostics.
2. Check the LowMC current against the analytic two-muon soft current at fixed
   photon phase-space points.
3. Run a small PHOTOS sample to confirm configuration, event-record output,
   cutoff behavior, and photon identification.
4. Compare normalized photon shapes in one soft window.
5. Compare the absolute accepted yield per decay.
6. Run the nested soft-limit scan with production statistics.
7. Repeat with independent seeds and record statistical uncertainties.
8. Only after the baseline passes, enable PHOTOS features beyond the shared
   model one at a time to quantify their effects.

## Success Criteria

The benchmark is successful when:

- common inputs and conventions are fully recorded
- LowMC agrees with the analytic two-muon current point by point
- no unexplained normalization factor remains between the programs
- LowMC/PHOTOS rate and shape ratios move toward unity as the upper photon
  energy decreases
- the two softest statistically useful windows are compatible with unity under
  a tolerance fixed after the pilot and before the production run
- recoil-sensitive differences behave consistently with the decreasing soft
  scale
- results reproduce with independent random seeds

Any residual discrepancy must first be classified as physics-model,
acceptance, normalization, recoil, event-record, or sampling related before
changing the LowMC density or defaults.

## Deliverables

- an optional PHOTOS-enabled benchmark executable or helper outside the
  dependency-free core
- machine-readable run configuration and summary output
- per-bin comparison data in CSV or an equivalent simple format
- rate, shape, ratio, and soft-limit-convergence plots
- a short recorded conclusion in `docs/STATUS.md` and the relevant algorithm
  or decision document

The PHOTOS dependency and production-statistics run should remain opt-in and
outside the default CTest path. Lightweight configuration and `--help` checks
may be included in the default test suite.

## Follow-Up Sampling Study

After the physical comparison is established, use the same two-muon benchmark
to test the PHOTOS-inspired analytic leg-centered proposal described in
`docs/ROADMAP.md`. That proposal must remain diagnostic-only until it has also
passed boosted-toy and realistic Pythia-event correctness and performance
checks. A Vincia-inspired charged-pair proposal is a later exploratory step.

## Initial Pilot Status

The first opt-in implementation is now available as
`lowmc_photos_zmumu_benchmark` when `LOWMC_ENABLE_PHOTOS=ON`. It constructs an
exact fixed `Z -> mu+ mu-` HepMC3 decay, runs PHOTOS 3.64, writes every PHOTOS
photon to CSV, and writes a one-row configuration and rate summary. It records
angles relative to both the common Born muons and PHOTOS's recoil-modified
muons.

The shared baseline explicitly uses interference on, exponentiation off,
double and higher bremsstrahlung off, pair emission off, and the exact
`Z` matrix-element correction off. The PHOTOS fixed-order implementation
aborts when its generation cutoff is pushed below roughly `3e-4` for this
decay. The pilot therefore uses `xphcut = 3e-4` (about `27.4 MeV`) and a common
rectangular acceptance `0.05 <= kT < 0.2 GeV`, `|eta| < 2`, and full `phi`.

With `10000` Born decays and seeds `1802, 9373`, PHOTOS generated exactly one
photon in `4219` events and none in `5781`; no multiphoton event occurred. In
the common acceptance it produced `791` photons, or
`0.0791 +/- 0.0028` photons per decay. The L9 factorized LowMC integration is
`0.0801377`, with an L8/L9 relative difference of `7.2e-5`, giving
`PHOTOS / LowMC = 0.987`. This is agreement at the statistical precision of
the pilot, not yet the planned production soft-limit scan.

The pointwise LowMC current agrees with the independent analytic two-muon
expression to a maximum relative difference of `7.5e-10` over generated
photons. PHOTOS conserves four-momentum with a maximum component residual of
`4.5e-7 GeV` in this sample. Validation stages 1 through 3 and the first pilot
absolute-rate check are therefore complete. The production soft-window scan
and its final independent-seed statistics remain open.

The benchmark now also writes an independent weighted LowMC shape sample. The
first `100000`-trial sample has mean weight `0.0801946`, consistent with the
deterministic integral, and effective sample fraction `0.959`. The companion
`scripts/plot_photos_zmumu_comparison.py` produces normalized `kT`, `eta`,
`phi`, and nearest-Born-muon angle distributions, binwise PHOTOS/LowMC ratios,
and machine-readable bin and approximate chi-square summaries.

For the first PHOTOS seed, the approximate chi-square per degree of freedom is
`1.38`, `1.09`, `1.26`, and `1.73` for `kT`, `eta`, `phi`, and nearest-muon
angle respectively. A second independent PHOTOS run with seeds
`24680, 13579` gives `797` accepted photons, rate `0.0797 +/- 0.0028`, and
`PHOTOS / LowMC = 0.995`; its corresponding shape values are `0.90`, `0.56`,
`0.64`, and `0.72`. The plots therefore show no reproducible shape discrepancy
at current PHOTOS statistics. These chi-square values are diagnostic only:
the histograms are normalized, bin correlations are not included, and a
half-count uncertainty floor is used for empty PHOTOS bins.

Validation stage 4 and an initial independent-seed check are now complete.
The next step is a higher-statistics, predeclared energy-window scan rather
than further interpretation of the present statistics-limited ratio bins.

## Initial Soft-Energy Scan

The benchmark now accepts a common `--energy-max` cut in the fixed `Z` rest
frame. The LowMC reference integrates the resulting non-rectangular boundary
exactly in the soft `kT` scaling: at each angular point the allowed logarithmic
`kT` interval ends at `min(kT_max, E_gamma_max / cosh(eta))`. The remaining
angular integral uses the same L9 adaptive map with an L8 convergence check.

The predeclared scan uses
`E_gamma_max = 0.06, 0.075, 0.10, 0.15, 0.25, 0.50 GeV`, three PHOTOS seed
pairs, and `100000` Born decays per seed and energy window. The aggregate
results are:

| `E_gamma_max` [GeV] | PHOTOS rate | LowMC rate | LowMC / PHOTOS |
| ---: | ---: | ---: | ---: |
| 0.060 | 0.008873 | 0.008825 | 0.995 +/- 0.019 |
| 0.075 | 0.020723 | 0.020934 | 1.010 +/- 0.013 |
| 0.100 | 0.036637 | 0.037106 | 1.013 +/- 0.010 |
| 0.150 | 0.059933 | 0.060326 | 1.007 +/- 0.008 |
| 0.250 | 0.078103 | 0.078844 | 1.009 +/- 0.007 |
| 0.500 | 0.079330 | 0.080056 | 1.009 +/- 0.007 |

Every point is statistically compatible with unity, and the softest point is
within `0.5%` of unity. The L8/L9 integration difference stays below `1e-4`
throughout. The scan therefore supports the expected LowMC/PHOTOS convergence
in the shared leading-soft domain and finds no unexplained normalization
factor at present precision.

The nested points are strongly correlated because the same PHOTOS seeds and
successively inclusive energy cuts are used. Consequently, the small
point-to-point changes must not be interpreted as a fitted monotonic trend.
The appropriate conclusion is compatibility with the soft-limit prediction,
not a measured slope. Validation stages 5 through 7 are now complete for this
baseline configuration. A later production refinement could add substantially
more statistics to the softest window or use disjoint energy bins and a full
covariance matrix, but that is not required to establish the present baseline.
