# Scripts

Optional helper scripts for diagnostics live here. They are not part of the
LowMC core and may depend on Python packages that are not required for building
or testing the C++ library.

## Plot the two-particle diagnostic

The diagnostic uses a transverse back-to-back charged-pion toy event by default,
with particles pointing near `(eta, phi) = (0, pi/2)` and `(0, 3pi/2)`.
The default momentum is `1 GeV`, where the finite-mass dead-cone scale is
visible in the heatmap.

Generate the CSV with CTest or by running the diagnostic directly:

```sh
cmake --build build
./build/lowmc_two_particle_diagnostic build/two_particle_eta_phi.csv
```

Create a PNG heatmap:

```sh
mkdir -p build/matplotlib
MPLCONFIGDIR=build/matplotlib \
python3 scripts/plot_eta_phi.py build/two_particle_eta_phi.csv build/two_particle_eta_phi.png --log
```

## Momentum comparison

Generate comparison CSV files for several charged-pion momenta:

```sh
for p in 1 10 100 1000; do
  ./build/lowmc_two_particle_diagnostic build/two_particle_eta_phi_p${p}.csv --momentum ${p}
done
```

Create matching PNG heatmaps:

```sh
mkdir -p build/matplotlib
for p in 1 10 100 1000; do
  MPLCONFIGDIR=build/matplotlib \
  python3 scripts/plot_eta_phi.py build/two_particle_eta_phi_p${p}.csv build/two_particle_eta_phi_p${p}.png --log
done
```

Suggested interpretation:

- `p = 1 GeV`: visual sanity check with visible finite-mass scale
- `p = 10 GeV`: moderately boosted case
- `p = 100 GeV`: high-energy case
- `p = 1000 GeV`: stress test for angular resolution

The plotting script requires matplotlib. The C++ core does not.

The fixed grid has configurable point counts. The defaults are unchanged at
`101 x 121`; use finer grids for convergence checks:

```sh
./build/lowmc_two_particle_diagnostic build/two_particle_eta_phi_p10_fine.csv \
  --momentum 10 \
  --eta-bins 401 \
  --phi-bins 481
```

## DELPHI hadronic bremsstrahlung summary

Aggregate one or more event-level CSV files from
`lowmc_pythia8_delphi_hadronic_benchmark` and compare the unscaled result with
the published DELPHI Table 4 curve:

```sh
MPLCONFIGDIR=build-pythia8/matplotlib \
python3 scripts/summarize_delphi_hadronic_benchmark.py \
  build-pythia8/delphi_hadronic_seed*.csv \
  --summary build-pythia8/delphi_hadronic_summary.csv \
  --plot build-pythia8/delphi_hadronic.png \
  --joint-inputs build-pythia8/delphi_hadronic_seed*_joint.csv \
  --theta-summary build-pythia8/delphi_hadronic_theta_summary.csv \
  --theta-plot build-pythia8/delphi_hadronic_theta.png
```

The summary keeps event-to-event jackknife uncertainty separate from the
conditional LowMC photon-integration uncertainty. Its diagonal chi-square is
a compact diagnostic only because DELPHI and generator-event correlations are
not available. The angular reference consists of vector-digitized Figure 4b
markers; it is used for an unscaled shape and normalization cross-check, not a
formal chi-square.

## Direction and acceptance examples

The first toy particle direction can be changed with `--particle-eta` and
`--particle-phi`. The second particle is constructed back-to-back in 3D.

Default midrapidity transverse pair:

```sh
./build/lowmc_two_particle_diagnostic build/two_particle_midrapidity.csv \
  --momentum 1 \
  --particle-eta 0 \
  --particle-phi 1.5707963267948966
```

Forward particle at `eta = 4`:

```sh
./build/lowmc_two_particle_diagnostic build/two_particle_forward_eta4.csv \
  --momentum 100 \
  --particle-eta 4 \
  --particle-phi 1.5707963267948966 \
  --eta-min 2 \
  --eta-max 5 \
  --phi-min 0 \
  --phi-max 6.283185307179586
```

Plot either CSV with the same `plot_eta_phi.py` command shown above.

## Weighted sampler study

The fixed-`k_T` weighted-sample diagnostic writes individual sampled photons:

```sh
./build/lowmc_fixed_kt_weighted_samples build/fixed_kt_weighted_samples.csv \
  --samples 10000 \
  --momentum 10
```

Pass `--proposal adaptive` to use the first adaptive fixed-`k_T` `eta-phi`
proposal instead of the uniform baseline. The adaptive proposal is a mixture of
a uniform floor over the full acceptance and local boxes around the selected
charged-particle directions. Its diagnostic parameters can be changed with
`--adaptive-uniform-fraction`, `--adaptive-width-scale`,
`--adaptive-min-half-width`, and `--adaptive-max-half-width`:

```sh
./build/lowmc_fixed_kt_weighted_samples build/fixed_kt_adaptive_samples.csv \
  --samples 10000 \
  --momentum 10 \
  --proposal adaptive \
  --adaptive-width-scale 10 \
  --adaptive-min-half-width 0.02
```

Pass `--proposal smooth` to use the first smooth local-kernel diagnostic
proposal. It uses compact triangular kernels around charged-particle
directions plus a uniform floor:

```sh
./build/lowmc_fixed_kt_weighted_samples build/fixed_kt_smooth_samples.csv \
  --samples 10000 \
  --momentum 10 \
  --proposal smooth \
  --smooth-width-scale 10 \
  --smooth-min-half-width 0.02
```

Pass `--proposal leg` to use the analytic finite-mass leg-centred diagnostic.
It is normalized on the full photon sphere and keeps the exact coherent Low
density in every nonzero sample weight:

```sh
./build/lowmc_fixed_kt_weighted_samples build/fixed_kt_leg_samples.csv \
  --samples 10000 \
  --momentum 100 \
  --proposal leg \
  --leg-uniform-fraction 0.1
```

The reproducible multi-fixture comparison with the tuned grid proposal is:

```sh
python3 scripts/compare_leg_centered_proposal.py \
  build/leg_centered_proposal_comparison.csv \
  --benchmark build/lowmc_fixed_kt_weighted_samples
```

Summarize independent calibration/trial runs of the Pythia analytic-envelope
diagnostic with:

```sh
python3 scripts/summarize_pythia_leg_envelope.py \
  build-pythia8/leg_envelope_*.csv \
  --summary build-pythia8/leg_envelope_summary.csv \
  --plot build-pythia8/leg_envelope_study.png
```

The sampler-study diagnostic writes one summary row per toy charged-pion
momentum, including effective sample size:

```sh
./build/lowmc_weighted_sampler_study build/weighted_sampler_study.csv \
  --samples 10000
```

Compare against the adaptive proposal:

```sh
./build/lowmc_weighted_sampler_study build/weighted_sampler_study_adaptive.csv \
  --samples 10000 \
  --proposal adaptive
```

The default momenta are `1`, `10`, `100`, and `1000 GeV`. Repeating
`--momentum` overrides that list:

```sh
./build/lowmc_weighted_sampler_study build/weighted_sampler_forward.csv \
  --samples 10000 \
  --particle-eta 4 \
  --eta-min 2 \
  --eta-max 5 \
  --momentum 10 \
  --momentum 100
```

Plot the sampler-study summary:

```sh
mkdir -p build/matplotlib
MPLCONFIGDIR=build/matplotlib \
python3 scripts/plot_sampler_study.py \
  build/weighted_sampler_study.csv \
  build/weighted_sampler_study.png
```

The plot shows mean weight, weight RMS, and effective sample fraction versus
charged-pion momentum. The momentum axis is logarithmic by default; pass
`--linear-x` for a linear axis.

Run the first adaptive tuning scan:

```sh
python3 scripts/scan_adaptive_tuning.py \
  build/adaptive_tuning_summary.csv \
  --samples 50000
```

By default this compares the uniform baseline with adaptive variants using
`widthScale = 2, 5, 10, 20`, `minHalfWidth = 0.01, 0.02, 0.05`, and
momenta `10, 100, 1000 GeV`. The output CSV includes mean weight, weight RMS,
effective sample fraction, and pull against fine-grid references. Use
`--skip-reference` for a faster ESS-only scan.

Plot and summarize the scan:

```sh
mkdir -p build/matplotlib
MPLCONFIGDIR=build/matplotlib \
python3 scripts/plot_adaptive_tuning.py \
  build/adaptive_tuning_summary.csv \
  build/adaptive_tuning_summary.png \
  --summary-md build/adaptive_tuning_summary.md
```

The PNG shows weight RMS, effective sample fraction, and absolute pull versus
momentum for all scanned variants. The Markdown summary lists the preferred
tested variant per momentum, using `|pull| <= 3` as the default sanity cut.

Run the focused `p = 100 GeV` benchmark check:

```sh
python3 scripts/run_p100_focus.py \
  build/p100_focus_summary.csv \
  --summary-md build/p100_focus_summary.md \
  --samples 100000
```

This performs grid-reference checks at increasingly fine resolutions and then
compares a targeted set of wider local-box adaptive variants against the
finest grid. It also includes first grid-adaptive proposal variants built from
cell-center `S(k)` estimates with high-mass-cell refinement. The default run
treats `p = 100 GeV` as the main high-energy tuning benchmark and keeps
`p = 1000 GeV` as an extreme stress case.

Scan grid-adaptive map parameters:

```sh
python3 scripts/scan_grid_tuning.py \
  build/grid_tuning_summary.csv \
  --summary-md build/grid_tuning_summary.md \
  --samples 100000
```

The scan varies the coarse grid size, refinement depth, refined-cell fraction,
and uniform floor. The output ranks variants by weight RMS, effective sample
fraction, and pull against the current best `p = 100 GeV` grid reference.

Check the `p = 100 GeV` reference quality:

```sh
python3 scripts/run_p100_reference.py \
  build/p100_reference_summary.csv \
  --summary-md build/p100_reference_summary.md \
  --samples 1000000
```

This compares fixed-grid integrals with high-stat weighted grid-adaptive
estimates. The weighted estimates are stochastic, but they use exact Low
weights and provide an independent check when fixed global grids are visibly
under-resolved.

Retabulate an existing grid-adaptive scan against that weighted reference:

```sh
python3 scripts/retabulate_grid_tuning.py \
  build/grid_tuning_summary.csv \
  build/p100_reference_summary.csv \
  build/grid_tuning_retabulated.csv \
  --summary-md build/grid_tuning_retabulated.md
```

This recomputes differences and pulls using an inverse-variance weighted
reference from the `weighted_grid_adaptive` rows. For the current `p = 100 GeV`
diagnostics, the retabulated reference is about `157.54 +/- 0.32`, while the
`1601 x 1921` fixed-grid integral remains a lower cross-check at about
`148.78`.

Compare the tuned `p = 100 GeV` grid map with the broad local-box baseline and
the first smooth triangular-kernel proposal:

```sh
python3 scripts/compare_p100_proposals.py \
  build/p100_proposal_comparison.csv \
  --summary-md build/p100_proposal_comparison.md \
  --samples 100000
```

The summary reports mean weight, weight RMS, effective sample fraction, and
pull against the high-stat weighted reference.

Scan smooth triangular-kernel parameters for the same benchmark:

```sh
python3 scripts/scan_smooth_tuning.py \
  build/smooth_tuning_summary.csv \
  --summary-md build/smooth_tuning_summary.md \
  --samples 100000
```

The scan varies `smoothWidthScale`, `smoothMinHalfWidth`, and
`smoothUniformFraction`, and includes the tuned grid map and broad local-box
baseline as fixed comparison rows.

Compare the current mass-refined grid map with the first variance-aware
grid-refinement diagnostic:

```sh
python3 scripts/compare_variance_grid.py \
  build/variance_grid_comparison.csv \
  --summary-md build/variance_grid_comparison.md \
  --samples 100000
```

The variance-aware grid mode can also be used directly through
`lowmc_weighted_sampler_study` with `--grid-refinement-score variance` and
`--grid-variance-probe-bins n`. The first score implemented is
`area^2 Var[S(k)]` estimated from a small lattice of probe points inside each
candidate cell.

Run the focused multi-seed mass-versus-variance grid comparison:

```sh
python3 scripts/multiseed_variance_grid.py \
  build/variance_grid_multiseed.csv \
  --aggregate-csv build/variance_grid_multiseed_aggregate.csv \
  --summary-md build/variance_grid_multiseed.md \
  --samples 100000 \
  --seeds 5
```

This writes both per-seed rows and aggregate seed-to-seed means and standard
deviations. The first five-seed result indicates that mass-refined and
variance-refined grids are statistically tied for the current `p = 100 GeV`
benchmark.

Compare the full log-uniform-`k_T` weighted path with uniform and tuned
grid-adaptive `eta-phi` proposals:

```sh
python3 scripts/compare_full_weighted.py \
  build/full_weighted_comparison.csv \
  --summary-md build/full_weighted_comparison.md \
  --samples 100000
```

The underlying diagnostic executable is `lowmc_full_weighted_sampler_study`.
By default it uses `k_T = 0.002..0.2 GeV`, samples `k_T` log-uniformly, and
builds the grid proposal once at the geometric mean of the `k_T` range.

Check the full grid-adaptive path over multiple seeds and `k_T` ranges:

```sh
python3 scripts/multirange_full_weighted.py \
  build/full_weighted_multirange.csv \
  --aggregate-csv build/full_weighted_multirange_aggregate.csv \
  --summary-md build/full_weighted_multirange.md \
  --samples 100000 \
  --seeds 5
```

This verifies the expected logarithmic scaling with
`ln(kTMax / kTMin)` and checks that building the grid at the geometric-mean
`k_T` is stable for the current toy benchmark.

## Visual sampler validation

Create a fixed-`k_T` visual comparison between the deterministic
two-particle grid, the tuned grid-adaptive proposal occupancy, and the weighted
sample estimate:

```sh
mkdir -p build/sampler_validation build/matplotlib
for p in 1 10 100; do
  ./build/lowmc_two_particle_diagnostic \
    build/sampler_validation/truth_p${p}.csv \
    --momentum ${p} \
    --eta-bins 401 \
    --phi-bins 481
  ./build/lowmc_fixed_kt_weighted_samples \
    build/sampler_validation/samples_p${p}.csv \
    --momentum ${p} \
    --samples 100000 \
    --proposal grid \
    --grid-eta-bins 32 \
    --grid-phi-bins 32 \
    --grid-refinement-levels 3 \
    --grid-refine-fraction 0.1 \
    --grid-uniform-fraction 0.1
  MPLCONFIGDIR=build/matplotlib \
  python3 scripts/plot_sampler_validation.py \
    build/sampler_validation/truth_p${p}.csv \
    build/sampler_validation/samples_p${p}.csv \
    build/sampler_validation/validation_p${p}.png \
    --title "p = ${p} GeV, fixed kT = 0.02 GeV" \
    --truth-label "401x481 grid"
done
```

The plot shows the truth grid shape, unweighted proposal occupancy, weighted
sample shape, and normalized one-dimensional projections in `eta` and `phi`.
It also adds zoomed projection panels around the strongest grid peak, which
makes the sharp boosted structures easier to compare by eye. The weighted
integral is estimated from the sample mean. The grid integral is only a
reliable reference when the grid resolves the narrow charged-particle
structures. In the current `401 x 481` visual check, `p = 1 GeV` and
`p = 10 GeV` are useful direct comparisons.

For the current `p = 100 GeV` stress plot, label the known high-stat weighted
reference explicitly and warn that the fixed grid is under-resolved:

```sh
MPLCONFIGDIR=build/matplotlib \
python3 scripts/plot_sampler_validation.py \
  build/sampler_validation/truth_p100.csv \
  build/sampler_validation/samples_p100.csv \
  build/sampler_validation/validation_p100.png \
  --title "p = 100 GeV, fixed kT = 0.02 GeV" \
  --truth-label "401x481 grid" \
  --reference-label "weighted reference" \
  --reference-integral 157.54 \
  --reference-uncertainty 0.32 \
  --warning "warning: fixed grid is under-resolved here; use as a stress visualization"
```

Run the minimal event-level weighted API on the two-pion toy event:

```sh
./build/lowmc_weighted_event_example \
  build/weighted_event_example_p100.csv \
  --samples 1000 \
  --seed 12345 \
  --momentum 100 \
  --kt-min 0.002 \
  --kt-max 0.2
```

This writes weighted photon four-vectors and weights. The executable uses the
generator-independent `generateWeightedPhotons` API; it does not read HepMC or
Pythia event records.
The CSV also includes sampled `k_T`, `eta`, `phi`, and proposal density for
validation. A compact C++ API snippet is available in `docs/API_EXAMPLE.md`.

## Compare grid and weighted diagnostics

At fixed `k_T`, the grid diagnostic can be integrated over `eta-phi` and
compared with the mean weight from the weighted-sample diagnostic:

```sh
./build/lowmc_two_particle_diagnostic build/two_particle_eta_phi.csv
./build/lowmc_fixed_kt_weighted_samples build/fixed_kt_weighted_samples.csv \
  --samples 100000
python3 scripts/compare_fixed_grid_weighted.py \
  build/two_particle_eta_phi.csv \
  build/fixed_kt_weighted_samples.csv
```

The script reports the grid integral, weighted-sample mean, standard error,
relative difference, and pull.

## Pythia-produced HepMC3 afterburner validation

For an opt-in manual validation of the Pythia8-to-HepMC3 afterburner path, use:

```sh
python3 scripts/validate_pythia_afterburner.py \
  --tool-dir build-pythia8-hepmc3/adapters/hepmc3 \
  --output-dir build-pythia8-hepmc3/pythia_afterburner_validation
```

The helper wraps `generate_pythia_hepmc3_sample.py`,
`lowmc_hepmc3_afterburner_diagnostic`, and
`lowmc_hepmc3_output_summary_diagnostic`. It writes generated input,
afterburner output, logs, `validation_summary.csv`, and
`validation_events.csv` under the selected output directory. The summary CSV
contains aggregate runtime and maximum-complexity fields. The per-event CSV
records selected charged multiplicity, requested and retained photon counts,
shape-unweighting diagnostics, output-summary comparison fields, true
`event_seconds` from the afterburner diagnostic when available, and the older
estimated afterburner time per event for compatibility. This can take minutes
for multi-event samples and is not part of the default CTest path beyond a
cheap `--help` test.

New afterburner diagnostics also write per-stage timing keys:
`select_seconds`, `yield_seconds`, `poisson_seconds`, `shape_seconds`, and
`write_seconds`. The validation helper preserves these fields in both
`validation_events.csv` and aggregate sums/maxima in `validation_summary.csv`.

To aggregate one or more helper runs by true per-event timing, use:

```sh
python3 scripts/summarize_pythia_afterburner_timing.py \
  build-pythia8-hepmc3/pythia_afterburner_timing_scan_seed24681 \
  build-pythia8-hepmc3/pythia_afterburner_timing_scan_seed24682 \
  --output-csv build-pythia8-hepmc3/pythia_afterburner_timing_scan_events.csv \
  --summary-md build-pythia8-hepmc3/pythia_afterburner_timing_scan_summary.md
```

The aggregator is read-only. It accepts either `validation_events.csv` paths or
directories containing that file, then reports total and maximum event time,
charged-multiplicity and photon-count ranges, stage-time totals when present,
simple correlations with event time, and the slowest event rows. These
correlations are intended as compact manual diagnostics, not as a production
performance model.

## Grid-vs-weighted stress study

The stress-study helper runs matching fixed-grid and weighted-sample
diagnostics for default midrapidity momenta `1`, `10`, `100`, and `1000 GeV`,
plus one forward case at `eta = 4`:

```sh
python3 scripts/run_grid_weighted_stress.py \
  build/grid_weighted_stress_summary.csv \
  --samples 100000
```

It writes per-case input CSV files under `build/grid_weighted_stress/` and
a summary CSV with the grid integral, weighted mean, standard error, pull, and
effective sample fraction.

This is a diagnostic stress test for the current uniform fixed-`k_T`
`eta-phi` proposal and the coarse fixed grid. Large pulls or small effective
sample fractions in boosted cases should be treated as evidence that the
diagnostic resolution or proposal is inadequate, not as a change to the
underlying Low-current formula.

For grid-convergence checks, keep the weighted sample fixed and increase the
fixed-grid point counts:

```sh
python3 scripts/run_grid_weighted_stress.py \
  build/grid_weighted_stress_p10_fine.csv \
  --case mid_p10 \
  --samples 100000 \
  --eta-bins 401 \
  --phi-bins 481
```
