# Analytic Leg-Centred Proposal

## Purpose

This workstream tests whether a PHOTOS-inspired analytic angular proposal can
reduce LowMC's dependence on event-specific adaptive grids. It is a sampling
study, not an extension of the PHOTOS physics comparison. The exact coherent
Low density remains the target for every weighted sample.

The prototype is diagnostic-only. It does not replace the current grid path,
change the shape-unweighting envelope, or alter reduced-probe defaults while
their independent review is pending.

## Proposal

For charged leg `i`, define `beta_i = |p_i| / E_i` and the normalized
full-sphere component

```text
g_i(Omega) = beta_i / [2 pi log((1 + beta_i)/(1 - beta_i))]
             x 1 / (1 - beta_i cos(theta_i_gamma)).
```

The `beta -> 0` limit is `1 / (4 pi)`. The proposal mixes these components
using absolute particle charge and adds a uniform angular floor:

```text
q_Omega = f_uniform / (4 pi)
          + (1 - f_uniform) sum_i w_i g_i(Omega),

w_i = |Q_i| / sum_j |Q_j|.
```

Since `d Omega = sech^2(eta) d eta d phi`, the proposal density passed to the
existing fixed-`kT` weighted sampler is

```text
q_eta_phi(eta, phi) = q_Omega sech^2(eta).
```

The component polar angle is sampled by analytic inversion and its azimuth
around the charged leg is uniform. The physical weight is always calculated
from the exact coherent current:

```text
w_gamma = [alpha / (2 pi)^2] S_exact(k) / q_eta_phi(eta, phi).
```

The proposal is normalized on the full photon sphere. Samples outside the
configured `eta, phi` acceptance receive zero target weight. This preserves an
analytic normalization and makes acceptance loss explicit. The uniform floor
guarantees support. `minimumOneMinusBeta` is a proposal-only numerical
regulator for inputs that are massless at machine precision; it does not alter
the Low current.

## First Validation

Unit coverage checks the full-sphere normalization numerically, analytic
inverse sampling of the finite-mass denominator, and use of the exact Low
target density in the weighted sample.

The first three-seed comparison uses `100000` samples per seed at fixed
`kT = 0.02 GeV`, a leg-proposal uniform fraction of `0.10`, and the tuned
`32 x 32`, three-refinement-level mass grid as the reference proposal.

| Fixture | Proposal | Mean weight | ESS/N | Maximum weight |
| --- | --- | ---: | ---: | ---: |
| `event4` | leg | `272.094 +/- 0.181` | `0.883` | `520` |
|  | grid | `272.067 +/- 0.405` | `0.600` | `5637` |
| `clustered8` | leg | `436.635 +/- 0.474` | `0.739` | `1203` |
|  | grid | `436.002 +/- 1.146` | `0.325` | `45676` |
| `incoming_outgoing4` | leg | `265.894 +/- 0.227` | `0.821` | `486` |
|  | grid | `265.693 +/- 0.240` | `0.804` | `2744` |
| `hepmc3_validation` | leg | `18.896 +/- 0.036` | `0.484` | `84` |
|  | grid | `18.838 +/- 0.007` | `0.957` | `53` |

The independently sampled means agree within about two combined standard
errors in every fixture. The analytic proposal is substantially better for
the boosted and clustered fixtures, while the grid is clearly better for the
low-boost HepMC3 validation fixture. The result therefore supports a useful
analytic diagnostic path, not a universal replacement for the grid.

## Gate Before Any Envelope Study

Before using `q_eta_phi` in a rejection envelope, the residual
`S_exact / q_eta_phi` must be studied on independent Pythia-produced events.
The required diagnostics are upper-tail quantiles and maxima, charged
multiplicity and boost dependence, acceptance loss, total runtime, and
proposal evaluation cost. No production envelope or hybrid policy should be
implemented until those tails are understood.

## Independent Pythia Tail Gate

The first independent event study is now complete. It uses three Pythia seeds,
twenty valid events per seed, and `20000` weighted samples per proposal and
event at `kT = 0.02 GeV`. Two event classes are included:

- `e+ e- -> gamma*/Z -> u,d,s,c,b` at the Pythia `Z` mass
- non-diffractive `pp` SoftQCD at `sqrt(s) = 13 TeV`

Both use the existing adapter selection of incoming beams plus final stable
charged particles. The comparison records per-event residual quantiles and
maxima, angular-acceptance loss, selected charged multiplicity, outgoing boost
and rapidity markers, proposal construction time, and sampling time.

| Mode | Proposal | Median ESS/N | Minimum ESS/N | Largest max/mean | Effective samples/s |
| --- | --- | ---: | ---: | ---: | ---: |
| `ee-hadrons` | leg | `0.334` | `0.188` | `17.4` | `2.41e5` |
|  | grid | `0.868` | `0.00040` | `7090` | `5.89e4` |
| `softqcd` | leg | `0.323` | `0.174` | `35.0` | `8.83e4` |
|  | grid | `0.932` | `0.253` | `150` | `5.30e4` |

The grid usually gives a higher per-sample effective fraction, especially in
SoftQCD. The analytic proposal is nevertheless fast enough to yield a higher
median effective-sample throughput in both modes. More importantly for an
envelope study, its observed residual tail is much more stable: no analytic
event exceeds `max(weight) / mean(weight) = 35.1`, whereas rare grid events
reach several thousand.

The median paired yield difference is `-0.30%` in `ee-hadrons` and `+0.18%` in
SoftQCD. One nominal greater-than-three-standard-error pair occurs in each
mode out of sixty. Focused `200000`-sample reruns show why ordinary pull
counting is insufficient for the outliers: the problematic grid samples have
extreme, still-growing maxima and underestimated ordinary standard errors.
For the worst SoftQCD pair, three independent analytic reruns give means
`364.90`, `362.90`, and `365.21`, while their maximum-to-mean ratios remain
near `22`. The corresponding grid estimates range from `352.28` to `360.99`
and develop maximum-to-mean ratios from `1805` to `3665`.

This passes the intended gate for implementing a separate diagnostic analytic
envelope prototype. It does not establish a guaranteed majorant: the observed
maximum is empirical, and any next prototype must retain strict failure
detection and independent validation. The current production grid and
reduced-probe defaults remain unchanged.

## Empirical Analytic-Envelope Prototype

The first strict diagnostic envelope is now implemented. For each event it:

1. draws `20000` calibration points from the analytic leg-centred proposal;
2. sets a global residual bound to a predeclared safety factor times the
   largest calibration weight;
3. draws `100000` independent trial points; and
4. declares the envelope failed if even one trial weight exceeds the bound.

No photon is emitted by this diagnostic and an overshoot is never silently
clipped into a valid result. The predeclared safety-factor scan uses `1.2`,
`1.5`, and `2.0`, three independent Pythia seeds, and twenty events per seed
in each of the two event classes used for the tail gate.

| Mode | Safety factor | Events | Trial points | Failures | Largest trial/bound | Mean expected acceptance |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `ee-hadrons` | `1.2` | `60` | `6.0e6` | `0` | `0.844` | `7.97%` |
|  | `1.5` | `60` | `6.0e6` | `0` | `0.675` | `6.38%` |
|  | `2.0` | `60` | `6.0e6` | `0` | `0.507` | `4.78%` |
| `softqcd` | `1.2` | `60` | `6.0e6` | `0` | `0.915` | `6.58%` |
|  | `1.5` | `60` | `6.0e6` | `0` | `0.732` | `5.26%` |
|  | `2.0` | `60` | `6.0e6` | `0` | `0.549` | `3.95%` |

The observed acceptance agrees with the expected acceptance to within the
trial statistics. The proposal itself is essentially free to construct on
this scale; the mean `20000`-point calibration time is about `51 ms` for
hadronic Z events and `130 ms` for SoftQCD on the machine used for this scan.

The `1.2` candidate therefore passes this finite diagnostic sample and is the
best of the predeclared factors on acceptance. This remains an empirical
majorant, not a mathematical upper bound. It is suitable as an opt-in
reference for further stress tests, with mandatory overshoot failure, but is
not ready to replace the production cell-local grid envelope. In particular,
the present evidence does not cover arbitrary multiplicity, angular windows,
`kT`, or event generators.
