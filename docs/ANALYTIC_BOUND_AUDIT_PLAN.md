# Analytic Low-Bound Audit Plan

## Scope

Test the closed-form global and rectangular cell-local bounds proposed in
`docs/analytic-dipole-envelope-proposal-2026-07-23.md` without using them for
sampling or changing the current shape-unweighting envelope.

The audit covers the complete importance weight

```text
w = alpha/(2 pi)^2 * log(kTmax/kTmin) * g(eta,phi) / q_eta_phi
```

and not only the angular function `g`.

## Fixed held-out Pythia block

After bounded exploratory development runs, use:

- Pythia seeds `46801`, `46802`, and `46803`;
- five selected events per seed;
- both hadronic Z and 13 TeV SoftQCD modes;
- the default `0.002 < kT < 0.2 GeV`, `|eta| < 5`, full-phi phase space;
- the existing tuned grid, full `11 x 11` base probes, boosted-cell policy,
  `1000` calibration trials, and `10000` rejection trials;
- analytic global and cell bounds as read-only diagnostics.

This gives thirty held-out generated events.

## Correctness gate

Every event must have:

- zero invalid analytic legs;
- zero fallback cells;
- zero probe/calibration points above the analytic cell bound;
- zero independent rejection trials above the analytic cell bound;
- maximum observed/bound and trial/bound ratios no greater than one.

Any failure keeps the analytic construction out of further envelope work.

## Practical gates

A standalone analytic replacement for the current empirical envelope is worth
prototyping only if, separately in each process:

- the median analytic/current cell-envelope integral ratio is at most `2`;
- no event has a ratio above `10`; and
- the analytic construction has no material runtime regression relative to
  current base probing.

If the safety gate passes but the tightness gate fails, retain the analytic
bound only as a possible audit, fallback, or branch-and-bound ingredient. Do
not implement the cosh-cos inverse-CDF sampler or replace the current envelope.

The global absolute-value bound is reported for completeness but is not a
standalone candidate unless its scale is unexpectedly competitive.

## Rapidity check

True rapidity is the exact fixed-kT center of the denominator. Its effect on
proposal construction must be tested as a separate A/B study; this audit does
not change the current pseudorapidity-centered refinement policy.

## Held-out result

The fixed held-out block was run on 2026-07-23. All thirty events passed the
correctness gate:

| Process | Events | Charged particles (min/median/max) | Invalid legs | Fallback cells | Probe/calibration violations | Trial violations | Maximum observed/bound | Maximum trial/bound |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| hadronic Z | 15 | 16 / 22 / 36 | 0 | 0 | 0 | 0 | 0.1278 | 0.1212 |
| SoftQCD | 15 | 38 / 112 / 324 | 0 | 0 | 0 | 0 | 0.03021 | 0.03010 |

The practical gates failed:

| Process | Median cell/current | 90th percentile | Maximum | Median current acceptance | Median analytic acceptance | Analytic/base-probe construction time |
|---|---:|---:|---:|---:|---:|---:|
| hadronic Z | 8.34 | 13.98 | 86.96 | 0.102 | 0.0123 | 0.094 |
| SoftQCD | 40.52 | 69.39 | 1339.38 | 0.0672 | 0.00156 | 0.322 |

The median global-bound/current-envelope ratios were approximately
`4.6e9` and `6.2e9`, respectively. The cell construction itself is relatively
cheap, but its acceptance loss makes a standalone analytic rejection sampler
uncompetitive.

The controlled two-pion and text-fixture checks also had zero violations.
Several modest fixtures produced analytic/current ratios below one, while the
boosted and more complex cases became looser. This is consistent with the
held-out conclusion rather than evidence for a broadly useful replacement.

Decision: retain the formulas and diagnostics as an audit or possible future
cell-certification/branch-and-bound ingredient. Do not replace the current
sampler and do not implement the cosh-cos inverse-CDF proposal on this basis.
