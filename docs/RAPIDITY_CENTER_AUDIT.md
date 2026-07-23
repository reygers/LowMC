# Charged-Particle Eta-Center Audit

## Question

The exact fixed-`kT` denominator

```text
D_i = mT_i cosh(eta - y_i) - pT_i cos(phi - phi_i)
```

is centered at true rapidity `y_i`, whereas LowMC's charged-direction targeting
historically used pseudorapidity. This audit asks whether moving only those
targeting centers to true rapidity improves the current shape-unweighting
envelope.

The adaptive grid itself is built from exact Low intensities and is not
particle-centered. The center convention enters only optional targeted
refinement, targeted base probes, boosted-cell identification, and
charged-near audit labels. It does not enter `LowCurrent`, the physical photon
density, event yield, multiplicity, or the acceptance formula.

## Diagnostic contract

`ShapeUnweightingConfig::chargedParticleEtaCenter` selects
`Pseudorapidity` or `Rapidity`. Pseudorapidity remains the default. Paired runs
use identical particles, phase space, seeds, calibration counts, and trial
counts. Unit tests evaluate the exact Low density at a fixed point set and
confirm identical values.

The generated-event block uses Pythia seeds `57911`, `57912`, and `57913`,
five events per seed, for both hadronic Z and SoftQCD production. Each event is
run once with each center convention, using `1000` calibration trials and
`10000` rejection trials. This gives thirty paired events.

## Controlled result

For slow charged pions with `p = 0.1 GeV` and pseudorapidity `+/-2`, the
particle pseudorapidity and true rapidity differ strongly. With the current
boosted-cell policy:

| Center | Envelope integral | Expected acceptance | Promoted cells | Probe evaluations | Failures |
|---|---:|---:|---:|---:|---:|
| pseudorapidity | 0.0423414 | 0.2695 | 46 | 293383 | 0 |
| rapidity | 0.0465221 | 0.2456 | 48 | 294265 | 0 |

True rapidity moves the conservative factor-20 boosted-cell treatment onto the
exact peak and makes the envelope about `9.9%` larger. This is safer targeting,
not more efficient targeting.

With two levels of explicit local refinement and the fixed probe policy, the
same slow fixture gives envelope integrals `0.0381817` and `0.0381731`: a
`0.02%` difference. For ordinary `p = 10 GeV` pions at pseudorapidity `+/-2`,
both modes are numerically identical in the tested configurations.

## Paired Pythia result

| Process | Pairs | Median particle `|eta-y|` | Maximum event `|eta-y|` | Median rapidity/pseudorapidity envelope | Maximum ratio | Median probe difference | Median promoted-cell difference | Failures pseudo/rapid | Median runtime ratio |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| hadronic Z | 15 | 0.00441 | 1.516 | 1.00031 | 1.01215 | 0 | 0 | 0 / 0 | 0.998 |
| SoftQCD | 15 | 0.0973 | 3.815 | 1.01765 | 1.06795 | 6174 | 14 | 0 / 0 | 1.013 |

The envelope ratio implies essentially unchanged acceptance for hadronic Z and
a median `1.7%` acceptance loss for SoftQCD. Direct expected- and
observed-acceptance estimates fluctuate at the sub-percent level because the
two envelopes sample different trial points; the envelope integral is the
more stable paired comparison.

## Decision

True rapidity is the correct mathematical center of an individual massive-leg
denominator, but changing the center alone does not improve LowMC's current
grid/envelope algorithm. It is neutral for relativistic hadronic-Z particles
and modestly worsens the conservative SoftQCD envelope.

Keep pseudorapidity as the production default. Retain true rapidity as an
explicit diagnostic option and use it in future analytic constructions where
the exact denominator center matters. A future redesign of the boosted-cell
safety policy may revisit this result; the present audit does not justify a
default change.
