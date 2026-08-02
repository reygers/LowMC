# WA102 Fixed-Target Soft-Photon Benchmark

## Purpose

This benchmark compares LowMC's coherent Low current with the inner
bremsstrahlung calculation published by WA102 for minimum-bias fixed-target
`pp` interactions. The direct-integral/LowMC agreement is the implementation
validation target. WA102 provides a complementary external shape comparison,
but its absolute normalization is generator- and reconstruction-dependent.
The substantially larger background-subtracted photon signal is included for
context and is not a LowMC pass criterion.

Primary reference:

- A. Belogianni et al., *Observation of a soft photon signal in excess of QED
  expectations in pp interactions*, Phys. Lett. B 548 (2002) 129-139,
  DOI: 10.1016/S0370-2693(02)02837-X.

The digitized Figure 2d values are stored in
`fixtures/reference/wa102_figure2d_pt_digitized.csv`. The filled
`Data - BG` markers and their statistical errors were digitized independently.
The open `Brems` triangles were digitized from the supplied high-resolution
figure. Both are labelled as digitized whenever plotted.

## Published Target

The WA102 result uses:

- a `450 GeV/c` proton beam on a proton at rest;
- minimum-bias interaction triggers;
- fewer than eight reconstructed charged tracks;
- photons with `0.2 < E_gamma < 1.0 GeV`;
- photon polar angle `theta_gamma < 20 mrad` relative to the beam;
- ten `2 MeV/c`-wide transverse-momentum bins from `0` to `20 MeV/c`;
- normalization per selected interaction.

Table 1 quotes an inner-bremsstrahlung integral of

```text
0.0116 +/- 0.0002 photons / event
```

with a `5%` systematic uncertainty. The background-subtracted signal is

```text
0.0473 +/- 0.0018 +/- 0.0091 photons / event.
```

The digitized data markers are put on the per-event scale by requiring their
ten-bin sum to equal the published `0.0473`. The same single scale factor is
applied to the digitized bremsstrahlung markers. Their displayed ten-bin sum
is about `0.01055 photons/event`; it is not renormalized to Table 1. The
difference is consistent with finite transverse-momentum resolution moving
events across the displayed `20 MeV/c` boundary.

## Generator-Level Approximation

WA102 used a data-tuned FRITIOF sample. The reproducible LowMC benchmark uses
Pythia 8 with `SoftQCD:inelastic` at the same fixed-target beam momentum and
disables generator QED radiation. It applies `N_ch < 8` to all final stable
charged particles as a truth-level proxy for the reconstructed-track cut.

The coherent charged system is:

```text
incoming beam proton + incoming target proton
    + all final stable charged particles
```

This is charge conserving event by event. LowMC evaluates the exact Low
current. WA102 used the Haissinski form for incomplete reconstructed events
but reports agreement with the exact Low formula for complete simulated
events.

The generator replacement and truth-level multiplicity definition are
material approximations. Event-to-event uncertainty and generator dependence
must not be confused with LowMC photon-integration uncertainty.

## Approximate Isolation Model

The published calculation includes the experimental isolation requirement:
the photon conversion point in the lead sheet must be more than `3 mm` from
every charged track. The paper specifies a `50 x 50 cm2` lead sheet located
`66 cm` downstream of the target centre and a `1.2 T` magnetic field
perpendicular to the beam.

The benchmark implements the following explicit approximation:

1. the interaction is at the target centre;
2. the lead sheet is a plane `66 cm` downstream;
3. the magnetic field is uniform at `1.2 T` from the vertex to the sheet;
4. photons propagate in straight lines;
5. final stable charged particles start at the primary vertex and follow
   helical trajectories;
6. only the first charged-track intersection inside the `50 x 50 cm2` sheet
   is considered;
7. the isolation distance is Euclidean distance in the sheet plane;
8. a photon is retained only when its nearest track is farther than `3 mm`.

The lead distance, field, isolation radius, and sheet half-size are command
line parameters. The nominal geometry should be accompanied by `36 cm` and
`96 cm` lead-distance variations, corresponding to the ends of the `60 cm`
target, and by a zero-field diagnostic. These variations probe the missing
vertex distribution and field map; they are not a detector simulation.

## Transverse-Momentum Resolution

WA102 quotes an average `pT` resolution of `2.9 MeV/c`. The benchmark records
truth-level, isolated, and isolated-plus-smeared spectra. Smearing is applied
as an analytic folded-Gaussian bin-response matrix, so reconstructed `pT`
remains non-negative and no additional Monte Carlo randomness is introduced.
Probability migrating above `20 MeV/c` is retained as overflow rather than
renormalized back into the displayed range.

## Validation Contract

For identical selected Pythia events, the executable compares:

- direct numerical integration of the Low density;
- weighted LowMC photon sampling;
- each method before isolation, after isolation, and after isolation plus
  resolution smearing.

The direct and weighted estimates must agree within their combined numerical
uncertainties. The paper-facing external comparison uses the no-isolation
spectrum convolved with the quoted average `2.9 MeV/c` resolution. This is
explicitly not the same selection as the published WA102 markers, which
include isolation. The primary shape window is `2 < pT < 18 MeV/c`: the first
bin is isolation-sensitive, while the final digitized bremsstrahlung marker is
too close to zero for a stable ratio.

The approximate isolation result is retained as a separate sensitivity
diagnostic. It is not applied as a correction and is not shown as the primary
LowMC result.

Agreement with the anomalous `Data - BG` markers is neither expected nor used
as a validation criterion.

## Six-Seed Production Result

The frozen production block uses seeds `24680`, `13579`, `97531`, `86420`,
`11223`, and `44556`, with `300` selected events and `200000` weighted photon
samples per event and seed. The combined `1800`-event result is:

```text
direct truth                         0.00885941 photons / event
direct truth + average smearing      0.00884536 +/- 0.000269 photons / event
weighted LowMC truth + smearing      0.00884820 +/- 0.000269 photons / event
LowMC integration error              0.00000244 photons / event
LowMC / direct                       1.000321
direct coarse/fine shift             0.284%
WA102 Table 1                        0.0116 photons / event
WA102 Figure 2d displayed sum        0.0105461 photons / event (digitized)
```

The paired numerical check succeeds: direct integration and LowMC sampling
agree much more closely than the event-to-event uncertainty. In the paper
window `2 < pT < 18 MeV/c`, LowMC/direct is `1.000268`, the direct coarse/fine
shift is `0.247%`, and the conditional LowMC integration error is
`1.72e-6 photons/event`.

The normalized WA102 shape agrees well in that window: the shape L1 distance
is `0.03096`. The absolute Pythia prediction is lower: the direct integral is
`0.7388` times the digitized WA102 bremsstrahlung sum in the same window. This
normalization difference is not attributed to LowMC. WA102 used data-tuned
FRITIOF events and fewer than eight reconstructed tracks, whereas the
reproduction uses Pythia 8 and fewer than eight truth-level final charged
particles. The first `100`-event pilot happened to lie high; the six-seed
result supersedes it and demonstrates why event-to-event statistics must be
kept separate from photon-integration precision.

## Isolation Sensitivity Diagnostic

The literal all-truth-track isolation proxy is not an adequate reproduction
of the experimental isolation efficiency. In the original `100`-event
geometry scan with seed `24680`, the direct isolated-plus-smeared result was:

```text
lead distance 0.36 m                 0.00307854 photons / event
lead distance 0.66 m                 0.00605932 photons / event
lead distance 0.96 m                 0.00785087 photons / event
zero field at 0.66 m                 0.00600128 photons / event
```

This strong vertex-distance dependence and the remaining suppression relative
to WA102 show that reconstructed-track acceptance and the actual target/field
geometry matter. The isolated curve must therefore remain labelled as an
approximate sensitivity study. It should not replace the unisolated comparison
as evidence about the Low current itself, and the isolation model must not be
tuned to the published rate without independent detector information.
