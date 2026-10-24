# API Example

This example uses the generator-independent core API directly. This is the
simplest path for hand-made toy events: provide the charged particles, configure
the photon phase space, and call the weighted generator.

The optional HepMC3 and Pythia8 adapters select charged particles before
calling the same core. Adapter-support selection types are not needed for
manual toy events.

## Two Charged Pions

```cpp
#include "lowmc/LowMC.h"

#include <cmath>
#include <iostream>
#include <vector>

const double pionMass = 0.13957039;
const double momentum = 10.0;
const double energy = std::sqrt(momentum * momentum + pionMass * pionMass);

std::vector<lowmc::ChargedParticle> selectedParticles = {
    {lowmc::FourVector{energy, 0.0, momentum, 0.0}, +1.0, lowmc::InOutSign::Outgoing},
    {lowmc::FourVector{energy, 0.0, -momentum, 0.0}, -1.0, lowmc::InOutSign::Outgoing},
};

lowmc::WeightedEventGeneratorConfig config{
    lowmc::PhotonPhaseSpace{
        0.002,
        0.2,
        -5.0,
        5.0,
        0.0,
        2.0 * std::acos(-1.0),
    },
    1000,
    12345,
    lowmc::tunedMassGridAdaptiveEtaPhiConfig(),
    1, // oversamplingFactor
};

lowmc::WeightedEventGenerationResult result =
    lowmc::generateWeightedPhotons(selectedParticles, config);

std::cout << "weighted photons: " << result.samples.size() << "\n";
std::cout << "mean event weight: " << result.weightSummary.meanWeight << "\n";
std::cout << "effective sample size: "
          << result.weightSummary.effectiveSampleSize << "\n";
std::cout << "signed charge sum: "
          << result.chargeConservation.netCharge << "\n";

for (const lowmc::WeightedEventPhotonSample& sample : result.samples) {
    const lowmc::Photon& photon = sample.photon;

    // Fill histograms with photon.weight.
    // sample.kT, sample.eta, sample.phi, and sample.proposalDensity are
    // retained for diagnostics and validation.
}
```

## Input Conventions

Each `ChargedParticle` uses:

- four-vector order `(E, px, py, pz)`, in GeV
- charge in units of the positron charge, so `pi+` has `+1`
- `inOutSign = Incoming` for incoming external particles
- `inOutSign = Outgoing` for outgoing external particles

The manual example above uses two outgoing charged pions pointing along
opposite transverse directions. Incoming beams can be included by adding
particles with `lowmc::InOutSign::Incoming`.

## Output Conventions

The returned `photons` vector is a convenience view of `result.samples[i].photon`.
Use `result.samples` when proposal density, sampled `k_T`, `eta`, or `phi` are
needed for validation.

`Photon::weight` is the weighted event contribution in the configured photon
phase-space region. Histograms filled with this weight estimate the Low
photon spectrum for the selected event. Averaging such histograms over many
input generator events gives the inclusive prediction.

`WeightedEventGeneratorConfig::oversamplingFactor` can be used to reduce
histogram noise for weighted analyses. If it is larger than one, LowMC returns
`samples * oversamplingFactor` photon samples and divides each photon weight by
`oversamplingFactor`. This preserves the expected event contribution while
providing more histogram entries. It should not be interpreted as a physical
photon-yield enhancement.

`result.chargeConservation.netCharge` reports
`sum_i inOutSign_i Q_i`. A nonzero value is not silently corrected; it should be
understood as either a deliberate selection choice or a sign/selection mistake.

## Existing Executable

The same direct-input pattern is available as a tiny compiled example:

```sh
./build/lowmc_example_weighted_photons
```

Its source is
`examples/basic/generate_weighted_photons.cpp`. It prints a compact summary
with the sample count, mean weight, effective sample size, signed charge sum,
and the first weighted photon.

A CSV-writing diagnostic executable is also available:

```sh
./build/lowmc_weighted_event_example \
  build/weighted_event_example_p10.csv \
  --samples 1000 \
  --oversampling-factor 4 \
  --seed 12345 \
  --momentum 10 \
  --kt-min 0.002 \
  --kt-max 0.2
```

It writes a CSV with photon four-vectors, sampled `k_T`, `eta`, `phi`, proposal
density, and weight.

## Weighted Soft Dielectrons

The leading virtual-photon extension uses the same selected charged particles.
It samples pair mass, pair transverse momentum, pair rapidity, azimuth, and the
two decay angles:

```cpp
lowmc::WeightedDielectronGeneratorConfig pairConfig{
    lowmc::DielectronPhaseSpace{
        0.002,                  // M_min [GeV]
        0.3,                    // M_max [GeV]
        0.002,                  // qT_min [GeV]
        0.2,                    // qT_max [GeV]
        -1.0,                   // y_min
        1.0,                    // y_max
        0.0,
        2.0 * std::acos(-1.0),
    },
    1000,
    12345,
    1,                         // oversamplingFactor
};

lowmc::WeightedDielectronGenerationResult pairs =
    lowmc::generateWeightedDielectrons(selectedParticles, pairConfig);

for (const lowmc::WeightedDielectronSample& sample : pairs.samples) {
    // Fill pair or single-lepton histograms with sample.dielectron.weight.
    const lowmc::FourVector& electron = sample.dielectron.electron;
    const lowmc::FourVector& positron = sample.dielectron.positron;
}
```

The weight is defined for the physical measure
`dM^2 qT dqT dy dphi` and includes the leading virtual Low density, the
angle-integrated Kroll--Wada density, the exact proposal density, and a
normalized polarization-dependent decay-angle correction. The electron and
positron sum exactly to `sample.virtualPhoton` up to floating-point precision.

This mode is a weighted leading-soft baseline. It does not apply recoil to the
hadronic event and does not implement Lichard finite-$q$ terms, form factors,
or subleading process-dependent radiation. A charge-conserving selected
current is required for its Ward identity and polarization interpretation.

The compiled example is:

```sh
./build/lowmc_example_weighted_dielectrons
```

## Shape-Unweighted Candidate API

The shape-unweighting API uses the same direct charged-particle input, but
returns accepted candidates without per-photon analysis weights. These
candidates have the Low-photon shape only; multiplicity and event-record output
remain separate layers.

```cpp
lowmc::ShapeUnweightingConfig unweighting{
    lowmc::PhotonPhaseSpace{
        0.002,
        0.2,
        -5.0,
        5.0,
        0.0,
        2.0 * std::acos(-1.0),
    },
    1000,   // trial photons
    12345,  // seed
    lowmc::tunedMassGridAdaptiveEtaPhiConfig(),
    1000,   // calibrationTrials
    11,     // probeBins
    lowmc::ShapeUnweightingProbeRule::Interior,
    3.0,    // safetyFactor
    lowmc::EnvelopeFailurePolicy::Strict,
};

lowmc::ShapeUnweightingResult unweighted =
    lowmc::generateShapeUnweightedPhotonCandidates(selectedParticles, unweighting);

if (unweighted.status == lowmc::ShapeUnweightingStatus::Ok) {
    for (const lowmc::ShapeUnweightedPhotonCandidate& candidate : unweighted.candidates) {
        // candidate.photon has an ordinary four-momentum and no analysis weight.
        // candidate.kT, candidate.eta, and candidate.phi are retained for validation.
    }
}

std::cout << "accepted candidates: "
          << unweighted.envelope.acceptedCount << "\n";
std::cout << "acceptance efficiency: "
          << unweighted.envelope.observedAcceptanceEfficiency << "\n";
std::cout << "envelope failures: "
          << unweighted.envelope.failedEnvelopeCount << "\n";
```

These candidates only have the correct sampled Low-photon shape. They are not,
by themselves, a complete physical photon multiplicity model. The provisional
HepMC3 afterburner composes the separately implemented event-yield, Poisson
multiplicity, shape-unweighting, and event-record layers.

With `EnvelopeFailurePolicy::Strict`, any observed envelope violation clears
the returned candidates and sets `status` to `EnvelopeFailed`. Use
`EnvelopeFailurePolicy::Diagnostic` for sampler studies that should retain
accepted candidates together with the failure diagnostics.

A compact diagnostic executable exercises this API without writing event
records:

```sh
./build/lowmc_shape_unweighting_diagnostic \
  --trials 1000 \
  --calibration-trials 1000 \
  --momentum 100 \
  --safety-factor 3 \
  --probe-bins 11 \
  --failure-policy diagnostic
```

It prints one machine-readable key-value summary line with the status,
accepted candidate count, envelope-failure count, acceptance efficiency, probe
evaluation count, and related diagnostics.
