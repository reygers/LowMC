#include "lowmc/LowMC.h"

#include <cmath>
#include <iostream>
#include <vector>

int main() {
    const double pionMass = 0.13957039;
    const double momentum = 10.0;
    const double energy = std::sqrt(momentum * momentum + pionMass * pionMass);
    const std::vector<lowmc::ChargedParticle> particles = {
        {lowmc::FourVector{energy, 0.0, momentum, 0.0}, +1.0, lowmc::InOutSign::Outgoing},
        {lowmc::FourVector{energy, 0.0, -momentum, 0.0}, -1.0, lowmc::InOutSign::Outgoing},
    };

    lowmc::ShapeUnweightingConfig config;
    config.phaseSpace = lowmc::PhotonPhaseSpace{
        0.002,
        0.2,
        -5.0,
        5.0,
        0.0,
        2.0 * std::acos(-1.0),
    };
    config.trials = 5000;
    config.seed = 12345;
    config.calibrationTrials = 1000;
    config.probePolicy = lowmc::ShapeUnweightingProbePolicy::BoostedCells;
    config.boostedProbeBins = 21;
    config.boostedCellSafetyFactor = 5.0;
    config.boostedCellNeighborScale = 0.0;
    config.targetedRefinementLevels = 5;
    config.targetedRefinementNeighborScale = 0.5;
    config.failurePolicy = lowmc::EnvelopeFailurePolicy::Strict;

    const lowmc::ShapeUnweightingResult result =
        lowmc::generateShapeUnweightedPhotonCandidates(particles, config);

    std::cout << "Example: generate unweighted Low-shape candidates\n";
    std::cout << "status_ok="
              << (result.status == lowmc::ShapeUnweightingStatus::Ok) << '\n';
    std::cout << "accepted_candidates=" << result.candidates.size() << '\n';
    std::cout << "envelope_failures="
              << result.envelope.failedEnvelopeCount << '\n';
    std::cout << "max_trial_bound_ratio="
              << result.envelope.maxTrialBoundRatio << '\n';
    std::cout << "acceptance_efficiency="
              << result.envelope.observedAcceptanceEfficiency << '\n';
    std::cout << "note=shape_candidates_are_not_a_physical_multiplicity\n";

    const bool valid =
        result.status == lowmc::ShapeUnweightingStatus::Ok &&
        result.envelope.failedEnvelopeCount == 0 &&
        !result.candidates.empty();
    return valid ? 0 : 1;
}
