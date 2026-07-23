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
    const lowmc::PhotonPhaseSpace phaseSpace{
        0.002,
        0.2,
        -5.0,
        5.0,
        0.0,
        2.0 * std::acos(-1.0),
    };

    lowmc::AdaptiveEventYieldIntegrationConfig integration;
    integration.gridConfig.refinementLevels = 7;
    const lowmc::EventYieldEstimate yield =
        lowmc::estimateAdaptiveEventYield(particles, phaseSpace, integration);

    // This is one physical Poisson draw: a result of zero is entirely normal
    // when lambda_event is below one.
    const lowmc::PhotonMultiplicityDiagnostic multiplicity =
        lowmc::samplePhotonMultiplicity(
            yield.lambdaEvent,
            lowmc::PhotonMultiplicityConfig{1.0, 12345});

    std::cout << "Example: integrate the event yield and sample multiplicity\n";
    std::cout << "integration_status=" << yield.status << '\n';
    std::cout << "lambda_event=" << yield.lambdaEvent << '\n';
    std::cout << "coarse_fine_relative_difference="
              << yield.relativeDifference << '\n';
    std::cout << "sampled_photons=" << multiplicity.sampledPhotons << '\n';
    std::cout << "sample_weight=" << multiplicity.sampleWeight << '\n';

    return std::isfinite(yield.lambdaEvent) && yield.lambdaEvent >= 0.0 ? 0 : 1;
}
