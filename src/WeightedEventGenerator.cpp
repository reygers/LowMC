#include "lowmc/WeightedEventGenerator.h"

#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>

#include "lowmc/WeightedPhotonGenerator.h"

namespace lowmc {

WeightedEventGenerationResult generateWeightedPhotons(
    const std::vector<ChargedParticle>& selectedParticles,
    const WeightedEventGeneratorConfig& config,
    double alpha) {
    if (!validatePhotonPhaseSpace(config.phaseSpace).valid()) {
        throw std::invalid_argument("invalid photon phase space");
    }
    if (config.samples < 0) {
        throw std::invalid_argument("weighted event generator sample count must be non-negative");
    }
    if (config.oversamplingFactor <= 0) {
        throw std::invalid_argument("weighted event generator oversampling factor must be positive");
    }
    if (config.samples > 0 &&
        config.samples > std::numeric_limits<int>::max() / config.oversamplingFactor) {
        throw std::invalid_argument("weighted event generator oversampled sample count is too large");
    }
    validateGridAdaptiveEtaPhiProposalConfig(config.gridConfig);

    const int totalSamples = config.samples * config.oversamplingFactor;
    const double inverseOversamplingFactor =
        1.0 / static_cast<double>(config.oversamplingFactor);

    const LowCurrent current(selectedParticles);
    WeightedEventGenerationResult result;
    result.chargeConservation = current.chargeConservation();
    result.gridBuildKT = std::sqrt(config.phaseSpace.kTMin * config.phaseSpace.kTMax);
    result.samples.reserve(static_cast<std::size_t>(totalSamples));
    result.photons.reserve(static_cast<std::size_t>(totalSamples));

    if (totalSamples == 0) {
        result.weightSummary = summarizePhotons(result.photons);
        return result;
    }

    const GridAdaptiveEtaPhiProposal gridProposal =
        makeGridAdaptiveEtaPhiProposal(
            current,
            config.phaseSpace,
            result.gridBuildKT,
            config.gridConfig);

    std::mt19937_64 rng(config.seed);
    for (int index = 0; index < totalSamples; ++index) {
        const double unitKT = std::generate_canonical<double, 53>(rng);
        const double unitEta = std::generate_canonical<double, 53>(rng);
        const double unitPhi = std::generate_canonical<double, 53>(rng);
        const double unitComponent = std::generate_canonical<double, 53>(rng);
        const double unitCell = std::generate_canonical<double, 53>(rng);

        const GridAdaptiveWeightedPhotonSample sample =
            makeGridAdaptiveWeightedPhotonSample(
                current,
                config.phaseSpace,
                gridProposal,
                unitKT,
                unitComponent,
                unitCell,
                unitEta,
                unitPhi,
                alpha);
        Photon photon = sample.weightedPhoton.photon;
        photon.weight *= inverseOversamplingFactor;
        result.samples.push_back(WeightedEventPhotonSample{
            photon,
            sample.kTProposal.kT,
            sample.etaPhiProposal.eta,
            sample.etaPhiProposal.phi,
            sample.proposalDensity,
        });
        result.photons.push_back(photon);
    }

    result.weightSummary = summarizePhotons(result.photons);
    return result;
}

}  // namespace lowmc
