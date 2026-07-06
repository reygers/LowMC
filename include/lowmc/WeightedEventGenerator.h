#pragma once

#include <cstdint>
#include <vector>

#include "lowmc/ChargedParticle.h"
#include "lowmc/GridEtaPhiProposal.h"
#include "lowmc/LowCurrent.h"
#include "lowmc/Photon.h"
#include "lowmc/PhotonPhaseSpace.h"
#include "lowmc/SamplerDiagnostics.h"

namespace lowmc {

struct WeightedEventGeneratorConfig {
    PhotonPhaseSpace phaseSpace;
    int samples = 1;
    std::uint64_t seed = 12345;
    GridAdaptiveEtaPhiProposalConfig gridConfig = tunedMassGridAdaptiveEtaPhiConfig();
    int oversamplingFactor = 1;
};

struct WeightedEventPhotonSample {
    Photon photon;
    double kT = 0.0;
    double eta = 0.0;
    double phi = 0.0;
    double proposalDensity = 0.0;
};

struct WeightedEventGenerationResult {
    std::vector<WeightedEventPhotonSample> samples;
    std::vector<Photon> photons;
    WeightSummary weightSummary;
    ChargeConservationDiagnostic chargeConservation;
    double gridBuildKT = 0.0;
};

[[nodiscard]] WeightedEventGenerationResult generateWeightedPhotons(
    const std::vector<ChargedParticle>& selectedParticles,
    const WeightedEventGeneratorConfig& config,
    double alpha = FineStructureConstant);

}  // namespace lowmc
