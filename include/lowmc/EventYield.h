#pragma once

#include <cstdint>
#include <vector>

#include "lowmc/ChargedParticle.h"
#include "lowmc/GridEtaPhiProposal.h"
#include "lowmc/LowCurrent.h"
#include "lowmc/PhotonPhaseSpace.h"

namespace lowmc {

struct EventYieldIntegrationConfig {
    int logKTBins = 24;
    int etaBins = 80;
    int phiBins = 96;
};

struct AdaptiveEventYieldIntegrationConfig {
    int logKTBins = 24;
    int etaSubBins = 2;
    int phiSubBins = 2;
    GridAdaptiveEtaPhiProposalConfig gridConfig = tunedMassGridAdaptiveEtaPhiConfig();
};

struct EventYieldEstimate {
    double lambdaEvent = 0.0;
    double coarseLambdaEvent = 0.0;
    double absoluteDifference = 0.0;
    double relativeDifference = 0.0;
    int logKTBins = 0;
    int etaBins = 0;
    int phiBins = 0;
    int evaluations = 0;
    int gridCells = 0;
    const char* method = "log-kt-eta-phi-midpoint";
    const char* status = "ok";
};

struct PhotonMultiplicityConfig {
    double oversamplingFactor = 1.0;
    std::uint64_t seed = 12345;
};

struct PhotonMultiplicityDiagnostic {
    double lambdaEvent = 0.0;
    double poissonMean = 0.0;
    int sampledPhotons = 0;
    double sampleWeight = 1.0;
    double expectedWeightedCount = 0.0;
    double oversamplingFactor = 1.0;
    std::uint64_t seed = 0;
};

[[nodiscard]] EventYieldEstimate estimateEventYield(
    const LowCurrent& current,
    const PhotonPhaseSpace& phaseSpace,
    const EventYieldIntegrationConfig& config = {},
    double alpha = FineStructureConstant);

[[nodiscard]] EventYieldEstimate estimateEventYield(
    const std::vector<ChargedParticle>& selectedParticles,
    const PhotonPhaseSpace& phaseSpace,
    const EventYieldIntegrationConfig& config = {},
    double alpha = FineStructureConstant);

[[nodiscard]] EventYieldEstimate estimateAdaptiveEventYield(
    const LowCurrent& current,
    const PhotonPhaseSpace& phaseSpace,
    const AdaptiveEventYieldIntegrationConfig& config = {},
    double alpha = FineStructureConstant);

[[nodiscard]] EventYieldEstimate estimateAdaptiveEventYield(
    const std::vector<ChargedParticle>& selectedParticles,
    const PhotonPhaseSpace& phaseSpace,
    const AdaptiveEventYieldIntegrationConfig& config = {},
    double alpha = FineStructureConstant);

[[nodiscard]] EventYieldEstimate estimateFactorizedAdaptiveEventYield(
    const LowCurrent& current,
    const PhotonPhaseSpace& phaseSpace,
    const AdaptiveEventYieldIntegrationConfig& config = {},
    double alpha = FineStructureConstant);

[[nodiscard]] EventYieldEstimate estimateFactorizedAdaptiveEventYield(
    const std::vector<ChargedParticle>& selectedParticles,
    const PhotonPhaseSpace& phaseSpace,
    const AdaptiveEventYieldIntegrationConfig& config = {},
    double alpha = FineStructureConstant);

[[nodiscard]] PhotonMultiplicityDiagnostic samplePhotonMultiplicity(
    double lambdaEvent,
    const PhotonMultiplicityConfig& config = {});

}  // namespace lowmc
