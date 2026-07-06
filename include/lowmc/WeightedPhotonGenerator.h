#pragma once

#include "lowmc/EtaPhiProposal.h"
#include "lowmc/FixedKTWeightedSampler.h"
#include "lowmc/GridEtaPhiProposal.h"
#include "lowmc/KTProposal.h"
#include "lowmc/LowCurrent.h"
#include "lowmc/PhotonPhaseSpace.h"
#include "lowmc/WeightedPhotonSampler.h"

namespace lowmc {

struct MinimalWeightedPhotonSample {
    LogUniformKTProposalPoint kTProposal;
    UniformEtaPhiProposalPoint etaPhiProposal;
    double proposalDensity = 0.0;
    WeightedPhotonSample weightedPhoton;
};

struct GridAdaptiveWeightedPhotonSample {
    LogUniformKTProposalPoint kTProposal;
    GridAdaptiveEtaPhiProposalPoint etaPhiProposal;
    double proposalDensity = 0.0;
    WeightedPhotonSample weightedPhoton;
};

[[nodiscard]] MinimalWeightedPhotonSample makeMinimalWeightedPhotonSample(
    const LowCurrent& current,
    const PhotonPhaseSpace& phaseSpace,
    double unitKT,
    double unitEta,
    double unitPhi,
    double alpha = FineStructureConstant);

[[nodiscard]] GridAdaptiveWeightedPhotonSample makeGridAdaptiveWeightedPhotonSample(
    const LowCurrent& current,
    const PhotonPhaseSpace& phaseSpace,
    const GridAdaptiveEtaPhiProposal& etaPhiProposal,
    double unitKT,
    double unitComponent,
    double unitCell,
    double unitEta,
    double unitPhi,
    double alpha = FineStructureConstant);

}  // namespace lowmc
