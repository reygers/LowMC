#pragma once

#include "lowmc/FourVector.h"
#include "lowmc/LowCurrent.h"
#include "lowmc/Photon.h"
#include "lowmc/PhotonPhaseSpace.h"

namespace lowmc {

struct WeightedPhotonProposal {
    FourVector momentum;
    double proposalDensity = 0.0;
};

struct WeightedPhotonSample {
    Photon photon;
    double proposalDensity = 0.0;
};

[[nodiscard]] WeightedPhotonSample makeWeightedPhotonSample(
    const LowCurrent& current,
    const PhotonPhaseSpace& phaseSpace,
    const WeightedPhotonProposal& proposal,
    double alpha = FineStructureConstant);

}  // namespace lowmc
