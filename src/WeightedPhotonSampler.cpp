#include "lowmc/WeightedPhotonSampler.h"

#include <cmath>
#include <stdexcept>

namespace lowmc {

WeightedPhotonSample makeWeightedPhotonSample(
    const LowCurrent& current,
    const PhotonPhaseSpace& phaseSpace,
    const WeightedPhotonProposal& proposal,
    double alpha) {
    if (!validatePhotonPhaseSpace(phaseSpace).valid()) {
        throw std::invalid_argument("invalid photon phase space");
    }

    if (!std::isfinite(proposal.proposalDensity) || proposal.proposalDensity <= 0.0) {
        throw std::invalid_argument("proposal density must be positive and finite");
    }

    const double density = lowPhotonDensity(current, proposal.momentum, alpha);
    return WeightedPhotonSample{
        Photon{proposal.momentum, density / proposal.proposalDensity},
        proposal.proposalDensity,
    };
}

}  // namespace lowmc
