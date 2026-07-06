#include "lowmc/WeightedPhotonGenerator.h"

#include <stdexcept>

namespace lowmc {

MinimalWeightedPhotonSample makeMinimalWeightedPhotonSample(
    const LowCurrent& current,
    const PhotonPhaseSpace& phaseSpace,
    double unitKT,
    double unitEta,
    double unitPhi,
    double alpha) {
    const LogUniformKTProposalPoint kTProposal = sampleLogUniformKT(phaseSpace, unitKT);
    const UniformEtaPhiProposalPoint etaPhiProposal =
        sampleUniformEtaPhi(phaseSpace, unitEta, unitPhi);
    const FourVector momentum =
        photonMomentumFromKTEtaPhi(kTProposal.kT, etaPhiProposal.eta, etaPhiProposal.phi);
    const double proposalDensity =
        kTProposal.densityPhysicalMeasure * etaPhiProposal.density;
    const WeightedPhotonSample weightedPhoton = makeWeightedPhotonSample(
        current,
        phaseSpace,
        WeightedPhotonProposal{momentum, proposalDensity},
        alpha);

    return MinimalWeightedPhotonSample{
        kTProposal,
        etaPhiProposal,
        proposalDensity,
        weightedPhoton,
    };
}

GridAdaptiveWeightedPhotonSample makeGridAdaptiveWeightedPhotonSample(
    const LowCurrent& current,
    const PhotonPhaseSpace& phaseSpace,
    const GridAdaptiveEtaPhiProposal& gridProposal,
    double unitKT,
    double unitComponent,
    double unitCell,
    double unitEta,
    double unitPhi,
    double alpha) {
    if (gridProposal.phaseSpace.kTMin != phaseSpace.kTMin ||
        gridProposal.phaseSpace.kTMax != phaseSpace.kTMax ||
        gridProposal.phaseSpace.etaMin != phaseSpace.etaMin ||
        gridProposal.phaseSpace.etaMax != phaseSpace.etaMax ||
        gridProposal.phaseSpace.phiMin != phaseSpace.phiMin ||
        gridProposal.phaseSpace.phiMax != phaseSpace.phiMax) {
        throw std::invalid_argument("grid eta-phi proposal phase space does not match");
    }

    const LogUniformKTProposalPoint kTProposal = sampleLogUniformKT(phaseSpace, unitKT);
    const GridAdaptiveEtaPhiProposalPoint etaPhiProposal =
        sampleGridAdaptiveEtaPhi(
            gridProposal,
            unitComponent,
            unitCell,
            unitEta,
            unitPhi);
    const FourVector momentum =
        photonMomentumFromKTEtaPhi(kTProposal.kT, etaPhiProposal.eta, etaPhiProposal.phi);
    const double proposalDensity =
        kTProposal.densityPhysicalMeasure * etaPhiProposal.density;
    const WeightedPhotonSample weightedPhoton = makeWeightedPhotonSample(
        current,
        phaseSpace,
        WeightedPhotonProposal{momentum, proposalDensity},
        alpha);

    return GridAdaptiveWeightedPhotonSample{
        kTProposal,
        etaPhiProposal,
        proposalDensity,
        weightedPhoton,
    };
}

}  // namespace lowmc
