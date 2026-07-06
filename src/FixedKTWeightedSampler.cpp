#include "lowmc/FixedKTWeightedSampler.h"

#include <cmath>
#include <stdexcept>

namespace lowmc {

FourVector photonMomentumFromKTEtaPhi(double kT, double eta, double phi) {
    const double px = kT * std::cos(phi);
    const double py = kT * std::sin(phi);
    const double pz = kT * std::sinh(eta);
    const double energy = kT * std::cosh(eta);
    return FourVector{energy, px, py, pz};
}

FixedKTWeightedSample makeFixedKTUniformEtaPhiWeightedSample(
    const LowCurrent& current,
    const PhotonPhaseSpace& phaseSpace,
    double kT,
    double unitEta,
    double unitPhi,
    double alpha) {
    if (!validatePhotonPhaseSpace(phaseSpace).valid()) {
        throw std::invalid_argument("invalid photon phase space");
    }

    if (!std::isfinite(kT) || kT < phaseSpace.kTMin || kT > phaseSpace.kTMax) {
        throw std::invalid_argument("fixed kT must be finite and inside the phase-space range");
    }

    const UniformEtaPhiProposalPoint etaPhiProposal =
        sampleUniformEtaPhi(phaseSpace, unitEta, unitPhi);
    const FourVector momentum =
        photonMomentumFromKTEtaPhi(kT, etaPhiProposal.eta, etaPhiProposal.phi);
    const WeightedPhotonSample weightedPhoton = makeWeightedPhotonSample(
        current,
        phaseSpace,
        WeightedPhotonProposal{momentum, etaPhiProposal.density},
        alpha);

    return FixedKTWeightedSample{
        kT,
        etaPhiProposal,
        weightedPhoton,
    };
}

FixedKTAdaptiveWeightedSample makeFixedKTAdaptiveEtaPhiWeightedSample(
    const LowCurrent& current,
    const PhotonPhaseSpace& phaseSpace,
    double kT,
    double unitComponent,
    double unitEta,
    double unitPhi,
    const AdaptiveEtaPhiProposalConfig& config,
    double alpha) {
    if (!validatePhotonPhaseSpace(phaseSpace).valid()) {
        throw std::invalid_argument("invalid photon phase space");
    }

    if (!std::isfinite(kT) || kT < phaseSpace.kTMin || kT > phaseSpace.kTMax) {
        throw std::invalid_argument("fixed kT must be finite and inside the phase-space range");
    }

    const AdaptiveEtaPhiProposal proposal =
        makeAdaptiveEtaPhiProposal(phaseSpace, current.particles(), config);
    const AdaptiveEtaPhiProposalPoint etaPhiProposal =
        sampleAdaptiveEtaPhi(proposal, unitComponent, unitEta, unitPhi);
    const FourVector momentum =
        photonMomentumFromKTEtaPhi(kT, etaPhiProposal.eta, etaPhiProposal.phi);
    const WeightedPhotonSample weightedPhoton = makeWeightedPhotonSample(
        current,
        phaseSpace,
        WeightedPhotonProposal{momentum, etaPhiProposal.density},
        alpha);

    return FixedKTAdaptiveWeightedSample{
        kT,
        etaPhiProposal,
        weightedPhoton,
    };
}

FixedKTSmoothWeightedSample makeFixedKTSmoothEtaPhiWeightedSample(
    const LowCurrent& current,
    const PhotonPhaseSpace& phaseSpace,
    double kT,
    double unitComponent,
    double unitEta,
    double unitPhi,
    const SmoothEtaPhiProposalConfig& config,
    double alpha) {
    if (!validatePhotonPhaseSpace(phaseSpace).valid()) {
        throw std::invalid_argument("invalid photon phase space");
    }

    if (!std::isfinite(kT) || kT < phaseSpace.kTMin || kT > phaseSpace.kTMax) {
        throw std::invalid_argument("fixed kT must be finite and inside the phase-space range");
    }

    const SmoothEtaPhiProposal proposal =
        makeSmoothEtaPhiProposal(phaseSpace, current.particles(), config);
    const SmoothEtaPhiProposalPoint etaPhiProposal =
        sampleSmoothEtaPhi(proposal, unitComponent, unitEta, unitPhi);
    const FourVector momentum =
        photonMomentumFromKTEtaPhi(kT, etaPhiProposal.eta, etaPhiProposal.phi);
    const WeightedPhotonSample weightedPhoton = makeWeightedPhotonSample(
        current,
        phaseSpace,
        WeightedPhotonProposal{momentum, etaPhiProposal.density},
        alpha);

    return FixedKTSmoothWeightedSample{
        kT,
        etaPhiProposal,
        weightedPhoton,
    };
}

FixedKTGridAdaptiveWeightedSample makeFixedKTGridAdaptiveEtaPhiWeightedSample(
    const LowCurrent& current,
    const PhotonPhaseSpace& phaseSpace,
    double kT,
    double unitComponent,
    double unitCell,
    double unitEta,
    double unitPhi,
    const GridAdaptiveEtaPhiProposalConfig& config,
    double alpha) {
    if (!validatePhotonPhaseSpace(phaseSpace).valid()) {
        throw std::invalid_argument("invalid photon phase space");
    }

    if (!std::isfinite(kT) || kT < phaseSpace.kTMin || kT > phaseSpace.kTMax) {
        throw std::invalid_argument("fixed kT must be finite and inside the phase-space range");
    }

    const GridAdaptiveEtaPhiProposal proposal =
        makeGridAdaptiveEtaPhiProposal(current, phaseSpace, kT, config);
    const GridAdaptiveEtaPhiProposalPoint etaPhiProposal =
        sampleGridAdaptiveEtaPhi(proposal, unitComponent, unitCell, unitEta, unitPhi);
    const FourVector momentum =
        photonMomentumFromKTEtaPhi(kT, etaPhiProposal.eta, etaPhiProposal.phi);
    const WeightedPhotonSample weightedPhoton = makeWeightedPhotonSample(
        current,
        phaseSpace,
        WeightedPhotonProposal{momentum, etaPhiProposal.density},
        alpha);

    return FixedKTGridAdaptiveWeightedSample{
        kT,
        etaPhiProposal,
        weightedPhoton,
    };
}

}  // namespace lowmc
