#pragma once

#include <vector>

#include "lowmc/LowCurrent.h"
#include "lowmc/PhotonPhaseSpace.h"

namespace lowmc {

enum class GridAdaptiveRefinementScore {
    TargetMass,
    Variance,
};

struct GridAdaptiveEtaPhiProposalConfig {
    int etaBins = 16;
    int phiBins = 16;
    int refinementLevels = 2;
    double refineFraction = 0.1;
    double uniformFraction = 0.2;
    GridAdaptiveRefinementScore refinementScore = GridAdaptiveRefinementScore::TargetMass;
    int varianceProbeBins = 3;
};

struct GridAdaptiveEtaPhiCell {
    double etaMin = 0.0;
    double etaMax = 0.0;
    double phiMin = 0.0;
    double phiMax = 0.0;
    double targetMass = 0.0;
    double refinementScore = 0.0;
    double mixtureWeight = 0.0;
};

struct GridAdaptiveEtaPhiProposal {
    PhotonPhaseSpace phaseSpace;
    GridAdaptiveEtaPhiProposalConfig config;
    std::vector<GridAdaptiveEtaPhiCell> cells;
};

struct GridAdaptiveEtaPhiProposalPoint {
    double eta = 0.0;
    double phi = 0.0;
    double density = 0.0;
};

void validateGridAdaptiveEtaPhiProposalConfig(
    const GridAdaptiveEtaPhiProposalConfig& config);

[[nodiscard]] constexpr GridAdaptiveEtaPhiProposalConfig tunedMassGridAdaptiveEtaPhiConfig() {
    return GridAdaptiveEtaPhiProposalConfig{
        32,
        32,
        3,
        0.1,
        0.1,
        GridAdaptiveRefinementScore::TargetMass,
        3,
    };
}

[[nodiscard]] GridAdaptiveEtaPhiProposal makeGridAdaptiveEtaPhiProposal(
    const LowCurrent& current,
    const PhotonPhaseSpace& phaseSpace,
    double kT,
    const GridAdaptiveEtaPhiProposalConfig& config = {});

[[nodiscard]] double gridAdaptiveEtaPhiDensity(
    const GridAdaptiveEtaPhiProposal& proposal,
    double eta,
    double phi);

[[nodiscard]] GridAdaptiveEtaPhiProposalPoint sampleGridAdaptiveEtaPhi(
    const GridAdaptiveEtaPhiProposal& proposal,
    double unitComponent,
    double unitCell,
    double unitEta,
    double unitPhi);

}  // namespace lowmc
