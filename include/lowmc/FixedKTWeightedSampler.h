#pragma once

#include "lowmc/EtaPhiProposal.h"
#include "lowmc/FourVector.h"
#include "lowmc/GridEtaPhiProposal.h"
#include "lowmc/LegCenteredEtaPhiProposal.h"
#include "lowmc/LowCurrent.h"
#include "lowmc/PhotonPhaseSpace.h"
#include "lowmc/WeightedPhotonSampler.h"

namespace lowmc {

struct FixedKTWeightedSample {
    double kT = 0.0;
    UniformEtaPhiProposalPoint etaPhiProposal;
    WeightedPhotonSample weightedPhoton;
};

struct FixedKTAdaptiveWeightedSample {
    double kT = 0.0;
    AdaptiveEtaPhiProposalPoint etaPhiProposal;
    WeightedPhotonSample weightedPhoton;
};

struct FixedKTSmoothWeightedSample {
    double kT = 0.0;
    SmoothEtaPhiProposalPoint etaPhiProposal;
    WeightedPhotonSample weightedPhoton;
};

struct FixedKTGridAdaptiveWeightedSample {
    double kT = 0.0;
    GridAdaptiveEtaPhiProposalPoint etaPhiProposal;
    WeightedPhotonSample weightedPhoton;
};

struct FixedKTLegCenteredWeightedSample {
    double kT = 0.0;
    LegCenteredEtaPhiProposalPoint etaPhiProposal;
    WeightedPhotonSample weightedPhoton;
};

[[nodiscard]] FourVector photonMomentumFromKTEtaPhi(double kT, double eta, double phi);

[[nodiscard]] FixedKTWeightedSample makeFixedKTUniformEtaPhiWeightedSample(
    const LowCurrent& current,
    const PhotonPhaseSpace& phaseSpace,
    double kT,
    double unitEta,
    double unitPhi,
    double alpha = FineStructureConstant);

[[nodiscard]] FixedKTAdaptiveWeightedSample makeFixedKTAdaptiveEtaPhiWeightedSample(
    const LowCurrent& current,
    const PhotonPhaseSpace& phaseSpace,
    double kT,
    double unitComponent,
    double unitEta,
    double unitPhi,
    const AdaptiveEtaPhiProposalConfig& config = {},
    double alpha = FineStructureConstant);

[[nodiscard]] FixedKTSmoothWeightedSample makeFixedKTSmoothEtaPhiWeightedSample(
    const LowCurrent& current,
    const PhotonPhaseSpace& phaseSpace,
    double kT,
    double unitComponent,
    double unitEta,
    double unitPhi,
    const SmoothEtaPhiProposalConfig& config = {},
    double alpha = FineStructureConstant);

[[nodiscard]] FixedKTGridAdaptiveWeightedSample makeFixedKTGridAdaptiveEtaPhiWeightedSample(
    const LowCurrent& current,
    const PhotonPhaseSpace& phaseSpace,
    double kT,
    double unitComponent,
    double unitCell,
    double unitEta,
    double unitPhi,
    const GridAdaptiveEtaPhiProposalConfig& config = {},
    double alpha = FineStructureConstant);

[[nodiscard]] FixedKTLegCenteredWeightedSample makeFixedKTLegCenteredEtaPhiWeightedSample(
    const LowCurrent& current,
    double kT,
    const LegCenteredEtaPhiProposal& proposal,
    double unitComponent,
    double unitPolar,
    double unitAzimuth,
    double alpha = FineStructureConstant);

}  // namespace lowmc
