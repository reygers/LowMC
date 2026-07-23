#pragma once

#include <vector>

#include "lowmc/ChargedParticle.h"
#include "lowmc/PhotonPhaseSpace.h"

namespace lowmc {

// Diagnostic analytic proposal inspired by the leg-centred crude kernels used
// in soft-photon afterburners. Each component is proportional to
// 1 / (1 - beta_i cos(theta_i_gamma)) on the full photon sphere. The returned
// density is with respect to d eta d phi.
struct LegCenteredEtaPhiProposalConfig {
    double uniformFraction = 0.2;
    double minimumOneMinusBeta = 1.0e-12;
};

struct LegCenteredEtaPhiProposalComponent {
    double nx = 0.0;
    double ny = 0.0;
    double nz = 0.0;
    double beta = 0.0;
    double sphereNormalization = 0.0;
    double mixtureWeight = 0.0;
};

struct LegCenteredEtaPhiProposal {
    PhotonPhaseSpace phaseSpace;
    LegCenteredEtaPhiProposalConfig config;
    std::vector<LegCenteredEtaPhiProposalComponent> components;
};

struct LegCenteredEtaPhiProposalPoint {
    double eta = 0.0;
    double phi = 0.0;
    double density = 0.0;
    bool insideAcceptance = false;
};

void validateLegCenteredEtaPhiProposalConfig(
    const LegCenteredEtaPhiProposalConfig& config);

[[nodiscard]] LegCenteredEtaPhiProposal makeLegCenteredEtaPhiProposal(
    const PhotonPhaseSpace& phaseSpace,
    const std::vector<ChargedParticle>& particles,
    const LegCenteredEtaPhiProposalConfig& config = {});

[[nodiscard]] double legCenteredEtaPhiDensity(
    const LegCenteredEtaPhiProposal& proposal,
    double eta,
    double phi);

[[nodiscard]] LegCenteredEtaPhiProposalPoint sampleLegCenteredEtaPhi(
    const LegCenteredEtaPhiProposal& proposal,
    double unitComponent,
    double unitPolar,
    double unitAzimuth);

}  // namespace lowmc
