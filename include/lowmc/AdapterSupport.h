#pragma once

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

#include "lowmc/ChargedParticle.h"
#include "lowmc/FourVector.h"
#include "lowmc/PDGCharge.h"

namespace lowmc {

enum class ParticleSelection {
    IncludeIncomingBeams,
    IncludeFinalStableCharged,
    IncludeUnstableIntermediateCharged,
    IncludeBeamRemnants,
};

struct ParticleSelectionConfig {
    bool includeIncomingBeams = true;
    bool includeFinalStableCharged = true;
    bool includeUnstableIntermediateCharged = false;
    bool includeBeamRemnants = false;

    [[nodiscard]] constexpr bool enabled(ParticleSelection selection) const {
        switch (selection) {
        case ParticleSelection::IncludeIncomingBeams:
            return includeIncomingBeams;
        case ParticleSelection::IncludeFinalStableCharged:
            return includeFinalStableCharged;
        case ParticleSelection::IncludeUnstableIntermediateCharged:
            return includeUnstableIntermediateCharged;
        case ParticleSelection::IncludeBeamRemnants:
            return includeBeamRemnants;
        }

        return false;
    }

    constexpr void set(ParticleSelection selection, bool value) {
        switch (selection) {
        case ParticleSelection::IncludeIncomingBeams:
            includeIncomingBeams = value;
            break;
        case ParticleSelection::IncludeFinalStableCharged:
            includeFinalStableCharged = value;
            break;
        case ParticleSelection::IncludeUnstableIntermediateCharged:
            includeUnstableIntermediateCharged = value;
            break;
        case ParticleSelection::IncludeBeamRemnants:
            includeBeamRemnants = value;
            break;
        }
    }
};

enum class RejectedParticleReason {
    Neutral,
    UnstableIntermediate,
    BeamRemnant,
    OutsideSelection,
    UnsupportedCharge,
    UnsupportedStatus,
    AmbiguousRecord,
    InvalidMomentum,
    Unknown,
};

struct RejectedChargedCandidate {
    RejectedParticleReason reason = RejectedParticleReason::Unknown;
    double charge = 0.0;
    std::string sourceLabel;
};

struct AdapterSelectionDiagnostic {
    std::size_t selectedIncomingCharged = 0;
    std::size_t selectedOutgoingCharged = 0;
    std::size_t ignoredNeutralCandidates = 0;
    double signedChargeSum = 0.0;
    double ordinaryChargeSum = 0.0;
    std::vector<RejectedChargedCandidate> rejectedChargedCandidates;

    [[nodiscard]] std::size_t selectedChargedCount() const {
        return selectedIncomingCharged + selectedOutgoingCharged;
    }

    [[nodiscard]] std::size_t rejectedChargedCount() const {
        return rejectedChargedCandidates.size();
    }
};

[[nodiscard]] inline AdapterSelectionDiagnostic summarizeAdapterSelection(
    const std::vector<ChargedParticle>& selectedParticles,
    std::vector<RejectedChargedCandidate> rejectedChargedCandidates = {}) {
    AdapterSelectionDiagnostic diagnostic;
    diagnostic.rejectedChargedCandidates = std::move(rejectedChargedCandidates);

    for (const ChargedParticle& particle : selectedParticles) {
        if (particle.inOutSign == InOutSign::Incoming) {
            ++diagnostic.selectedIncomingCharged;
        } else {
            ++diagnostic.selectedOutgoingCharged;
        }

        diagnostic.signedChargeSum +=
            static_cast<double>(signValue(particle.inOutSign)) * particle.charge;
        diagnostic.ordinaryChargeSum += particle.charge;
    }

    return diagnostic;
}

struct AdapterParticleCandidate {
    FourVector momentum;
    int pdgId = 0;
    bool isIncomingBeam = false;
    bool isFinalStable = false;
    bool isBeamRemnant = false;
    std::string sourceLabel;
};

struct AdapterParticleSelectionResult {
    std::vector<ChargedParticle> selectedParticles;
    AdapterSelectionDiagnostic diagnostic;
};

[[nodiscard]] inline RejectedChargedCandidate rejectedCandidate(
    RejectedParticleReason reason,
    double charge,
    const AdapterParticleCandidate& candidate) {
    return RejectedChargedCandidate{reason, charge, candidate.sourceLabel};
}

[[nodiscard]] inline AdapterParticleSelectionResult selectAdapterParticles(
    const std::vector<AdapterParticleCandidate>& candidates,
    const ParticleSelectionConfig& config = {}) {
    std::vector<ChargedParticle> selectedParticles;
    std::vector<RejectedChargedCandidate> rejectedChargedCandidates;
    std::size_t ignoredNeutralCandidates = 0;

    for (const AdapterParticleCandidate& candidate : candidates) {
        const PDGChargeLookupResult chargeResult = chargeTimesThreeFromPDG(candidate.pdgId);
        if (!chargeResult.known) {
            rejectedChargedCandidates.push_back(rejectedCandidate(
                RejectedParticleReason::UnsupportedCharge,
                0.0,
                candidate));
            continue;
        }

        if (chargeResult.chargeTimesThree == 0) {
            ++ignoredNeutralCandidates;
            continue;
        }

        const double charge = static_cast<double>(chargeResult.chargeTimesThree) / 3.0;
        if (candidate.isBeamRemnant && !config.includeBeamRemnants) {
            rejectedChargedCandidates.push_back(rejectedCandidate(
                RejectedParticleReason::BeamRemnant,
                charge,
                candidate));
            continue;
        }

        const bool selectIncoming =
            candidate.isIncomingBeam && config.includeIncomingBeams;
        const bool selectOutgoing =
            (candidate.isFinalStable && config.includeFinalStableCharged) ||
            (!candidate.isIncomingBeam &&
             !candidate.isFinalStable &&
             !candidate.isBeamRemnant &&
             config.includeUnstableIntermediateCharged);

        if (selectIncoming && selectOutgoing) {
            rejectedChargedCandidates.push_back(rejectedCandidate(
                RejectedParticleReason::AmbiguousRecord,
                charge,
                candidate));
            continue;
        }

        if (selectIncoming || selectOutgoing) {
            selectedParticles.push_back(ChargedParticle{
                candidate.momentum,
                charge,
                selectIncoming ? InOutSign::Incoming : InOutSign::Outgoing,
            });
            continue;
        }

        RejectedParticleReason reason = RejectedParticleReason::OutsideSelection;
        if (!candidate.isIncomingBeam && !candidate.isFinalStable && !candidate.isBeamRemnant) {
            reason = RejectedParticleReason::UnstableIntermediate;
        } else if (candidate.isBeamRemnant) {
            reason = RejectedParticleReason::BeamRemnant;
        }

        rejectedChargedCandidates.push_back(rejectedCandidate(reason, charge, candidate));
    }

    AdapterSelectionDiagnostic diagnostic =
        summarizeAdapterSelection(selectedParticles, std::move(rejectedChargedCandidates));
    diagnostic.ignoredNeutralCandidates = ignoredNeutralCandidates;

    return AdapterParticleSelectionResult{std::move(selectedParticles), std::move(diagnostic)};
}

}  // namespace lowmc
