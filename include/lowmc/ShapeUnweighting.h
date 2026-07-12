#pragma once

#include <cstdint>
#include <vector>

#include "lowmc/ChargedParticle.h"
#include "lowmc/GridEtaPhiProposal.h"
#include "lowmc/LowCurrent.h"
#include "lowmc/Photon.h"
#include "lowmc/PhotonPhaseSpace.h"

namespace lowmc {

enum class ShapeUnweightingProbeRule {
    Interior,
    Edge,
};

enum class ShapeUnweightingProbePolicy {
    Fixed,
    BoostedCells,
};

enum class EnvelopeFailurePolicy {
    Strict,
    Diagnostic,
};

enum class ShapeUnweightingStatus {
    Ok,
    EnvelopeFailed,
    InvalidSelection,
    NoTrials,
};

struct ShapeUnweightingConfig {
    PhotonPhaseSpace phaseSpace;
    int trials = 1;
    std::uint64_t seed = 12345;
    GridAdaptiveEtaPhiProposalConfig gridConfig = tunedMassGridAdaptiveEtaPhiConfig();
    int calibrationTrials = 1000;
    int probeBins = 11;
    int baseProbeBins = 0;
    bool auditFullBaseProbes = false;
    bool targetedBaseProbes = false;
    double targetedBaseProbeNeighborScale = 0.0;
    bool evidenceTargetedBaseProbes = false;
    double evidenceTargetedBaseProbeNeighborScale = 0.0;
    ShapeUnweightingProbeRule probeRule = ShapeUnweightingProbeRule::Interior;
    ShapeUnweightingProbePolicy probePolicy = ShapeUnweightingProbePolicy::Fixed;
    int boostedProbeBins = 21;
    double boostedCellSafetyFactor = 20.0;
    double boostedCellNeighborScale = 1.0;
    int targetedRefinementLevels = 0;
    double targetedRefinementNeighborScale = 0.0;
    double safetyFactor = 3.0;
    EnvelopeFailurePolicy failurePolicy = EnvelopeFailurePolicy::Strict;
};

struct ShapeUnweightedPhotonCandidate {
    Photon photon;
    double kT = 0.0;
    double eta = 0.0;
    double phi = 0.0;
};

struct ShapeUnweightingEnvelopeDiagnostics {
    int calibrationCount = 0;
    int trialCount = 0;
    int acceptedCount = 0;
    int failedEnvelopeCount = 0;
    int fallbackCount = 0;
    int suspectCellCount = 0;
    int promotedCellCount = 0;
    std::size_t gridCellCount = 0;
    std::size_t probeEvaluations = 0;
    double gridSeconds = 0.0;
    double baseProbeSeconds = 0.0;
    double calibrationSeconds = 0.0;
    double probePlanSeconds = 0.0;
    double boostedProbeSeconds = 0.0;
    double envelopeBuildSeconds = 0.0;
    double trialSeconds = 0.0;
    double finalizationSeconds = 0.0;
    double maxTrialBoundRatio = 0.0;
    double observedAcceptanceEfficiency = 0.0;
    double expectedAcceptanceEfficiency = 0.0;
    double envelopeIntegral = 0.0;
    double safetyFactor = 0.0;
    double boostedCellSafetyFactor = 0.0;
    double boostedCellNeighborScale = 0.0;
    int targetedRefinementLevels = 0;
    double targetedRefinementNeighborScale = 0.0;
    int probeBins = 0;
    int baseProbeBins = 0;
    int auditCellCount = 0;
    int auditUnderBoundedCellCount = 0;
    double auditMaxFullOverBaseProbeRatio = 0.0;
    int auditRatioAbove110CellCount = 0;
    int auditRatioAbove150CellCount = 0;
    int auditRatioAbove200CellCount = 0;
    int fullBaseProbeCellCount = 0;
    double targetedBaseProbeNeighborScale = 0.0;
    int evidenceBaseProbeCellCount = 0;
    double evidenceTargetedBaseProbeNeighborScale = 0.0;
    int boostedProbeBins = 0;
    ShapeUnweightingProbeRule probeRule = ShapeUnweightingProbeRule::Interior;
    ShapeUnweightingProbePolicy probePolicy = ShapeUnweightingProbePolicy::Fixed;
};

struct ShapeUnweightingResult {
    std::vector<ShapeUnweightedPhotonCandidate> candidates;
    ShapeUnweightingEnvelopeDiagnostics envelope;
    ChargeConservationDiagnostic chargeConservation;
    double gridBuildKT = 0.0;
    ShapeUnweightingStatus status = ShapeUnweightingStatus::Ok;
};

struct ShapeUnweightingIndexConsistencyDiagnostic {
    bool consistent = false;
    std::size_t gridCellCount = 0;
    std::size_t pointLookupChecks = 0;
    std::size_t envelopeSamplingChecks = 0;
    std::size_t knownCellProbeWeightChecks = 0;
};

[[nodiscard]] ShapeUnweightingResult generateShapeUnweightedPhotonCandidates(
    const std::vector<ChargedParticle>& selectedParticles,
    const ShapeUnweightingConfig& config,
    double alpha = FineStructureConstant);

[[nodiscard]] ShapeUnweightingIndexConsistencyDiagnostic
validateShapeUnweightingIndexConsistency(
    const std::vector<ChargedParticle>& selectedParticles,
    const ShapeUnweightingConfig& config,
    double alpha = FineStructureConstant);

}  // namespace lowmc
