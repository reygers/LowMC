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
    double reducedProbeCompensationFactor = 1.0;
    bool auditFullBaseProbes = false;
    double auditSevereRatioThreshold = 1.1;
    int auditSevereCellReportMaxRows = 0;
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
    bool auditAnalyticBounds = false;
};

struct ShapeUnweightedPhotonCandidate {
    Photon photon;
    double kT = 0.0;
    double eta = 0.0;
    double phi = 0.0;
};

struct ShapeUnweightingAuditCellDiagnostic {
    std::size_t cellIndex = 0;
    double etaCenter = 0.0;
    double phiCenter = 0.0;
    double etaWidth = 0.0;
    double phiWidth = 0.0;
    double area = 0.0;
    double proposalMass = 0.0;
    double reducedMaxWeight = 0.0;
    double fullMaxWeight = 0.0;
    double fullOverReducedRatio = 0.0;
    double calibrationMaxWeight = 0.0;
    double calibrationOverReducedRatio = 0.0;
    bool calibrationSuspect = false;
    bool chargedNear = false;
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
    double reducedProbeCompensationFactor = 1.0;
    int reducedProbeCompensatedCellCount = 0;
    int auditCellCount = 0;
    int auditUnderBoundedCellCount = 0;
    double auditMaxFullOverBaseProbeRatio = 0.0;
    int auditRatioAbove110CellCount = 0;
    int auditRatioAbove150CellCount = 0;
    int auditRatioAbove200CellCount = 0;
    double auditSeconds = 0.0;
    double auditSevereRatioThreshold = 1.1;
    int auditSevereCellCount = 0;
    int auditSevereCalibrationSuspectCellCount = 0;
    int auditSevereChargedNearCellCount = 0;
    double auditSevereProposalMassSum = 0.0;
    double auditSevereMeanCellArea = 0.0;
    double auditSevereMaxCellArea = 0.0;
    double auditSevereMeanEtaWidth = 0.0;
    double auditSevereMeanPhiWidth = 0.0;
    double auditSevereMaxCalibrationOverBaseProbeRatio = 0.0;
    bool analyticBoundAuditEnabled = false;
    bool analyticBoundAvailable = false;
    int analyticBoundInvalidParticleCount = 0;
    int analyticBoundCellCount = 0;
    int analyticBoundFallbackCellCount = 0;
    int analyticBoundObservedViolationCount = 0;
    int analyticBoundTrialViolationCount = 0;
    std::size_t analyticBoundPairEvaluations = 0;
    double analyticBoundSeconds = 0.0;
    double analyticGlobalAngularBound = 0.0;
    double analyticGlobalWeightBound = 0.0;
    double analyticCellEnvelopeIntegral = 0.0;
    double analyticCellEnvelopeOverCurrent = 0.0;
    double analyticExpectedAcceptanceEfficiency = 0.0;
    double analyticMaxObservedOverBound = 0.0;
    double analyticMaxTrialOverBound = 0.0;
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
    std::vector<ShapeUnweightingAuditCellDiagnostic> auditCells;
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
