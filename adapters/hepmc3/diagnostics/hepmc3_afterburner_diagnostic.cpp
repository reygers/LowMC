#include "lowmc/LowMC.h"
#include "lowmc/hepmc3/HepMC3Adapter.h"

#include <HepMC3/GenEvent.h>
#include <HepMC3/ReaderAscii.h>
#include <HepMC3/Units.h>
#include <HepMC3/WriterAscii.h>

#include <algorithm>
#include <chrono>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <cstdio>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

enum class FailurePolicy {
    Skip,
    Fail,
};

struct Options {
    std::string inputPath;
    std::string outputPath;
    int maxEvents = std::numeric_limits<int>::max();
    double enhancementFactor = 1.0;
    std::uint64_t multiplicitySeed = 12345;
    std::uint64_t shapeSeed = 67890;
    int minShapeTrials = 10000;
    int shapeTrialsPerPhoton = 2000;
    int shapeBaseProbeBins = 11;
    bool auditFullBaseProbes = false;
    bool targetedBaseProbes = false;
    double targetedBaseProbeNeighborScale = 0.0;
    bool evidenceTargetedBaseProbes = false;
    double evidenceTargetedBaseProbeNeighborScale = 0.0;
    bool disableBoostedCellProbes = false;
    double kTMin = 0.002;
    double kTMax = 0.2;
    double etaMin = -5.0;
    double etaMax = 5.0;
    double phiMin = 0.0;
    double phiMax = 2.0 * std::acos(-1.0);
    double shapeSafetyFactor = 3.0;
    double boostedShapeSafetyFactor = 5.0;
    FailurePolicy failurePolicy = FailurePolicy::Skip;
    bool atomicOutput = false;
};

void printUsage(std::ostream& out) {
    out << "Usage: lowmc_hepmc3_afterburner_diagnostic [options] <input.hepmc3> <output.hepmc3>\n"
        << "\n"
        << "Reads HepMC3 ASCII events, runs the current LowMC count-plus-shape path,\n"
        << "and writes a new HepMC3 ASCII file with Low photons appended\n"
        << "for events that pass all gates. Failed events are written unchanged by\n"
        << "default and reported in the text diagnostics.\n"
        << "\n"
        << "Options:\n"
        << "  --max-events N                  Maximum events to process (default: all)\n"
        << "  --enhancement-factor F          Photon-count enhancement factor (default: 1)\n"
        << "  --multiplicity-seed N           Poisson RNG seed (default: 12345)\n"
        << "  --shape-seed N                  Shape RNG seed base (default: 67890)\n"
        << "  --min-shape-trials N            Minimum shape-unweighting trials (default: 10000)\n"
        << "  --shape-trials-per-photon N     Extra trials per requested photon (default: 2000)\n"
        << "  --shape-base-probe-bins N       Base probe bins before boosted cells (default: 11)\n"
        << "  --audit-full-base-probes        Compare reduced and full base probe maxima\n"
        << "  --targeted-base-probes          Keep full base probes near charged particles\n"
        << "  --targeted-base-probe-neighbor-scale X  Neighbor scale for targeted base probes (default: 0)\n"
        << "  --evidence-targeted-base-probes Keep full base probes in calibration-suspect cells\n"
        << "  --evidence-targeted-base-probe-neighbor-scale X  Neighbor scale for evidence base probes (default: 0)\n"
        << "  --disable-boosted-cell-probes   Disable diagnostic boosted-cell reprobing\n"
        << "  --kt-min GeV                    Minimum photon k_T (default: 0.002)\n"
        << "  --kt-max GeV                    Maximum photon k_T (default: 0.2)\n"
        << "  --eta-min X                     Minimum photon eta (default: -5)\n"
        << "  --eta-max X                     Maximum photon eta (default: 5)\n"
        << "  --phi-min X                     Minimum photon phi (default: 0)\n"
        << "  --phi-max X                     Maximum photon phi (default: 2 pi)\n"
        << "  --shape-safety-factor X         Diagnostic envelope safety factor (default: 3)\n"
        << "  --boosted-shape-safety-factor X  Diagnostic boosted-cell safety factor (default: 5)\n"
        << "  --failure-policy skip|fail      Failed-event policy (default: skip)\n"
        << "  --atomic-output                 Write to <output>.tmp and publish only on success\n";
}

std::string readStringOption(int& index, int argc, char** argv, const std::string& optionName) {
    if (index + 1 >= argc) {
        throw std::invalid_argument(optionName + " requires a value");
    }
    return argv[++index];
}

int parsePositiveInt(const std::string& text, const std::string& optionName) {
    std::size_t parsed = 0;
    const long value = std::stol(text, &parsed);
    if (parsed != text.size() || value <= 0 || value > std::numeric_limits<int>::max()) {
        throw std::invalid_argument(optionName + " must be a positive integer");
    }
    return static_cast<int>(value);
}

std::uint64_t parseUInt64(const std::string& text, const std::string& optionName) {
    std::size_t parsed = 0;
    const unsigned long long value = std::stoull(text, &parsed);
    if (parsed != text.size()) {
        throw std::invalid_argument(optionName + " must be an unsigned integer");
    }
    return static_cast<std::uint64_t>(value);
}

double parseDouble(const std::string& text, const std::string& optionName) {
    std::size_t parsed = 0;
    const double value = std::stod(text, &parsed);
    if (parsed != text.size() || !std::isfinite(value)) {
        throw std::invalid_argument(optionName + " must be finite");
    }
    return value;
}

double parsePositiveDouble(const std::string& text, const std::string& optionName) {
    const double value = parseDouble(text, optionName);
    if (value <= 0.0) {
        throw std::invalid_argument(optionName + " must be positive");
    }
    return value;
}

FailurePolicy parseFailurePolicy(const std::string& text) {
    if (text == "skip") {
        return FailurePolicy::Skip;
    }
    if (text == "fail") {
        return FailurePolicy::Fail;
    }
    throw std::invalid_argument("--failure-policy must be skip or fail");
}

Options parseOptions(int argc, char** argv) {
    Options options;
    std::vector<std::string> paths;

    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument == "--help" || argument == "-h") {
            printUsage(std::cout);
            std::exit(EXIT_SUCCESS);
        }
        if (argument == "--max-events") {
            options.maxEvents =
                parsePositiveInt(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--enhancement-factor") {
            options.enhancementFactor =
                parsePositiveDouble(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--multiplicity-seed") {
            options.multiplicitySeed =
                parseUInt64(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--shape-seed") {
            options.shapeSeed =
                parseUInt64(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--min-shape-trials") {
            options.minShapeTrials =
                parsePositiveInt(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--shape-trials-per-photon") {
            options.shapeTrialsPerPhoton =
                parsePositiveInt(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--shape-base-probe-bins") {
            options.shapeBaseProbeBins =
                parsePositiveInt(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--audit-full-base-probes") {
            options.auditFullBaseProbes = true;
            continue;
        }
        if (argument == "--targeted-base-probes") {
            options.targetedBaseProbes = true;
            continue;
        }
        if (argument == "--evidence-targeted-base-probes") {
            options.evidenceTargetedBaseProbes = true;
            continue;
        }
        if (argument == "--evidence-targeted-base-probe-neighbor-scale") {
            options.evidenceTargetedBaseProbeNeighborScale =
                parseDouble(readStringOption(index, argc, argv, argument), argument);
            if (options.evidenceTargetedBaseProbeNeighborScale < 0.0) {
                throw std::invalid_argument(
                    "--evidence-targeted-base-probe-neighbor-scale must be non-negative");
            }
            continue;
        }
        if (argument == "--disable-boosted-cell-probes") {
            options.disableBoostedCellProbes = true;
            continue;
        }
        if (argument == "--targeted-base-probe-neighbor-scale") {
            options.targetedBaseProbeNeighborScale =
                parseDouble(readStringOption(index, argc, argv, argument), argument);
            if (options.targetedBaseProbeNeighborScale < 0.0) {
                throw std::invalid_argument(
                    "--targeted-base-probe-neighbor-scale must be non-negative");
            }
            continue;
        }
        if (argument == "--kt-min") {
            options.kTMin = parseDouble(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--kt-max") {
            options.kTMax = parseDouble(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--eta-min") {
            options.etaMin = parseDouble(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--eta-max") {
            options.etaMax = parseDouble(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--phi-min") {
            options.phiMin = parseDouble(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--phi-max") {
            options.phiMax = parseDouble(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--shape-safety-factor") {
            options.shapeSafetyFactor =
                parsePositiveDouble(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--boosted-shape-safety-factor") {
            options.boostedShapeSafetyFactor =
                parsePositiveDouble(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--failure-policy") {
            options.failurePolicy = parseFailurePolicy(readStringOption(index, argc, argv, argument));
            continue;
        }
        if (argument == "--atomic-output") {
            options.atomicOutput = true;
            continue;
        }

        if (!argument.empty() && argument.front() == '-') {
            throw std::invalid_argument("unknown option: " + argument);
        }
        paths.push_back(argument);
    }

    if (paths.size() != 2) {
        throw std::invalid_argument("expected exactly input and output HepMC3 paths");
    }
    options.inputPath = paths[0];
    options.outputPath = paths[1];

    const lowmc::PhotonPhaseSpace phaseSpace{
        options.kTMin,
        options.kTMax,
        options.etaMin,
        options.etaMax,
        options.phiMin,
        options.phiMax,
    };
    if (!lowmc::validatePhotonPhaseSpace(phaseSpace).valid()) {
        throw std::invalid_argument("invalid photon phase space");
    }

    return options;
}

lowmc::AdaptiveEventYieldIntegrationConfig makeYieldConfig() {
    lowmc::AdaptiveEventYieldIntegrationConfig config;
    config.logKTBins = 16;
    config.etaSubBins = 2;
    config.phiSubBins = 2;
    config.gridConfig.etaBins = 32;
    config.gridConfig.phiBins = 32;
    config.gridConfig.refineFraction = 0.1;
    config.gridConfig.uniformFraction = 0.1;
    config.gridConfig.refinementLevels = 7;
    return config;
}

lowmc::ShapeUnweightingConfig makeShapeConfig(
    const lowmc::PhotonPhaseSpace& phaseSpace,
    const Options& options,
    int trials,
    int eventIndex) {
    lowmc::ShapeUnweightingConfig config;
    config.phaseSpace = phaseSpace;
    config.trials = trials;
    config.seed = options.shapeSeed + static_cast<std::uint64_t>(eventIndex);
    config.gridConfig = lowmc::tunedMassGridAdaptiveEtaPhiConfig();
    config.calibrationTrials = std::max(1000, std::min(20000, trials));
    config.probeBins = 11;
    config.baseProbeBins = options.shapeBaseProbeBins;
    config.auditFullBaseProbes = options.auditFullBaseProbes;
    config.targetedBaseProbes = options.targetedBaseProbes;
    config.targetedBaseProbeNeighborScale = options.targetedBaseProbeNeighborScale;
    config.evidenceTargetedBaseProbes = options.evidenceTargetedBaseProbes;
    config.evidenceTargetedBaseProbeNeighborScale =
        options.evidenceTargetedBaseProbeNeighborScale;
    config.probeRule = lowmc::ShapeUnweightingProbeRule::Interior;
    config.probePolicy =
        options.disableBoostedCellProbes
            ? lowmc::ShapeUnweightingProbePolicy::Fixed
            : lowmc::ShapeUnweightingProbePolicy::BoostedCells;
    config.boostedProbeBins = 21;
    config.boostedCellSafetyFactor = options.boostedShapeSafetyFactor;
    config.boostedCellNeighborScale = 0.0;
    config.targetedRefinementLevels = 5;
    config.targetedRefinementNeighborScale = 0.5;
    config.safetyFactor = options.shapeSafetyFactor;
    config.failurePolicy = lowmc::EnvelopeFailurePolicy::Diagnostic;
    return config;
}

const char* shapeStatusName(lowmc::ShapeUnweightingStatus status) {
    switch (status) {
        case lowmc::ShapeUnweightingStatus::Ok:
            return "ok";
        case lowmc::ShapeUnweightingStatus::EnvelopeFailed:
            return "envelope_failed";
        case lowmc::ShapeUnweightingStatus::InvalidSelection:
            return "invalid_selection";
        case lowmc::ShapeUnweightingStatus::NoTrials:
            return "no_trials";
    }
    return "unknown";
}

void printAdapterDiagnosticPrefix(
    int eventIndex,
    HepMC3::Units::MomentumUnit momentumUnit,
    const lowmc::AdapterSelectionDiagnostic& diagnostic) {
    std::cout << "event_index=" << eventIndex
              << " momentum_unit=" << HepMC3::Units::name(momentumUnit)
              << " selected_incoming_charged=" << diagnostic.selectedIncomingCharged
              << " selected_outgoing_charged=" << diagnostic.selectedOutgoingCharged
              << " ignored_neutral_candidates=" << diagnostic.ignoredNeutralCandidates
              << " rejected_charged_candidates=" << diagnostic.rejectedChargedCount()
              << " signed_charge_sum=" << diagnostic.signedChargeSum
              << " ordinary_charge_sum=" << diagnostic.ordinaryChargeSum;
}

bool failOrSkip(const Options& options, const std::string& reason) {
    if (options.failurePolicy == FailurePolicy::Fail) {
        throw std::runtime_error(reason);
    }
    return false;
}

std::string temporaryOutputPath(const std::string& outputPath) {
    return outputPath + ".tmp";
}

void removeFileIfPresent(const std::string& path) {
    if (std::remove(path.c_str()) != 0 && errno != ENOENT) {
        throw std::runtime_error("could not remove file: " + path);
    }
}

void publishAtomicOutput(const std::string& temporaryPath, const std::string& outputPath) {
    removeFileIfPresent(outputPath);
    if (std::rename(temporaryPath.c_str(), outputPath.c_str()) != 0) {
        throw std::runtime_error(
            "could not publish atomic HepMC3 output: " + temporaryPath + " -> " + outputPath);
    }
}

double secondsSince(std::chrono::steady_clock::time_point start) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}

}  // namespace

int main(int argc, char** argv) {
    std::string atomicTemporaryPath;
    try {
        const Options options = parseOptions(argc, argv);
        HepMC3::ReaderAscii reader(options.inputPath);
        if (reader.failed()) {
            throw std::runtime_error("could not open HepMC3 input file: " + options.inputPath);
        }

        const std::string writerPath =
            options.atomicOutput ? temporaryOutputPath(options.outputPath) : options.outputPath;
        if (options.atomicOutput) {
            atomicTemporaryPath = writerPath;
            removeFileIfPresent(writerPath);
        }

        HepMC3::WriterAscii writer(writerPath);
        if (writer.failed()) {
            throw std::runtime_error("could not open HepMC3 output file: " + writerPath);
        }

        const lowmc::PhotonPhaseSpace phaseSpace{
            options.kTMin,
            options.kTMax,
            options.etaMin,
            options.etaMax,
            options.phiMin,
            options.phiMax,
        };

        int eventsRead = 0;
        int eventsWritten = 0;
        int totalPhotonsAdded = 0;
        while (eventsRead < options.maxEvents && !reader.failed()) {
            HepMC3::GenEvent event;
            if (!reader.read_event(event)) {
                break;
            }
            if (event.particles().empty() && event.vertices().empty() && event.beams().empty()) {
                break;
            }

            const auto eventStart = std::chrono::steady_clock::now();
            const auto selectStart = std::chrono::steady_clock::now();
            const lowmc::AdapterParticleSelectionResult selection =
                lowmc::hepmc3::selectChargedParticles(event);
            const double selectSeconds = secondsSince(selectStart);
            printAdapterDiagnosticPrefix(eventsRead, event.momentum_unit(), selection.diagnostic);

            int photonsAdded = 0;
            std::string status = "generated";
            double lambdaEvent = 0.0;
            int requestedPhotons = 0;
            int retainedPhotons = 0;
            double maxTrialBoundRatio = 0.0;
            const char* shapeStatus = "not_run";
            double yieldSeconds = 0.0;
            double poissonSeconds = 0.0;
            double shapeSeconds = 0.0;
            double shapeGridSeconds = 0.0;
            double shapeBaseProbeSeconds = 0.0;
            double shapeCalibrationSeconds = 0.0;
            double shapeProbePlanSeconds = 0.0;
            double shapeBoostedProbeSeconds = 0.0;
            double shapeEnvelopeBuildSeconds = 0.0;
            double shapeTrialSeconds = 0.0;
            double shapeFinalizationSeconds = 0.0;
            std::size_t shapeGridCellCount = 0;
            std::size_t shapeProbeEvaluations = 0;
            int shapeCalibrationCount = 0;
            int shapeTrialCount = 0;
            int shapeAcceptedCount = 0;
            int shapeFailedEnvelopeCount = 0;
            int shapeSuspectCellCount = 0;
            int shapePromotedCellCount = 0;
            int shapeBaseProbeBins = 0;
            int shapeAuditCellCount = 0;
            int shapeAuditUnderBoundedCellCount = 0;
            double shapeAuditMaxFullOverBaseProbeRatio = 0.0;
            int shapeAuditRatioAbove110CellCount = 0;
            int shapeAuditRatioAbove150CellCount = 0;
            int shapeAuditRatioAbove200CellCount = 0;
            double shapeAuditSeconds = 0.0;
            double shapeAuditSevereRatioThreshold = 0.0;
            int shapeAuditSevereCellCount = 0;
            int shapeAuditSevereCalibrationSuspectCellCount = 0;
            int shapeAuditSevereChargedNearCellCount = 0;
            double shapeAuditSevereProposalMassSum = 0.0;
            double shapeAuditSevereMeanCellArea = 0.0;
            double shapeAuditSevereMaxCellArea = 0.0;
            double shapeAuditSevereMeanEtaWidth = 0.0;
            double shapeAuditSevereMeanPhiWidth = 0.0;
            double shapeAuditSevereMaxCalibrationOverBaseProbeRatio = 0.0;
            int shapeFullBaseProbeCellCount = 0;
            double shapeTargetedBaseProbeNeighborScale = 0.0;
            int shapeEvidenceBaseProbeCellCount = 0;
            double shapeEvidenceTargetedBaseProbeNeighborScale = 0.0;
            double writeSeconds = 0.0;

            if (selection.selectedParticles.empty()) {
                status = "skipped_invalid_selection";
                (void)failOrSkip(options, "invalid charged-particle selection");
            } else if (std::abs(selection.diagnostic.signedChargeSum) > 1.0e-12) {
                status = "skipped_nonconserved_charge";
                (void)failOrSkip(options, "selected Low current is not charge conserving");
            } else {
                const auto yieldStart = std::chrono::steady_clock::now();
                const lowmc::EventYieldEstimate yield =
                    lowmc::estimateFactorizedAdaptiveEventYield(
                        selection.selectedParticles,
                        phaseSpace,
                        makeYieldConfig());
                yieldSeconds = secondsSince(yieldStart);
                lambdaEvent = yield.lambdaEvent;

                const bool yieldGateOk =
                    yield.method == std::string("factorized-adaptive-eta-phi-midpoint") &&
                    yield.gridCells > 0 &&
                    std::isfinite(yield.lambdaEvent) &&
                    yield.lambdaEvent >= 0.0;

                if (!yieldGateOk) {
                    status = "skipped_failed_yield_gate";
                    (void)failOrSkip(options, "adaptive yield gate failed");
                } else {
                    const auto poissonStart = std::chrono::steady_clock::now();
                    const lowmc::PhotonMultiplicityDiagnostic multiplicity =
                        lowmc::samplePhotonMultiplicity(
                            yield.lambdaEvent,
                            lowmc::PhotonMultiplicityConfig{
                                options.enhancementFactor,
                                options.multiplicitySeed + static_cast<std::uint64_t>(eventsRead),
                            });
                    poissonSeconds = secondsSince(poissonStart);
                    requestedPhotons = multiplicity.sampledPhotons;

                    if (requestedPhotons == 0) {
                        status = "ok_zero_photons";
                    } else {
                        const int shapeTrials = std::max(
                            options.minShapeTrials,
                            requestedPhotons * options.shapeTrialsPerPhoton);
                        const auto shapeStart = std::chrono::steady_clock::now();
                        const lowmc::ShapeUnweightingResult shape =
                            lowmc::generateShapeUnweightedPhotonCandidates(
                                selection.selectedParticles,
                                makeShapeConfig(phaseSpace, options, shapeTrials, eventsRead));
                        shapeSeconds = secondsSince(shapeStart);
                        shapeStatus = shapeStatusName(shape.status);
                        maxTrialBoundRatio = shape.envelope.maxTrialBoundRatio;
                        shapeGridSeconds = shape.envelope.gridSeconds;
                        shapeBaseProbeSeconds = shape.envelope.baseProbeSeconds;
                        shapeCalibrationSeconds = shape.envelope.calibrationSeconds;
                        shapeProbePlanSeconds = shape.envelope.probePlanSeconds;
                        shapeBoostedProbeSeconds = shape.envelope.boostedProbeSeconds;
                        shapeEnvelopeBuildSeconds = shape.envelope.envelopeBuildSeconds;
                        shapeTrialSeconds = shape.envelope.trialSeconds;
                        shapeFinalizationSeconds = shape.envelope.finalizationSeconds;
                        shapeGridCellCount = shape.envelope.gridCellCount;
                        shapeProbeEvaluations = shape.envelope.probeEvaluations;
                        shapeCalibrationCount = shape.envelope.calibrationCount;
                        shapeTrialCount = shape.envelope.trialCount;
                        shapeAcceptedCount = shape.envelope.acceptedCount;
                        shapeFailedEnvelopeCount = shape.envelope.failedEnvelopeCount;
                        shapeSuspectCellCount = shape.envelope.suspectCellCount;
                        shapePromotedCellCount = shape.envelope.promotedCellCount;
                        shapeBaseProbeBins = shape.envelope.baseProbeBins;
                        shapeAuditCellCount = shape.envelope.auditCellCount;
                        shapeAuditUnderBoundedCellCount =
                            shape.envelope.auditUnderBoundedCellCount;
                        shapeAuditMaxFullOverBaseProbeRatio =
                            shape.envelope.auditMaxFullOverBaseProbeRatio;
                        shapeAuditRatioAbove110CellCount =
                            shape.envelope.auditRatioAbove110CellCount;
                        shapeAuditRatioAbove150CellCount =
                            shape.envelope.auditRatioAbove150CellCount;
                        shapeAuditRatioAbove200CellCount =
                            shape.envelope.auditRatioAbove200CellCount;
                        shapeAuditSeconds = shape.envelope.auditSeconds;
                        shapeAuditSevereRatioThreshold =
                            shape.envelope.auditSevereRatioThreshold;
                        shapeAuditSevereCellCount = shape.envelope.auditSevereCellCount;
                        shapeAuditSevereCalibrationSuspectCellCount =
                            shape.envelope.auditSevereCalibrationSuspectCellCount;
                        shapeAuditSevereChargedNearCellCount =
                            shape.envelope.auditSevereChargedNearCellCount;
                        shapeAuditSevereProposalMassSum =
                            shape.envelope.auditSevereProposalMassSum;
                        shapeAuditSevereMeanCellArea =
                            shape.envelope.auditSevereMeanCellArea;
                        shapeAuditSevereMaxCellArea =
                            shape.envelope.auditSevereMaxCellArea;
                        shapeAuditSevereMeanEtaWidth =
                            shape.envelope.auditSevereMeanEtaWidth;
                        shapeAuditSevereMeanPhiWidth =
                            shape.envelope.auditSevereMeanPhiWidth;
                        shapeAuditSevereMaxCalibrationOverBaseProbeRatio =
                            shape.envelope.auditSevereMaxCalibrationOverBaseProbeRatio;
                        shapeFullBaseProbeCellCount = shape.envelope.fullBaseProbeCellCount;
                        shapeTargetedBaseProbeNeighborScale =
                            shape.envelope.targetedBaseProbeNeighborScale;
                        shapeEvidenceBaseProbeCellCount =
                            shape.envelope.evidenceBaseProbeCellCount;
                        shapeEvidenceTargetedBaseProbeNeighborScale =
                            shape.envelope.evidenceTargetedBaseProbeNeighborScale;
                        retainedPhotons =
                            std::min(requestedPhotons, static_cast<int>(shape.candidates.size()));
                        const bool envelopeOk =
                            shape.status == lowmc::ShapeUnweightingStatus::Ok &&
                            shape.envelope.failedEnvelopeCount == 0;
                        const bool enoughCandidates = retainedPhotons == requestedPhotons;

                        if (!envelopeOk) {
                            status = "skipped_failed_shape_envelope";
                            (void)failOrSkip(options, "shape envelope failed");
                        } else if (!enoughCandidates) {
                            status = "skipped_insufficient_shape_candidates";
                            (void)failOrSkip(options, "insufficient accepted shape candidates");
                        } else {
                            std::vector<lowmc::ShapeUnweightedPhotonCandidate> retained(
                                shape.candidates.begin(),
                                shape.candidates.begin() + retainedPhotons);
                            const lowmc::hepmc3::HepMC3LowPhotonOutputResult output =
                                lowmc::hepmc3::appendUnweightedLowPhotons(
                                    event,
                                    retained,
                                    lowmc::hepmc3::HepMC3LowPhotonOutputOptions{
                                        options.enhancementFactor,
                                        "",
                                        yield.lambdaEvent,
                                        true,
                                        "ordinary",
                                        shapeStatus,
                                        shape.envelope.maxTrialBoundRatio,
                                        true,
                                    });
                            photonsAdded = output.photonsAdded;
                            totalPhotonsAdded += photonsAdded;
                        }
                    }
                }
            }

            const auto writeStart = std::chrono::steady_clock::now();
            writer.write_event(event);
            if (writer.failed()) {
                throw std::runtime_error("failed while writing HepMC3 output file");
            }
            writeSeconds = secondsSince(writeStart);
            ++eventsWritten;
            const auto eventEnd = std::chrono::steady_clock::now();
            const double eventSeconds =
                std::chrono::duration<double>(eventEnd - eventStart).count();

            std::cout << " photon_status=" << status
                      << " photon_mode="
                      << (options.enhancementFactor == 1.0 ? "unweighted" : "enhanced")
                      << " enhancement_factor=" << options.enhancementFactor
                      << " photon_weight=" << (1.0 / options.enhancementFactor)
                      << " lambda_event=" << lambdaEvent
                      << " requested_photons=" << requestedPhotons
                      << " retained_photons=" << retainedPhotons
                      << " photons_added=" << photonsAdded
                      << " shape_status=" << shapeStatus
                      << " max_trial_bound_ratio=" << maxTrialBoundRatio
                      << " event_seconds=" << eventSeconds
                      << " select_seconds=" << selectSeconds
                      << " yield_seconds=" << yieldSeconds
                      << " poisson_seconds=" << poissonSeconds
                      << " shape_seconds=" << shapeSeconds
                      << " shape_grid_seconds=" << shapeGridSeconds
                      << " shape_base_probe_seconds=" << shapeBaseProbeSeconds
                      << " shape_calibration_seconds=" << shapeCalibrationSeconds
                      << " shape_probe_plan_seconds=" << shapeProbePlanSeconds
                      << " shape_boosted_probe_seconds=" << shapeBoostedProbeSeconds
                      << " shape_envelope_build_seconds=" << shapeEnvelopeBuildSeconds
                      << " shape_trial_seconds=" << shapeTrialSeconds
                      << " shape_finalization_seconds=" << shapeFinalizationSeconds
                      << " shape_grid_cell_count=" << shapeGridCellCount
                      << " shape_probe_evaluations=" << shapeProbeEvaluations
                      << " shape_calibration_count=" << shapeCalibrationCount
                      << " shape_trial_count=" << shapeTrialCount
                      << " shape_accepted_count=" << shapeAcceptedCount
                      << " shape_failed_envelope_count=" << shapeFailedEnvelopeCount
                      << " shape_suspect_cell_count=" << shapeSuspectCellCount
                      << " shape_promoted_cell_count=" << shapePromotedCellCount
                      << " shape_base_probe_bins=" << shapeBaseProbeBins
                      << " shape_audit_cell_count=" << shapeAuditCellCount
                      << " shape_audit_under_bounded_cell_count="
                      << shapeAuditUnderBoundedCellCount
                      << " shape_audit_max_full_over_base_probe_ratio="
                      << shapeAuditMaxFullOverBaseProbeRatio
                      << " shape_audit_ratio_above_110_cell_count="
                      << shapeAuditRatioAbove110CellCount
                      << " shape_audit_ratio_above_150_cell_count="
                      << shapeAuditRatioAbove150CellCount
                      << " shape_audit_ratio_above_200_cell_count="
                      << shapeAuditRatioAbove200CellCount
                      << " shape_audit_seconds=" << shapeAuditSeconds
                      << " shape_audit_severe_ratio_threshold="
                      << shapeAuditSevereRatioThreshold
                      << " shape_audit_severe_cell_count=" << shapeAuditSevereCellCount
                      << " shape_audit_severe_calibration_suspect_cell_count="
                      << shapeAuditSevereCalibrationSuspectCellCount
                      << " shape_audit_severe_charged_near_cell_count="
                      << shapeAuditSevereChargedNearCellCount
                      << " shape_audit_severe_proposal_mass_sum="
                      << shapeAuditSevereProposalMassSum
                      << " shape_audit_severe_mean_cell_area="
                      << shapeAuditSevereMeanCellArea
                      << " shape_audit_severe_max_cell_area="
                      << shapeAuditSevereMaxCellArea
                      << " shape_audit_severe_mean_eta_width="
                      << shapeAuditSevereMeanEtaWidth
                      << " shape_audit_severe_mean_phi_width="
                      << shapeAuditSevereMeanPhiWidth
                      << " shape_audit_severe_max_calibration_over_base_probe_ratio="
                      << shapeAuditSevereMaxCalibrationOverBaseProbeRatio
                      << " shape_full_base_probe_cell_count=" << shapeFullBaseProbeCellCount
                      << " shape_targeted_base_probe_neighbor_scale="
                      << shapeTargetedBaseProbeNeighborScale
                      << " shape_evidence_base_probe_cell_count="
                      << shapeEvidenceBaseProbeCellCount
                      << " shape_evidence_targeted_base_probe_neighbor_scale="
                      << shapeEvidenceTargetedBaseProbeNeighborScale
                      << " write_seconds=" << writeSeconds
                      << '\n';

            ++eventsRead;
        }

        writer.close();
        if (eventsRead == 0) {
            throw std::runtime_error("no events read from HepMC3 input file: " + options.inputPath);
        }
        if (options.atomicOutput) {
            publishAtomicOutput(writerPath, options.outputPath);
            atomicTemporaryPath.clear();
        }

        std::cerr << "events_read=" << eventsRead
                  << " events_written=" << eventsWritten
                  << " photons_added=" << totalPhotonsAdded
                  << " output=" << options.outputPath << '\n';
    } catch (const std::exception& error) {
        if (!atomicTemporaryPath.empty()) {
            try {
                removeFileIfPresent(atomicTemporaryPath);
            } catch (const std::exception& cleanupError) {
                std::cerr << "FAILED: could not remove atomic temporary output: "
                          << cleanupError.what() << '\n';
            }
        }
        std::cerr << "FAILED: " << error.what() << '\n';
        printUsage(std::cerr);
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
