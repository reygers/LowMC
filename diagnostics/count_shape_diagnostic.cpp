#include "ParticleInput.h"
#include "lowmc/LowMC.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

enum class ValidationClass {
    Ordinary,
    Stress,
};

struct Options {
    std::string outputPath = "count_shape_diagnostic.csv";
    std::string candidatesCsv;
    std::string particlesPath;
    std::string topology = "two-pion";
    double momentum = 10.0;
    double kTMin = 0.002;
    double kTMax = 0.2;
    double etaMin = -5.0;
    double etaMax = 5.0;
    double phiMin = 0.0;
    double phiMax = 2.0 * std::acos(-1.0);
    double oversamplingFactor = 1.0;
    std::uint64_t multiplicitySeed = 12345;
    std::uint64_t shapeSeed = 67890;
    int minShapeTrials = 10000;
    int shapeTrialsPerRequestedPhoton = 2000;
    ValidationClass validationClass = ValidationClass::Ordinary;
    bool validationClassExplicit = false;
};

void printUsage(const char* programName) {
    std::cout
        << "Usage: " << programName << " [summary.csv] [options]\n"
        << "Draw a diagnostic Poisson count and shape-unweighted photon candidates.\n\n"
        << "Options:\n"
        << "  --particles PATH                 Read charged particles from a text fixture\n"
        << "  --topology two-pion|event4        Toy topology if no fixture is given\n"
        << "  --momentum GeV                   Two-pion toy momentum\n"
        << "  --kt-min GeV                     Photon kT lower bound\n"
        << "  --kt-max GeV                     Photon kT upper bound\n"
        << "  --eta-min X                      Photon eta lower bound\n"
        << "  --eta-max X                      Photon eta upper bound\n"
        << "  --phi-min X                      Photon phi lower bound\n"
        << "  --phi-max X                      Photon phi upper bound\n"
        << "  --oversampling-factor F          Diagnostic oversampling factor\n"
        << "  --multiplicity-seed N            Poisson RNG seed\n"
        << "  --shape-seed N                   Shape-unweighting RNG seed\n"
        << "  --min-shape-trials N             Minimum shape-unweighting trials\n"
        << "  --shape-trials-per-photon N      Extra trials per requested photon\n"
        << "  --validation-class ordinary|stress  Yield gate class\n"
        << "  --candidates-csv PATH            Optional retained-candidate CSV\n";
}

std::string readStringOption(int& index, int argc, char** argv, const std::string& optionName) {
    if (index + 1 >= argc) {
        throw std::invalid_argument(optionName + " requires a value");
    }
    return argv[++index];
}

int parseNonNegativeInt(const std::string& text, const std::string& optionName) {
    std::size_t parsed = 0;
    const int value = std::stoi(text, &parsed);
    if (parsed != text.size() || value < 0) {
        throw std::invalid_argument(optionName + " must be a non-negative integer");
    }
    return value;
}

double parseDouble(const std::string& text, const std::string& optionName) {
    std::size_t parsed = 0;
    const double value = std::stod(text, &parsed);
    if (parsed != text.size() || !std::isfinite(value)) {
        throw std::invalid_argument(optionName + " must be a finite number");
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

std::uint64_t parseUInt64(const std::string& text, const std::string& optionName) {
    std::size_t parsed = 0;
    const unsigned long long value = std::stoull(text, &parsed);
    if (parsed != text.size()) {
        throw std::invalid_argument(optionName + " must be an unsigned integer");
    }
    return static_cast<std::uint64_t>(value);
}

const char* validationClassName(ValidationClass value) {
    switch (value) {
        case ValidationClass::Ordinary:
            return "ordinary";
        case ValidationClass::Stress:
            return "stress";
    }
    return "unknown";
}

ValidationClass parseValidationClass(const std::string& text) {
    if (text == "ordinary") {
        return ValidationClass::Ordinary;
    }
    if (text == "stress") {
        return ValidationClass::Stress;
    }
    throw std::invalid_argument("--validation-class must be ordinary or stress");
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

Options parseOptions(int argc, char** argv) {
    Options options;
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument == "--help" || argument == "-h") {
            printUsage(argv[0]);
            std::exit(EXIT_SUCCESS);
        }
        if (argument == "--particles") {
            options.particlesPath = readStringOption(index, argc, argv, argument);
            continue;
        }
        if (argument == "--topology") {
            options.topology = readStringOption(index, argc, argv, argument);
            if (options.topology != "two-pion" && options.topology != "event4") {
                throw std::invalid_argument("--topology must be two-pion or event4");
            }
            continue;
        }
        if (argument == "--momentum") {
            options.momentum =
                parsePositiveDouble(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--kt-min") {
            options.kTMin =
                parsePositiveDouble(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--kt-max") {
            options.kTMax =
                parsePositiveDouble(readStringOption(index, argc, argv, argument), argument);
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
        if (argument == "--oversampling-factor") {
            options.oversamplingFactor =
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
                parseNonNegativeInt(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--shape-trials-per-photon") {
            options.shapeTrialsPerRequestedPhoton =
                parseNonNegativeInt(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--validation-class") {
            options.validationClass =
                parseValidationClass(readStringOption(index, argc, argv, argument));
            options.validationClassExplicit = true;
            continue;
        }
        if (argument == "--candidates-csv") {
            options.candidatesCsv = readStringOption(index, argc, argv, argument);
            continue;
        }
        if (!argument.empty() && argument[0] == '-') {
            throw std::invalid_argument("unknown option: " + argument);
        }
        if (options.outputPath != "count_shape_diagnostic.csv") {
            throw std::invalid_argument("only one summary CSV path can be specified");
        }
        options.outputPath = argument;
    }
    return options;
}

lowmc::ChargedParticle makePion(double momentum, double eta, double phi, double charge) {
    const double mass = 0.13957039;
    const double pT = momentum / std::cosh(eta);
    const double px = pT * std::cos(phi);
    const double py = pT * std::sin(phi);
    const double pz = pT * std::sinh(eta);
    const double energy = std::sqrt(momentum * momentum + mass * mass);
    return {lowmc::FourVector{energy, px, py, pz}, charge, lowmc::InOutSign::Outgoing};
}

std::vector<lowmc::ChargedParticle> makeParticles(const Options& options) {
    if (!options.particlesPath.empty()) {
        return lowmc_diagnostics::readChargedParticlesFromTextFile(options.particlesPath);
    }
    if (options.topology == "event4") {
        return {
            makePion(80.0, 0.20, 0.70, +1.0),
            makePion(55.0, -0.80, 2.40, -1.0),
            makePion(25.0, 1.60, 4.30, +1.0),
            makePion(12.0, -2.20, 5.40, -1.0),
        };
    }
    const double phi = 0.5 * std::acos(-1.0);
    return {
        makePion(options.momentum, 0.0, phi, +1.0),
        makePion(options.momentum, 0.0, phi + std::acos(-1.0), -1.0),
    };
}

lowmc::AdaptiveEventYieldIntegrationConfig makeYieldConfig(ValidationClass validationClass) {
    lowmc::AdaptiveEventYieldIntegrationConfig config;
    config.logKTBins = 16;
    config.etaSubBins = 2;
    config.phiSubBins = 2;
    config.gridConfig.etaBins = 32;
    config.gridConfig.phiBins = 32;
    config.gridConfig.refineFraction = 0.1;
    config.gridConfig.uniformFraction = 0.1;
    config.gridConfig.refinementLevels =
        validationClass == ValidationClass::Stress ? 9 : 7;
    return config;
}

lowmc::ShapeUnweightingConfig makeShapeConfig(
    const lowmc::PhotonPhaseSpace& phaseSpace,
    const Options& options,
    int trials) {
    lowmc::ShapeUnweightingConfig config;
    config.phaseSpace = phaseSpace;
    config.trials = trials;
    config.seed = options.shapeSeed;
    config.gridConfig = lowmc::tunedMassGridAdaptiveEtaPhiConfig();
    config.calibrationTrials = std::max(1000, std::min(20000, trials));
    config.probeBins = 11;
    config.probeRule = lowmc::ShapeUnweightingProbeRule::Interior;
    config.probePolicy = lowmc::ShapeUnweightingProbePolicy::BoostedCells;
    config.boostedProbeBins = 21;
    config.boostedCellSafetyFactor = 5.0;
    config.boostedCellNeighborScale = 0.0;
    config.targetedRefinementLevels = 5;
    config.targetedRefinementNeighborScale = 0.5;
    config.safetyFactor = 3.0;
    config.failurePolicy = lowmc::EnvelopeFailurePolicy::Diagnostic;
    return config;
}

bool passesValidationGate(
    ValidationClass validationClass,
    const lowmc::EventYieldEstimate& yield) {
    if (validationClass == ValidationClass::Ordinary) {
        return yield.method == std::string("adaptive-log-kt-grid-eta-phi-midpoint") &&
               yield.gridCells > 0;
    }
    return yield.method == std::string("adaptive-log-kt-grid-eta-phi-midpoint") &&
           yield.gridCells > 0;
}

void writeCandidatesCsv(
    const std::string& path,
    const std::vector<lowmc::ShapeUnweightedPhotonCandidate>& candidates,
    int retainedCount,
    double analysisWeight) {
    std::ofstream output(path);
    if (!output) {
        throw std::runtime_error("could not open candidates CSV: " + path);
    }
    output << "index,kT,eta,phi,analysisWeight,E,px,py,pz\n";
    output << std::setprecision(17);
    const int count =
        std::min(retainedCount, static_cast<int>(candidates.size()));
    for (int index = 0; index < count; ++index) {
        const lowmc::ShapeUnweightedPhotonCandidate& candidate =
            candidates[static_cast<std::size_t>(index)];
        output << index << ','
               << candidate.kT << ','
               << candidate.eta << ','
               << candidate.phi << ','
               << analysisWeight << ','
               << candidate.photon.momentum.E << ','
               << candidate.photon.momentum.px << ','
               << candidate.photon.momentum.py << ','
               << candidate.photon.momentum.pz << '\n';
    }
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    try {
        options = parseOptions(argc, argv);
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << "\n\n";
        printUsage(argv[0]);
        return EXIT_FAILURE;
    }

    try {
        if (!options.validationClassExplicit && options.momentum >= 1000.0 &&
            options.particlesPath.empty()) {
            options.validationClass = ValidationClass::Stress;
        }

        const std::vector<lowmc::ChargedParticle> particles = makeParticles(options);
        const lowmc::PhotonPhaseSpace phaseSpace{
            options.kTMin,
            options.kTMax,
            options.etaMin,
            options.etaMax,
            options.phiMin,
            options.phiMax,
        };
        const lowmc::AdaptiveEventYieldIntegrationConfig yieldConfig =
            makeYieldConfig(options.validationClass);
        const lowmc::EventYieldEstimate yield =
            lowmc::estimateAdaptiveEventYield(particles, phaseSpace, yieldConfig);
        const lowmc::PhotonMultiplicityDiagnostic multiplicity =
            lowmc::samplePhotonMultiplicity(
                yield.lambdaEvent,
                lowmc::PhotonMultiplicityConfig{
                    options.oversamplingFactor,
                    options.multiplicitySeed,
                });

        const int requestedPhotons = multiplicity.sampledPhotons;
        const int shapeTrials = std::max(
            options.minShapeTrials,
            requestedPhotons * options.shapeTrialsPerRequestedPhoton);
        const lowmc::ShapeUnweightingResult shape =
            lowmc::generateShapeUnweightedPhotonCandidates(
                particles,
                makeShapeConfig(phaseSpace, options, shapeTrials));
        const int retainedPhotons = std::min(
            requestedPhotons,
            static_cast<int>(shape.candidates.size()));
        const bool enoughCandidates = retainedPhotons == requestedPhotons;
        const bool envelopeOk =
            shape.status == lowmc::ShapeUnweightingStatus::Ok &&
            shape.envelope.failedEnvelopeCount == 0;
        const bool gateOk = passesValidationGate(options.validationClass, yield);
        const bool eventOk = gateOk && envelopeOk && enoughCandidates;

        if (!options.candidatesCsv.empty()) {
            writeCandidatesCsv(
                options.candidatesCsv,
                shape.candidates,
                retainedPhotons,
                multiplicity.sampleWeight);
        }

        std::ofstream output(options.outputPath);
        if (!output) {
            std::cerr << "Could not open output CSV: " << options.outputPath << '\n';
            return EXIT_FAILURE;
        }

        output << "status,validation_class,validation_gate_ok,lambda_event,"
               << "yield_method,yield_log_kt_bins,yield_grid_cells,"
               << "multiplicity_seed,shape_seed,oversampling_factor,poisson_mean,"
               << "requested_photons,retained_photons,analysis_weight,"
               << "shape_status,shape_trials,shape_accepted_candidates,enough_candidates,"
               << "envelope_ok,failed_envelope_count,max_trial_bound_ratio,"
               << "envelope_integral,shape_grid_cells,particle_count,topology,particle_file\n";
        output << std::setprecision(17)
               << (eventOk ? "ok" : "incomplete") << ','
               << validationClassName(options.validationClass) << ','
               << (gateOk ? 1 : 0) << ','
               << yield.lambdaEvent << ','
               << yield.method << ','
               << yield.logKTBins << ','
               << yield.gridCells << ','
               << multiplicity.seed << ','
               << options.shapeSeed << ','
               << multiplicity.oversamplingFactor << ','
               << multiplicity.poissonMean << ','
               << requestedPhotons << ','
               << retainedPhotons << ','
               << multiplicity.sampleWeight << ','
               << shapeStatusName(shape.status) << ','
               << shape.envelope.trialCount << ','
               << shape.candidates.size() << ','
               << (enoughCandidates ? 1 : 0) << ','
               << (envelopeOk ? 1 : 0) << ','
               << shape.envelope.failedEnvelopeCount << ','
               << shape.envelope.maxTrialBoundRatio << ','
               << shape.envelope.envelopeIntegral << ','
               << shape.envelope.gridCellCount << ','
               << particles.size() << ','
               << (options.particlesPath.empty() ? options.topology : "particles") << ','
               << (options.particlesPath.empty() ? "none" : options.particlesPath)
               << '\n';

        std::cerr << "status=" << (eventOk ? "ok" : "incomplete")
                  << " validation_class=" << validationClassName(options.validationClass)
                  << " lambda_event=" << yield.lambdaEvent
                  << " requested_photons=" << requestedPhotons
                  << " retained_photons=" << retainedPhotons
                  << " shape_status=" << shapeStatusName(shape.status)
                  << " envelope_ok=" << (envelopeOk ? 1 : 0)
                  << " validation_gate_ok=" << (gateOk ? 1 : 0)
                  << '\n';
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
