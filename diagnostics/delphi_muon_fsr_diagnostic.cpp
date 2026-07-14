#include "lowmc/LowMC.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr double MuonMassGeV = 0.1056583755;

struct Options {
    std::string outputPath = "delphi_muon_fsr_samples.csv";
    int trials = 300000;
    int calibrationTrials = 1000;
    std::uint64_t seed = 12345;
    double gamma = 430.0;
    double thetaMaxMrad = 10.0;
    double phaseHalfWidthMrad = 12.0;
    double kTMin = 0.002;
    double kTMax = 0.2;
    double safetyFactor = 3.0;
    int probeBins = 11;
    lowmc::GridAdaptiveEtaPhiProposalConfig gridConfig =
        lowmc::tunedMassGridAdaptiveEtaPhiConfig();
};

void printUsage(std::ostream& out, const char* programName) {
    out << "Usage: " << programName << " [output.csv] [options]\n"
        << "\n"
        << "Generate unweighted LowMC photon-angle samples for a fixed back-to-back\n"
        << "mu+ mu- topology inspired by DELPHI Fig. 1a. The photon phase space is a\n"
        << "small eta-phi square around the positive muon direction; this diagnostic is\n"
        << "for shape comparison, not event-record production.\n"
        << "\n"
        << "Options:\n"
        << "  --trials N                  Trial photons for accept/reject (default: 300000)\n"
        << "  --calibration-trials N      Random calibration trials (default: 1000)\n"
        << "  --seed N                    Random seed (default: 12345)\n"
        << "  --gamma X                   Muon Lorentz gamma (default: 430)\n"
        << "  --theta-max-mrad X          Histogram-interest angle scale (default: 10)\n"
        << "  --phase-half-width-mrad X   Eta/phi half-width around muon (default: 12)\n"
        << "  --kt-min X                  Minimum photon k_T in GeV (default: 0.002)\n"
        << "  --kt-max X                  Maximum photon k_T in GeV (default: 0.2)\n"
        << "  --safety-factor X           Envelope safety factor (default: 3)\n"
        << "  --probe-bins N              Probe lattice bins per cell dimension (default: 11)\n"
        << "  --grid-eta-bins N           Initial eta bins for grid proposal (default: tuned)\n"
        << "  --grid-phi-bins N           Initial phi bins for grid proposal (default: tuned)\n"
        << "  --grid-refinement-levels N  Grid refinement levels (default: tuned)\n"
        << "  --grid-refine-fraction X    Fraction of cells refined per level (default: tuned)\n"
        << "  --grid-uniform-fraction X   Uniform proposal floor fraction (default: tuned)\n";
}

std::string readStringOption(int& index, int argc, char** argv, const std::string& optionName) {
    if (index + 1 >= argc) {
        throw std::invalid_argument(optionName + " requires a value");
    }
    return argv[++index];
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

int parseNonNegativeInt(const std::string& text, const std::string& optionName) {
    std::size_t parsed = 0;
    const int value = std::stoi(text, &parsed);
    if (parsed != text.size() || value < 0) {
        throw std::invalid_argument(optionName + " must be a non-negative integer");
    }
    return value;
}

int parsePositiveInt(const std::string& text, const std::string& optionName) {
    const int value = parseNonNegativeInt(text, optionName);
    if (value <= 0) {
        throw std::invalid_argument(optionName + " must be a positive integer");
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

Options parseOptions(int argc, char** argv) {
    Options options;
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument == "--help" || argument == "-h") {
            printUsage(std::cout, argv[0]);
            std::exit(EXIT_SUCCESS);
        }
        if (argument == "--trials") {
            options.trials = parsePositiveInt(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--calibration-trials") {
            options.calibrationTrials =
                parseNonNegativeInt(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--seed") {
            options.seed = parseUInt64(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--gamma") {
            options.gamma = parsePositiveDouble(readStringOption(index, argc, argv, argument), argument);
            if (options.gamma <= 1.0) {
                throw std::invalid_argument("--gamma must be greater than one");
            }
            continue;
        }
        if (argument == "--theta-max-mrad") {
            options.thetaMaxMrad =
                parsePositiveDouble(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--phase-half-width-mrad") {
            options.phaseHalfWidthMrad =
                parsePositiveDouble(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--kt-min") {
            options.kTMin = parsePositiveDouble(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--kt-max") {
            options.kTMax = parsePositiveDouble(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--safety-factor") {
            options.safetyFactor =
                parsePositiveDouble(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--probe-bins") {
            options.probeBins = parsePositiveInt(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--grid-eta-bins") {
            options.gridConfig.etaBins =
                parsePositiveInt(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--grid-phi-bins") {
            options.gridConfig.phiBins =
                parsePositiveInt(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--grid-refinement-levels") {
            options.gridConfig.refinementLevels =
                parseNonNegativeInt(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--grid-refine-fraction") {
            options.gridConfig.refineFraction =
                parsePositiveDouble(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--grid-uniform-fraction") {
            options.gridConfig.uniformFraction =
                parseDouble(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (!argument.empty() && argument[0] == '-') {
            throw std::invalid_argument("unknown option: " + argument);
        }
        if (options.outputPath != "delphi_muon_fsr_samples.csv") {
            throw std::invalid_argument("only one output CSV path can be specified");
        }
        options.outputPath = argument;
    }
    if (options.kTMax <= options.kTMin) {
        throw std::invalid_argument("--kt-max must be larger than --kt-min");
    }
    if (options.phaseHalfWidthMrad <= options.thetaMaxMrad) {
        throw std::invalid_argument("--phase-half-width-mrad should exceed --theta-max-mrad");
    }
    return options;
}

std::vector<lowmc::ChargedParticle> makeBackToBackMuons(double gamma) {
    const double energy = gamma * MuonMassGeV;
    const double momentum = MuonMassGeV * std::sqrt(gamma * gamma - 1.0);
    return {
        {lowmc::FourVector{energy, momentum, 0.0, 0.0}, +1.0, lowmc::InOutSign::Outgoing},
        {lowmc::FourVector{energy, -momentum, 0.0, 0.0}, -1.0, lowmc::InOutSign::Outgoing},
    };
}

double vectorMagnitude(const lowmc::FourVector& vector) {
    return std::sqrt(vector.px * vector.px + vector.py * vector.py + vector.pz * vector.pz);
}

double openingAngle(const lowmc::FourVector& left, const lowmc::FourVector& right) {
    const double leftMagnitude = vectorMagnitude(left);
    const double rightMagnitude = vectorMagnitude(right);
    if (!(leftMagnitude > 0.0) || !(rightMagnitude > 0.0)) {
        throw std::runtime_error("cannot compute opening angle with zero momentum");
    }
    const double cosine = std::clamp(
        (left.px * right.px + left.py * right.py + left.pz * right.pz) /
            (leftMagnitude * rightMagnitude),
        -1.0,
        1.0);
    return std::acos(cosine);
}

void writeSamplesCsv(
    const std::string& path,
    const std::vector<lowmc::ChargedParticle>& muons,
    const lowmc::ShapeUnweightingResult& result) {
    std::ofstream output(path);
    if (!output) {
        throw std::runtime_error("could not open output CSV: " + path);
    }

    output << "theta_rad,theta_mrad,theta_plus_rad,theta_minus_rad,nearest_muon,kT,eta,phi,E,px,py,pz\n";
    for (const lowmc::ShapeUnweightedPhotonCandidate& candidate : result.candidates) {
        const double thetaPlus = openingAngle(candidate.photon.momentum, muons[0].momentum);
        const double thetaMinus = openingAngle(candidate.photon.momentum, muons[1].momentum);
        const bool nearestPlus = thetaPlus <= thetaMinus;
        const double theta = nearestPlus ? thetaPlus : thetaMinus;
        output << theta << ','
               << theta * 1000.0 << ','
               << thetaPlus << ','
               << thetaMinus << ','
               << (nearestPlus ? "plus" : "minus") << ','
               << candidate.kT << ','
               << candidate.eta << ','
               << candidate.phi << ','
               << candidate.photon.momentum.E << ','
               << candidate.photon.momentum.px << ','
               << candidate.photon.momentum.py << ','
               << candidate.photon.momentum.pz << '\n';
    }
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Options options = parseOptions(argc, argv);
        const std::vector<lowmc::ChargedParticle> muons = makeBackToBackMuons(options.gamma);
        const double halfWidth = 0.001 * options.phaseHalfWidthMrad;

        lowmc::ShapeUnweightingConfig config;
        config.phaseSpace = lowmc::PhotonPhaseSpace{
            options.kTMin,
            options.kTMax,
            -halfWidth,
            halfWidth,
            -halfWidth,
            halfWidth,
        };
        config.trials = options.trials;
        config.seed = options.seed;
        config.gridConfig = options.gridConfig;
        config.calibrationTrials = options.calibrationTrials;
        config.probeBins = options.probeBins;
        config.safetyFactor = options.safetyFactor;
        config.failurePolicy = lowmc::EnvelopeFailurePolicy::Diagnostic;

        const lowmc::ShapeUnweightingResult result =
            lowmc::generateShapeUnweightedPhotonCandidates(muons, config);
        writeSamplesCsv(options.outputPath, muons, result);

        const double momentum = vectorMagnitude(muons[0].momentum);
        std::cout << "status=" << (result.status == lowmc::ShapeUnweightingStatus::Ok ? "ok" : "not_ok")
                  << " trials=" << result.envelope.trialCount
                  << " accepted_candidates=" << result.candidates.size()
                  << " failed_envelope_count=" << result.envelope.failedEnvelopeCount
                  << " max_trial_bound_ratio=" << result.envelope.maxTrialBoundRatio
                  << " observed_acceptance_efficiency="
                  << result.envelope.observedAcceptanceEfficiency
                  << " gamma=" << options.gamma
                  << " muon_momentum_GeV=" << momentum
                  << " expected_peak_mrad=" << std::sqrt(3.0) / options.gamma * 1000.0
                  << " theta_max_mrad=" << options.thetaMaxMrad
                  << " phase_half_width_mrad=" << options.phaseHalfWidthMrad
                  << " output=" << options.outputPath
                  << '\n';
    } catch (const std::exception& error) {
        std::cerr << "FAILED: " << error.what() << "\n\n";
        printUsage(std::cerr, argv[0]);
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
