#include "lowmc/LowMC.h"

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct Options {
    std::string outputPath = "weighted_sampler_study.csv";
    int samples = 1000;
    std::uint64_t seed = 12345;
    double kT = 0.02;
    double particleEta = 0.0;
    double particlePhi = 0.5 * std::acos(-1.0);
    double etaMin = -5.0;
    double etaMax = 5.0;
    double phiMin = 0.0;
    double phiMax = 2.0 * std::acos(-1.0);
    std::vector<double> momenta = {1.0, 10.0, 100.0, 1000.0};
    std::string proposal = "uniform";
    lowmc::AdaptiveEtaPhiProposalConfig adaptiveConfig;
    lowmc::SmoothEtaPhiProposalConfig smoothConfig;
    lowmc::LegCenteredEtaPhiProposalConfig legConfig;
    lowmc::GridAdaptiveEtaPhiProposalConfig gridConfig;
};

void printUsage(const char* programName) {
    std::cout << "Usage: " << programName
              << " [output.csv] [--samples n] [--seed n] [--kt GeV]\n"
              << "       [--momentum GeV]... [--particle-eta eta] [--particle-phi phi]\n"
              << "       [--eta-min eta] [--eta-max eta] [--phi-min phi] [--phi-max phi]\n"
              << "       [--proposal uniform|adaptive|smooth|leg|grid]\n\n"
              << "       [--adaptive-uniform-fraction f] [--adaptive-width-scale s]\n"
              << "       [--adaptive-min-half-width w] [--adaptive-max-half-width w]\n\n"
              << "       [--smooth-uniform-fraction f] [--smooth-width-scale s]\n"
              << "       [--smooth-min-half-width w] [--smooth-max-half-width w]\n\n"
              << "       [--leg-uniform-fraction f] [--leg-minimum-one-minus-beta x]\n\n"
              << "       [--grid-eta-bins n] [--grid-phi-bins n]\n"
              << "       [--grid-refinement-levels n] [--grid-refine-fraction f]\n"
              << "       [--grid-uniform-fraction f]\n"
              << "       [--grid-refinement-score mass|variance]\n"
              << "       [--grid-variance-probe-bins n]\n\n"
              << "Write one weight-summary row per toy charged-pion momentum.\n"
              << "The default fixed-kT eta-phi proposal is uniform and diagnostic-only.\n";
}

double parseDouble(const std::string& text, const std::string& optionName) {
    std::size_t parsed = 0;
    const double value = std::stod(text, &parsed);
    if (parsed != text.size()) {
        throw std::invalid_argument(optionName + " must be a number");
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

int parsePositiveInt(const std::string& text, const std::string& optionName) {
    std::size_t parsed = 0;
    const int value = std::stoi(text, &parsed);
    if (parsed != text.size() || value <= 0) {
        throw std::invalid_argument(optionName + " must be a positive integer");
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

std::uint64_t parseUInt64(const std::string& text, const std::string& optionName) {
    std::size_t parsed = 0;
    const unsigned long long value = std::stoull(text, &parsed);
    if (parsed != text.size()) {
        throw std::invalid_argument(optionName + " must be an unsigned integer");
    }
    return static_cast<std::uint64_t>(value);
}

std::string readStringOption(int& index, int argc, char** argv, const std::string& optionName) {
    if (index + 1 >= argc) {
        throw std::invalid_argument(optionName + " requires a value");
    }
    return argv[++index];
}

lowmc::GridAdaptiveRefinementScore parseGridRefinementScore(
    const std::string& text,
    const std::string& optionName) {
    if (text == "mass") {
        return lowmc::GridAdaptiveRefinementScore::TargetMass;
    }
    if (text == "variance") {
        return lowmc::GridAdaptiveRefinementScore::Variance;
    }
    throw std::invalid_argument(optionName + " must be mass or variance");
}

std::string gridRefinementScoreName(lowmc::GridAdaptiveRefinementScore score) {
    if (score == lowmc::GridAdaptiveRefinementScore::Variance) {
        return "variance";
    }
    return "mass";
}

double readDoubleOption(int& index, int argc, char** argv, const std::string& optionName) {
    return parseDouble(readStringOption(index, argc, argv, optionName), optionName);
}

Options parseOptions(int argc, char** argv) {
    Options options;
    bool customMomenta = false;

    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument == "--help" || argument == "-h") {
            printUsage(argv[0]);
            std::exit(EXIT_SUCCESS);
        }

        if (argument == "--samples") {
            options.samples =
                parsePositiveInt(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--seed") {
            options.seed = parseUInt64(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--kt") {
            options.kT =
                parsePositiveDouble(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--momentum") {
            if (!customMomenta) {
                options.momenta.clear();
                customMomenta = true;
            }
            options.momenta.push_back(
                parsePositiveDouble(readStringOption(index, argc, argv, argument), argument));
            continue;
        }
        if (argument == "--particle-eta") {
            options.particleEta = readDoubleOption(index, argc, argv, argument);
            continue;
        }
        if (argument == "--particle-phi") {
            options.particlePhi = readDoubleOption(index, argc, argv, argument);
            continue;
        }
        if (argument == "--eta-min") {
            options.etaMin = readDoubleOption(index, argc, argv, argument);
            continue;
        }
        if (argument == "--eta-max") {
            options.etaMax = readDoubleOption(index, argc, argv, argument);
            continue;
        }
        if (argument == "--phi-min") {
            options.phiMin = readDoubleOption(index, argc, argv, argument);
            continue;
        }
        if (argument == "--phi-max") {
            options.phiMax = readDoubleOption(index, argc, argv, argument);
            continue;
        }
        if (argument == "--proposal") {
            options.proposal = readStringOption(index, argc, argv, argument);
            if (options.proposal != "uniform" &&
                options.proposal != "adaptive" &&
                options.proposal != "smooth" &&
                options.proposal != "leg" &&
                options.proposal != "grid") {
                throw std::invalid_argument(
                    "--proposal must be uniform, adaptive, smooth, leg, or grid");
            }
            continue;
        }
        if (argument == "--adaptive-uniform-fraction") {
            options.adaptiveConfig.uniformFraction =
                readDoubleOption(index, argc, argv, argument);
            continue;
        }
        if (argument == "--adaptive-width-scale") {
            options.adaptiveConfig.widthScale = readDoubleOption(index, argc, argv, argument);
            continue;
        }
        if (argument == "--adaptive-min-half-width") {
            options.adaptiveConfig.minHalfWidth = readDoubleOption(index, argc, argv, argument);
            continue;
        }
        if (argument == "--adaptive-max-half-width") {
            options.adaptiveConfig.maxHalfWidth = readDoubleOption(index, argc, argv, argument);
            continue;
        }
        if (argument == "--smooth-uniform-fraction") {
            options.smoothConfig.uniformFraction =
                readDoubleOption(index, argc, argv, argument);
            continue;
        }
        if (argument == "--smooth-width-scale") {
            options.smoothConfig.widthScale = readDoubleOption(index, argc, argv, argument);
            continue;
        }
        if (argument == "--smooth-min-half-width") {
            options.smoothConfig.minHalfWidth = readDoubleOption(index, argc, argv, argument);
            continue;
        }
        if (argument == "--smooth-max-half-width") {
            options.smoothConfig.maxHalfWidth = readDoubleOption(index, argc, argv, argument);
            continue;
        }
        if (argument == "--leg-uniform-fraction") {
            options.legConfig.uniformFraction =
                readDoubleOption(index, argc, argv, argument);
            continue;
        }
        if (argument == "--leg-minimum-one-minus-beta") {
            options.legConfig.minimumOneMinusBeta =
                readDoubleOption(index, argc, argv, argument);
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
            options.gridConfig.refineFraction = readDoubleOption(index, argc, argv, argument);
            continue;
        }
        if (argument == "--grid-uniform-fraction") {
            options.gridConfig.uniformFraction = readDoubleOption(index, argc, argv, argument);
            continue;
        }
        if (argument == "--grid-refinement-score") {
            options.gridConfig.refinementScore =
                parseGridRefinementScore(
                    readStringOption(index, argc, argv, argument),
                    argument);
            continue;
        }
        if (argument == "--grid-variance-probe-bins") {
            options.gridConfig.varianceProbeBins =
                parsePositiveInt(readStringOption(index, argc, argv, argument), argument);
            continue;
        }

        if (!argument.empty() && argument[0] == '-') {
            throw std::invalid_argument("unknown option: " + argument);
        }

        if (options.outputPath != "weighted_sampler_study.csv") {
            throw std::invalid_argument("only one output CSV path can be specified");
        }
        options.outputPath = argument;
    }

    if (options.momenta.empty()) {
        throw std::invalid_argument("at least one --momentum value is required");
    }

    return options;
}

lowmc::PhotonPhaseSpace makeDiagnosticPhaseSpace(const Options& options) {
    return lowmc::PhotonPhaseSpace{
        options.kT,
        std::nextafter(options.kT, std::numeric_limits<double>::infinity()),
        options.etaMin,
        options.etaMax,
        options.phiMin,
        options.phiMax,
    };
}

lowmc::LowCurrent makeBackToBackToyCurrent(
    double momentum, double particleEta, double particlePhi) {
    const double mass = 0.13957039;
    const double pt = momentum / std::cosh(particleEta);
    const double px = pt * std::cos(particlePhi);
    const double py = pt * std::sin(particlePhi);
    const double pz = pt * std::sinh(particleEta);
    const double energy = std::sqrt(momentum * momentum + mass * mass);

    return lowmc::LowCurrent({
        {lowmc::FourVector{energy, px, py, pz}, +1.0, lowmc::InOutSign::Outgoing},
        {lowmc::FourVector{energy, -px, -py, -pz}, -1.0, lowmc::InOutSign::Outgoing},
    });
}

lowmc::WeightSummary runStudyPoint(
    const Options& options,
    const lowmc::PhotonPhaseSpace& phaseSpace,
    double momentum,
    std::uint64_t seed) {
    const lowmc::LowCurrent current =
        makeBackToBackToyCurrent(momentum, options.particleEta, options.particlePhi);
    std::mt19937_64 rng(seed);
    std::vector<lowmc::WeightedPhotonSample> samples;
    samples.reserve(static_cast<std::size_t>(options.samples));

    if (options.proposal == "grid") {
        const lowmc::GridAdaptiveEtaPhiProposal proposal =
            lowmc::makeGridAdaptiveEtaPhiProposal(
                current,
                phaseSpace,
                options.kT,
                options.gridConfig);
        for (int index = 0; index < options.samples; ++index) {
            const double unitEta = std::generate_canonical<double, 53>(rng);
            const double unitPhi = std::generate_canonical<double, 53>(rng);
            const double unitComponent = std::generate_canonical<double, 53>(rng);
            const double unitCell = std::generate_canonical<double, 53>(rng);
            const lowmc::GridAdaptiveEtaPhiProposalPoint etaPhiProposal =
                lowmc::sampleGridAdaptiveEtaPhi(
                    proposal,
                    unitComponent,
                    unitCell,
                    unitEta,
                    unitPhi);
            const lowmc::FourVector momentum =
                lowmc::photonMomentumFromKTEtaPhi(
                    options.kT,
                    etaPhiProposal.eta,
                    etaPhiProposal.phi);
            samples.push_back(lowmc::makeWeightedPhotonSample(
                current,
                phaseSpace,
                lowmc::WeightedPhotonProposal{momentum, etaPhiProposal.density}));
        }
        return lowmc::summarizeWeightedPhotonSamples(samples);
    }

    if (options.proposal == "leg") {
        const lowmc::LegCenteredEtaPhiProposal proposal =
            lowmc::makeLegCenteredEtaPhiProposal(
                phaseSpace,
                current.particles(),
                options.legConfig);
        for (int index = 0; index < options.samples; ++index) {
            const double unitComponent = std::generate_canonical<double, 53>(rng);
            const double unitPolar = std::generate_canonical<double, 53>(rng);
            const double unitAzimuth = std::generate_canonical<double, 53>(rng);
            samples.push_back(
                lowmc::makeFixedKTLegCenteredEtaPhiWeightedSample(
                    current,
                    options.kT,
                    proposal,
                    unitComponent,
                    unitPolar,
                    unitAzimuth)
                    .weightedPhoton);
        }
        return lowmc::summarizeWeightedPhotonSamples(samples);
    }

    for (int index = 0; index < options.samples; ++index) {
        const double unitEta = std::generate_canonical<double, 53>(rng);
        const double unitPhi = std::generate_canonical<double, 53>(rng);
        const double unitComponent = std::generate_canonical<double, 53>(rng);
        if (options.proposal == "adaptive") {
            const lowmc::FixedKTAdaptiveWeightedSample sample =
                lowmc::makeFixedKTAdaptiveEtaPhiWeightedSample(
                    current,
                    phaseSpace,
                    options.kT,
                    unitComponent,
                    unitEta,
                    unitPhi,
                    options.adaptiveConfig);
            samples.push_back(sample.weightedPhoton);
        } else if (options.proposal == "smooth") {
            const lowmc::FixedKTSmoothWeightedSample sample =
                lowmc::makeFixedKTSmoothEtaPhiWeightedSample(
                    current,
                    phaseSpace,
                    options.kT,
                    unitComponent,
                    unitEta,
                    unitPhi,
                    options.smoothConfig);
            samples.push_back(sample.weightedPhoton);
        } else {
            const lowmc::FixedKTWeightedSample sample =
                lowmc::makeFixedKTUniformEtaPhiWeightedSample(
                    current,
                    phaseSpace,
                    options.kT,
                    unitEta,
                    unitPhi);
            samples.push_back(sample.weightedPhoton);
        }
    }

    return lowmc::summarizeWeightedPhotonSamples(samples);
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

    const lowmc::PhotonPhaseSpace phaseSpace = makeDiagnosticPhaseSpace(options);
    if (!lowmc::validatePhotonPhaseSpace(phaseSpace).valid()) {
        std::cerr << "Invalid photon phase space\n";
        return EXIT_FAILURE;
    }

    std::ofstream output(options.outputPath);
    if (!output) {
        std::cerr << "Could not open output file: " << options.outputPath << '\n';
        return EXIT_FAILURE;
    }

    output << "proposal,adaptiveUniformFraction,adaptiveWidthScale,"
           << "adaptiveMinHalfWidth,adaptiveMaxHalfWidth,"
           << "smoothUniformFraction,smoothWidthScale,"
           << "smoothMinHalfWidth,smoothMaxHalfWidth,"
           << "legUniformFraction,legMinimumOneMinusBeta,"
           << "gridEtaBins,gridPhiBins,gridRefinementLevels,"
           << "gridRefineFraction,gridUniformFraction,"
           << "gridRefinementScore,gridVarianceProbeBins,"
           << "momentum,kT,samples,sumWeight,meanWeight,minWeight,maxWeight,"
           << "varianceWeight,effectiveSampleSize,effectiveSampleFraction\n";
    output << std::setprecision(17);

    for (std::size_t index = 0; index < options.momenta.size(); ++index) {
        const double momentum = options.momenta[index];
        const lowmc::WeightSummary summary =
            runStudyPoint(options, phaseSpace, momentum, options.seed + index);
        const double effectiveSampleFraction =
            summary.count > 0
                ? summary.effectiveSampleSize / static_cast<double>(summary.count)
                : 0.0;

        output << options.proposal << ','
               << options.adaptiveConfig.uniformFraction << ','
               << options.adaptiveConfig.widthScale << ','
               << options.adaptiveConfig.minHalfWidth << ','
               << options.adaptiveConfig.maxHalfWidth << ','
               << options.smoothConfig.uniformFraction << ','
               << options.smoothConfig.widthScale << ','
               << options.smoothConfig.minHalfWidth << ','
               << options.smoothConfig.maxHalfWidth << ','
               << options.legConfig.uniformFraction << ','
               << options.legConfig.minimumOneMinusBeta << ','
               << options.gridConfig.etaBins << ','
               << options.gridConfig.phiBins << ','
               << options.gridConfig.refinementLevels << ','
               << options.gridConfig.refineFraction << ','
               << options.gridConfig.uniformFraction << ','
               << gridRefinementScoreName(options.gridConfig.refinementScore) << ','
               << options.gridConfig.varianceProbeBins << ','
               << momentum << ','
               << options.kT << ','
               << summary.count << ','
               << summary.sumWeight << ','
               << summary.meanWeight << ','
               << summary.minWeight << ','
               << summary.maxWeight << ','
               << summary.varianceWeight << ','
               << summary.effectiveSampleSize << ','
               << effectiveSampleFraction << '\n';
    }

    if (!output) {
        std::cerr << "Failed while writing output file: " << options.outputPath << '\n';
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
