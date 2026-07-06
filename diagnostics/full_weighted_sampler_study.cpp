#include "ParticleInput.h"
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
    std::string outputPath = "full_weighted_sampler_study.csv";
    int samples = 1000;
    std::uint64_t seed = 12345;
    double kTMin = 0.002;
    double kTMax = 0.2;
    double particleEta = 0.0;
    double particlePhi = 0.5 * std::acos(-1.0);
    double etaMin = -5.0;
    double etaMax = 5.0;
    double phiMin = 0.0;
    double phiMax = 2.0 * std::acos(-1.0);
    std::vector<double> momenta = {100.0};
    std::string particlesPath;
    std::string proposal = "uniform";
    lowmc::GridAdaptiveEtaPhiProposalConfig gridConfig{
        32,
        32,
        3,
        0.1,
        0.1,
    };
};

void printUsage(const char* programName) {
    std::cout << "Usage: " << programName
              << " [output.csv] [--samples n] [--seed n]\n"
              << "       [--kt-min GeV] [--kt-max GeV]\n"
              << "       [--momentum GeV]... [--particle-eta eta] [--particle-phi phi]\n"
              << "       [--particles path]\n"
              << "       [--eta-min eta] [--eta-max eta] [--phi-min phi] [--phi-max phi]\n"
              << "       [--proposal uniform|grid]\n"
              << "       [--grid-eta-bins n] [--grid-phi-bins n]\n"
              << "       [--grid-refinement-levels n] [--grid-refine-fraction f]\n"
              << "       [--grid-uniform-fraction f]\n\n"
              << "Write full phase-space weighted sampler summaries for a two-pion toy event.\n";
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
        if (argument == "--momentum") {
            if (!customMomenta) {
                options.momenta.clear();
                customMomenta = true;
            }
            options.momenta.push_back(
                parsePositiveDouble(readStringOption(index, argc, argv, argument), argument));
            continue;
        }
        if (argument == "--particles") {
            options.particlesPath = readStringOption(index, argc, argv, argument);
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
            if (options.proposal != "uniform" && options.proposal != "grid") {
                throw std::invalid_argument("--proposal must be uniform or grid");
            }
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

        if (!argument.empty() && argument[0] == '-') {
            throw std::invalid_argument("unknown option: " + argument);
        }

        if (options.outputPath != "full_weighted_sampler_study.csv") {
            throw std::invalid_argument("only one output CSV path can be specified");
        }
        options.outputPath = argument;
    }

    if (options.momenta.empty()) {
        throw std::invalid_argument("at least one --momentum value is required");
    }
    if (options.kTMax <= options.kTMin) {
        throw std::invalid_argument("--kt-max must be greater than --kt-min");
    }

    return options;
}

lowmc::PhotonPhaseSpace makeDiagnosticPhaseSpace(const Options& options) {
    return lowmc::PhotonPhaseSpace{
        options.kTMin,
        options.kTMax,
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
    const lowmc::LowCurrent& current,
    std::uint64_t seed,
    double& gridBuildKT) {
    std::mt19937_64 rng(seed);
    std::vector<lowmc::WeightedPhotonSample> samples;
    samples.reserve(static_cast<std::size_t>(options.samples));

    gridBuildKT = std::sqrt(phaseSpace.kTMin * phaseSpace.kTMax);
    lowmc::GridAdaptiveEtaPhiProposal gridProposal;
    if (options.proposal == "grid") {
        gridProposal = lowmc::makeGridAdaptiveEtaPhiProposal(
            current,
            phaseSpace,
            gridBuildKT,
            options.gridConfig);
    }

    for (int index = 0; index < options.samples; ++index) {
        const double unitKT = std::generate_canonical<double, 53>(rng);
        const double unitEta = std::generate_canonical<double, 53>(rng);
        const double unitPhi = std::generate_canonical<double, 53>(rng);
        if (options.proposal == "grid") {
            const double unitComponent = std::generate_canonical<double, 53>(rng);
            const double unitCell = std::generate_canonical<double, 53>(rng);
            const lowmc::GridAdaptiveWeightedPhotonSample sample =
                lowmc::makeGridAdaptiveWeightedPhotonSample(
                    current,
                    phaseSpace,
                    gridProposal,
                    unitKT,
                    unitComponent,
                    unitCell,
                    unitEta,
                    unitPhi);
            samples.push_back(sample.weightedPhoton);
        } else {
            const lowmc::MinimalWeightedPhotonSample sample =
                lowmc::makeMinimalWeightedPhotonSample(
                    current,
                    phaseSpace,
                    unitKT,
                    unitEta,
                    unitPhi);
            samples.push_back(sample.weightedPhoton);
        }
    }

    return lowmc::summarizeWeightedPhotonSamples(samples);
}

lowmc::WeightSummary runToyStudyPoint(
    const Options& options,
    const lowmc::PhotonPhaseSpace& phaseSpace,
    double momentum,
    std::uint64_t seed,
    double& gridBuildKT) {
    const lowmc::LowCurrent current =
        makeBackToBackToyCurrent(momentum, options.particleEta, options.particlePhi);
    return runStudyPoint(options, phaseSpace, current, seed, gridBuildKT);
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

    output << "proposal,gridEtaBins,gridPhiBins,gridRefinementLevels,"
           << "gridRefineFraction,gridUniformFraction,"
           << "momentum,kTMin,kTMax,gridBuildKT,samples,sumWeight,meanWeight,"
           << "minWeight,maxWeight,varianceWeight,effectiveSampleSize,"
           << "effectiveSampleFraction\n";
    output << std::setprecision(17);

    if (!options.particlesPath.empty()) {
        double gridBuildKT = std::numeric_limits<double>::quiet_NaN();
        const lowmc::LowCurrent current(
            lowmc_diagnostics::readChargedParticlesFromTextFile(options.particlesPath));
        const lowmc::WeightSummary summary =
            runStudyPoint(options, phaseSpace, current, options.seed, gridBuildKT);
        const double effectiveSampleFraction =
            summary.count > 0
                ? summary.effectiveSampleSize / static_cast<double>(summary.count)
                : 0.0;

        output << options.proposal << ','
               << options.gridConfig.etaBins << ','
               << options.gridConfig.phiBins << ','
               << options.gridConfig.refinementLevels << ','
               << options.gridConfig.refineFraction << ','
               << options.gridConfig.uniformFraction << ','
               << std::numeric_limits<double>::quiet_NaN() << ','
               << options.kTMin << ','
               << options.kTMax << ','
               << gridBuildKT << ','
               << summary.count << ','
               << summary.sumWeight << ','
               << summary.meanWeight << ','
               << summary.minWeight << ','
               << summary.maxWeight << ','
               << summary.varianceWeight << ','
               << summary.effectiveSampleSize << ','
               << effectiveSampleFraction << '\n';
        return EXIT_SUCCESS;
    }

    for (std::size_t index = 0; index < options.momenta.size(); ++index) {
        const double momentum = options.momenta[index];
        double gridBuildKT = std::numeric_limits<double>::quiet_NaN();
        const lowmc::WeightSummary summary =
            runToyStudyPoint(options, phaseSpace, momentum, options.seed + index, gridBuildKT);
        const double effectiveSampleFraction =
            summary.count > 0
                ? summary.effectiveSampleSize / static_cast<double>(summary.count)
                : 0.0;

        output << options.proposal << ','
               << options.gridConfig.etaBins << ','
               << options.gridConfig.phiBins << ','
               << options.gridConfig.refinementLevels << ','
               << options.gridConfig.refineFraction << ','
               << options.gridConfig.uniformFraction << ','
               << momentum << ','
               << options.kTMin << ','
               << options.kTMax << ','
               << gridBuildKT << ','
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
