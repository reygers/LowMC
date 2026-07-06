#include "lowmc/LowMC.h"

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

struct Options {
    std::string outputPath = "weighted_event_example.csv";
    int samples = 1000;
    int oversamplingFactor = 1;
    std::uint64_t seed = 12345;
    double momentum = 100.0;
    double kTMin = 0.002;
    double kTMax = 0.2;
};

void printUsage(const char* programName) {
    std::cout << "Usage: " << programName
              << " [output.csv] [--samples n] [--oversampling-factor n] [--seed n]\n"
              << "       [--momentum GeV] [--kt-min GeV] [--kt-max GeV]\n"
              << "Generate weighted photons for a two-pion toy event with the event-level API.\n";
}

int parsePositiveInt(const std::string& text, const std::string& optionName) {
    std::size_t parsed = 0;
    const int value = std::stoi(text, &parsed);
    if (parsed != text.size() || value <= 0) {
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

double parsePositiveDouble(const std::string& text, const std::string& optionName) {
    std::size_t parsed = 0;
    const double value = std::stod(text, &parsed);
    if (parsed != text.size() || value <= 0.0) {
        throw std::invalid_argument(optionName + " must be a positive number");
    }
    return value;
}

std::string readStringOption(int& index, int argc, char** argv, const std::string& optionName) {
    if (index + 1 >= argc) {
        throw std::invalid_argument(optionName + " requires a value");
    }
    return argv[++index];
}

Options parseOptions(int argc, char** argv) {
    Options options;
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
        if (argument == "--oversampling-factor") {
            options.oversamplingFactor =
                parsePositiveInt(readStringOption(index, argc, argv, argument), argument);
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
        if (!argument.empty() && argument[0] == '-') {
            throw std::invalid_argument("unknown option: " + argument);
        }
        if (options.outputPath != "weighted_event_example.csv") {
            throw std::invalid_argument("only one output CSV path can be specified");
        }
        options.outputPath = argument;
    }
    if (options.kTMax <= options.kTMin) {
        throw std::invalid_argument("--kt-max must be greater than --kt-min");
    }
    return options;
}

std::vector<lowmc::ChargedParticle> makeBackToBackPions(double momentum) {
    const double mass = 0.13957039;
    const double energy = std::sqrt(momentum * momentum + mass * mass);
    return {
        {lowmc::FourVector{energy, 0.0, momentum, 0.0}, +1.0, lowmc::InOutSign::Outgoing},
        {lowmc::FourVector{energy, 0.0, -momentum, 0.0}, -1.0, lowmc::InOutSign::Outgoing},
    };
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

    const lowmc::WeightedEventGeneratorConfig config{
        lowmc::PhotonPhaseSpace{
            options.kTMin,
            options.kTMax,
            -5.0,
            5.0,
            0.0,
            2.0 * std::acos(-1.0),
        },
        options.samples,
        options.seed,
        lowmc::tunedMassGridAdaptiveEtaPhiConfig(),
        options.oversamplingFactor,
    };
    const lowmc::WeightedEventGenerationResult result =
        lowmc::generateWeightedPhotons(makeBackToBackPions(options.momentum), config);

    std::ofstream output(options.outputPath);
    if (!output) {
        std::cerr << "Could not open output file: " << options.outputPath << '\n';
        return EXIT_FAILURE;
    }

    output << "index,E,px,py,pz,kT,eta,phi,proposalDensity,weight\n";
    output << std::setprecision(17);
    for (std::size_t index = 0; index < result.samples.size(); ++index) {
        const lowmc::WeightedEventPhotonSample& sample = result.samples[index];
        const lowmc::Photon& photon = sample.photon;
        output << index << ','
               << photon.momentum.E << ','
               << photon.momentum.px << ','
               << photon.momentum.py << ','
               << photon.momentum.pz << ','
               << sample.kT << ','
               << sample.eta << ','
               << sample.phi << ','
               << sample.proposalDensity << ','
               << photon.weight << '\n';
    }

    if (!output) {
        std::cerr << "Failed while writing output file: " << options.outputPath << '\n';
        return EXIT_FAILURE;
    }

    std::cerr << "samples=" << result.weightSummary.count
              << " baseSamples=" << options.samples
              << " oversamplingFactor=" << options.oversamplingFactor
              << " meanWeight=" << result.weightSummary.meanWeight
              << " effectiveSampleSize=" << result.weightSummary.effectiveSampleSize
              << " gridBuildKT=" << result.gridBuildKT
              << '\n';
    return EXIT_SUCCESS;
}
