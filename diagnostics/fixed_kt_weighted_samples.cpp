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
    std::string outputPath = "fixed_kt_weighted_samples.csv";
    int samples = 1000;
    std::uint64_t seed = 12345;
    double kT = 0.02;
    double momentum = 1.0;
    double particleEta = 0.0;
    double particlePhi = 0.5 * std::acos(-1.0);
    double etaMin = -5.0;
    double etaMax = 5.0;
    double phiMin = 0.0;
    double phiMax = 2.0 * std::acos(-1.0);
    std::string proposal = "uniform";
    lowmc::AdaptiveEtaPhiProposalConfig adaptiveConfig;
    lowmc::SmoothEtaPhiProposalConfig smoothConfig;
    lowmc::GridAdaptiveEtaPhiProposalConfig gridConfig = lowmc::tunedMassGridAdaptiveEtaPhiConfig();
};

void printUsage(const char* programName) {
    std::cout << "Usage: " << programName
              << " [output.csv] [--samples n] [--seed n] [--kt GeV]\n"
              << "       [--momentum GeV] [--particle-eta eta] [--particle-phi phi]\n"
              << "       [--eta-min eta] [--eta-max eta] [--phi-min phi] [--phi-max phi]\n"
              << "       [--proposal uniform|adaptive|smooth|grid]\n\n"
              << "       [--adaptive-uniform-fraction f] [--adaptive-width-scale s]\n"
              << "       [--adaptive-min-half-width w] [--adaptive-max-half-width w]\n\n"
              << "       [--smooth-uniform-fraction f] [--smooth-width-scale s]\n"
              << "       [--smooth-min-half-width w] [--smooth-max-half-width w]\n\n"
              << "       [--grid-eta-bins n] [--grid-phi-bins n]\n"
              << "       [--grid-refinement-levels n] [--grid-refine-fraction f]\n"
              << "       [--grid-uniform-fraction f]\n\n"
              << "Write weighted fixed-kT photon samples for a back-to-back charged-pion toy event.\n"
              << "The default eta-phi proposal is uniform and diagnostic-only.\n";
}

double parseDouble(const std::string& text, const std::string& optionName) {
    std::size_t parsed = 0;
    const double value = std::stod(text, &parsed);
    if (parsed != text.size()) {
        throw std::invalid_argument(optionName + " must be a number");
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

double parsePositiveDouble(const std::string& text, const std::string& optionName) {
    const double value = parseDouble(text, optionName);
    if (value <= 0.0) {
        throw std::invalid_argument(optionName + " must be positive");
    }
    return value;
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
            options.momentum =
                parsePositiveDouble(readStringOption(index, argc, argv, argument), argument);
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
                options.proposal != "grid") {
                throw std::invalid_argument(
                    "--proposal must be uniform, adaptive, smooth, or grid");
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

        if (options.outputPath != "fixed_kt_weighted_samples.csv") {
            throw std::invalid_argument("only one output CSV path can be specified");
        }
        options.outputPath = argument;
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

    const lowmc::LowCurrent current =
        makeBackToBackToyCurrent(options.momentum, options.particleEta, options.particlePhi);
    lowmc::GridAdaptiveEtaPhiProposal gridProposal;
    if (options.proposal == "grid") {
        gridProposal =
            lowmc::makeGridAdaptiveEtaPhiProposal(
                current,
                phaseSpace,
                options.kT,
                options.gridConfig);
    }
    std::mt19937_64 rng(options.seed);
    std::vector<lowmc::WeightedPhotonSample> samples;
    samples.reserve(static_cast<std::size_t>(options.samples));

    output << "index,kT,eta,phi,E,px,py,pz,proposalDensity,weight\n";
    output << std::setprecision(17);

    for (int index = 0; index < options.samples; ++index) {
        const double unitEta = std::generate_canonical<double, 53>(rng);
        const double unitPhi = std::generate_canonical<double, 53>(rng);
        const double unitComponent = std::generate_canonical<double, 53>(rng);
        const double unitCell = std::generate_canonical<double, 53>(rng);
        lowmc::FourVector photon;
        double proposalDensity = 0.0;
        double weight = 0.0;
        double eta = 0.0;
        double phi = 0.0;

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
            photon = sample.weightedPhoton.photon.momentum;
            proposalDensity = sample.weightedPhoton.proposalDensity;
            weight = sample.weightedPhoton.photon.weight;
            eta = sample.etaPhiProposal.eta;
            phi = sample.etaPhiProposal.phi;
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
            photon = sample.weightedPhoton.photon.momentum;
            proposalDensity = sample.weightedPhoton.proposalDensity;
            weight = sample.weightedPhoton.photon.weight;
            eta = sample.etaPhiProposal.eta;
            phi = sample.etaPhiProposal.phi;
        } else if (options.proposal == "grid") {
            const lowmc::GridAdaptiveEtaPhiProposalPoint etaPhiProposal =
                lowmc::sampleGridAdaptiveEtaPhi(
                    gridProposal,
                    unitComponent,
                    unitCell,
                    unitEta,
                    unitPhi);
            photon = lowmc::photonMomentumFromKTEtaPhi(
                options.kT,
                etaPhiProposal.eta,
                etaPhiProposal.phi);
            const lowmc::WeightedPhotonSample sample =
                lowmc::makeWeightedPhotonSample(
                    current,
                    phaseSpace,
                    lowmc::WeightedPhotonProposal{photon, etaPhiProposal.density});
            samples.push_back(sample);
            proposalDensity = sample.proposalDensity;
            weight = sample.photon.weight;
            eta = etaPhiProposal.eta;
            phi = etaPhiProposal.phi;
        } else {
            const lowmc::FixedKTWeightedSample sample =
                lowmc::makeFixedKTUniformEtaPhiWeightedSample(
                    current,
                    phaseSpace,
                    options.kT,
                    unitEta,
                    unitPhi);
            samples.push_back(sample.weightedPhoton);
            photon = sample.weightedPhoton.photon.momentum;
            proposalDensity = sample.weightedPhoton.proposalDensity;
            weight = sample.weightedPhoton.photon.weight;
            eta = sample.etaPhiProposal.eta;
            phi = sample.etaPhiProposal.phi;
        }

        output << index << ','
               << options.kT << ','
               << eta << ','
               << phi << ','
               << photon.E << ','
               << photon.px << ','
               << photon.py << ','
               << photon.pz << ','
               << proposalDensity << ','
               << weight << '\n';
    }

    if (!output) {
        std::cerr << "Failed while writing output file: " << options.outputPath << '\n';
        return EXIT_FAILURE;
    }

    const lowmc::WeightSummary summary = lowmc::summarizeWeightedPhotonSamples(samples);
    std::cerr << "samples=" << summary.count
              << " sumWeight=" << summary.sumWeight
              << " meanWeight=" << summary.meanWeight
              << " minWeight=" << summary.minWeight
              << " maxWeight=" << summary.maxWeight
              << " effectiveSampleSize=" << summary.effectiveSampleSize
              << '\n';

    return EXIT_SUCCESS;
}
