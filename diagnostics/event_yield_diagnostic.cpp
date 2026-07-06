#include "ParticleInput.h"
#include "lowmc/LowMC.h"

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

struct Options {
    std::string outputPath = "event_yield_diagnostic.csv";
    std::string particlesPath;
    double momentum = 10.0;
    double kTMin = 0.002;
    double kTMax = 0.2;
    double etaMin = -5.0;
    double etaMax = 5.0;
    double phiMin = 0.0;
    double phiMax = 2.0 * std::acos(-1.0);
    int logKTBins = 16;
    int etaBins = 40;
    int phiBins = 48;
    int adaptiveEtaSubBins = 2;
    int adaptivePhiSubBins = 2;
    lowmc::GridAdaptiveEtaPhiProposalConfig adaptiveGridConfig =
        lowmc::tunedMassGridAdaptiveEtaPhiConfig();
    std::string method = "midpoint";
    double oversamplingFactor = 1.0;
    std::uint64_t seed = 12345;
};

void printUsage(const char* programName) {
    std::cout
        << "Usage: " << programName << " [output.csv] [options]\n"
        << "Estimate lambda_event and draw a diagnostic Poisson photon count.\n\n"
        << "Options:\n"
        << "  --particles path              Read charged particles from a text fixture\n"
        << "  --momentum GeV                Two-pion toy momentum if no fixture is given\n"
        << "  --kt-min GeV                  Photon kT lower bound\n"
        << "  --kt-max GeV                  Photon kT upper bound\n"
        << "  --eta-min value               Photon eta lower bound\n"
        << "  --eta-max value               Photon eta upper bound\n"
        << "  --phi-min value               Photon phi lower bound\n"
        << "  --phi-max value               Photon phi upper bound\n"
        << "  --log-kt-bins N               Midpoint bins in log(kT)\n"
        << "  --eta-bins N                  Midpoint bins in eta\n"
        << "  --phi-bins N                  Midpoint bins in phi\n"
        << "  --method midpoint|adaptive    Yield integration method\n"
        << "  --adaptive-eta-sub-bins N     Eta sub-bins per adaptive cell\n"
        << "  --adaptive-phi-sub-bins N     Phi sub-bins per adaptive cell\n"
        << "  --adaptive-grid-eta-bins N    Initial adaptive eta bins\n"
        << "  --adaptive-grid-phi-bins N    Initial adaptive phi bins\n"
        << "  --adaptive-refinement-levels N Adaptive refinement levels\n"
        << "  --adaptive-refine-fraction F  Fraction of cells refined per level\n"
        << "  --oversampling-factor F       Diagnostic Poisson oversampling factor\n"
        << "  --seed N                      Multiplicity RNG seed\n";
}

std::string readStringOption(int& index, int argc, char** argv, const std::string& optionName) {
    if (index + 1 >= argc) {
        throw std::invalid_argument(optionName + " requires a value");
    }
    return argv[++index];
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
        if (argument == "--log-kt-bins") {
            options.logKTBins =
                parsePositiveInt(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--eta-bins") {
            options.etaBins =
                parsePositiveInt(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--phi-bins") {
            options.phiBins =
                parsePositiveInt(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--method") {
            options.method = readStringOption(index, argc, argv, argument);
            if (options.method != "midpoint" && options.method != "adaptive") {
                throw std::invalid_argument("--method must be midpoint or adaptive");
            }
            continue;
        }
        if (argument == "--adaptive-eta-sub-bins") {
            options.adaptiveEtaSubBins =
                parsePositiveInt(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--adaptive-phi-sub-bins") {
            options.adaptivePhiSubBins =
                parsePositiveInt(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--adaptive-grid-eta-bins") {
            options.adaptiveGridConfig.etaBins =
                parsePositiveInt(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--adaptive-grid-phi-bins") {
            options.adaptiveGridConfig.phiBins =
                parsePositiveInt(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--adaptive-refinement-levels") {
            options.adaptiveGridConfig.refinementLevels =
                parsePositiveInt(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--adaptive-refine-fraction") {
            options.adaptiveGridConfig.refineFraction =
                parsePositiveDouble(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--oversampling-factor") {
            options.oversamplingFactor =
                parsePositiveDouble(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--seed") {
            options.seed = parseUInt64(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (!argument.empty() && argument[0] == '-') {
            throw std::invalid_argument("unknown option: " + argument);
        }
        if (options.outputPath != "event_yield_diagnostic.csv") {
            throw std::invalid_argument("only one output CSV path can be specified");
        }
        options.outputPath = argument;
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

    try {
        const std::vector<lowmc::ChargedParticle> particles =
            options.particlesPath.empty()
                ? makeBackToBackPions(options.momentum)
                : lowmc_diagnostics::readChargedParticlesFromTextFile(options.particlesPath);

        const lowmc::PhotonPhaseSpace phaseSpace{
            options.kTMin,
            options.kTMax,
            options.etaMin,
            options.etaMax,
            options.phiMin,
            options.phiMax,
        };
        lowmc::EventYieldEstimate yield;
        if (options.method == "adaptive") {
            const lowmc::AdaptiveEventYieldIntegrationConfig integrationConfig{
                options.logKTBins,
                options.adaptiveEtaSubBins,
                options.adaptivePhiSubBins,
                options.adaptiveGridConfig,
            };
            yield = lowmc::estimateAdaptiveEventYield(particles, phaseSpace, integrationConfig);
        } else {
            const lowmc::EventYieldIntegrationConfig integrationConfig{
                options.logKTBins,
                options.etaBins,
                options.phiBins,
            };
            yield = lowmc::estimateEventYield(particles, phaseSpace, integrationConfig);
        }
        const lowmc::PhotonMultiplicityDiagnostic multiplicity =
            lowmc::samplePhotonMultiplicity(
                yield.lambdaEvent,
                lowmc::PhotonMultiplicityConfig{options.oversamplingFactor, options.seed});
        const lowmc::ChargeConservationDiagnostic charge =
            lowmc::LowCurrent(particles).chargeConservation();

        std::ofstream output(options.outputPath);
        if (!output) {
            std::cerr << "Could not open output file: " << options.outputPath << '\n';
            return EXIT_FAILURE;
        }

        output << "method,status,particle_count,signed_charge_sum,charge_conserved,"
               << "kt_min,kt_max,eta_min,eta_max,phi_min,phi_max,"
               << "log_kt_bins,eta_bins,phi_bins,evaluations,grid_cells,"
               << "lambda_event,coarse_lambda_event,absolute_difference,relative_difference,"
               << "seed,oversampling_factor,poisson_mean,sampled_photons,sample_weight,"
               << "weighted_count\n";
        output << std::setprecision(17)
               << yield.method << ','
               << yield.status << ','
               << particles.size() << ','
               << charge.netCharge << ','
               << (charge.conserved ? 1 : 0) << ','
               << phaseSpace.kTMin << ','
               << phaseSpace.kTMax << ','
               << phaseSpace.etaMin << ','
               << phaseSpace.etaMax << ','
               << phaseSpace.phiMin << ','
               << phaseSpace.phiMax << ','
               << yield.logKTBins << ','
               << yield.etaBins << ','
               << yield.phiBins << ','
               << yield.evaluations << ','
               << yield.gridCells << ','
               << yield.lambdaEvent << ','
               << yield.coarseLambdaEvent << ','
               << yield.absoluteDifference << ','
               << yield.relativeDifference << ','
               << multiplicity.seed << ','
               << multiplicity.oversamplingFactor << ','
               << multiplicity.poissonMean << ','
               << multiplicity.sampledPhotons << ','
               << multiplicity.sampleWeight << ','
               << multiplicity.expectedWeightedCount << '\n';

        if (!output) {
            std::cerr << "Failed while writing output file: " << options.outputPath << '\n';
            return EXIT_FAILURE;
        }

        std::cerr << "method=" << yield.method
                  << " status=" << yield.status
                  << " particle_count=" << particles.size()
                  << " signed_charge_sum=" << charge.netCharge
                  << " lambda_event=" << yield.lambdaEvent
                  << " relative_difference=" << yield.relativeDifference
                  << " evaluations=" << yield.evaluations
                  << " grid_cells=" << yield.gridCells
                  << " poisson_mean=" << multiplicity.poissonMean
                  << " sampled_photons=" << multiplicity.sampledPhotons
                  << " sample_weight=" << multiplicity.sampleWeight
                  << '\n';
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
