#include "lowmc/LowMC.h"
#include "lowmc/pythia8/Pythia8Adapter.h"

#include <Pythia8/Pythia.h>

#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

#ifndef LOWMC_PYTHIA8_XMLDOC
#define LOWMC_PYTHIA8_XMLDOC "../share/Pythia8/xmldoc"
#endif

namespace {

struct Options {
    std::string outputPath = "pythia8_yield_convergence.csv";
    int events = 10;
    std::uint64_t seed = 35791;
};

std::string readOption(
    int& index, int argc, char** argv, const std::string& option) {
    if (index + 1 >= argc) {
        throw std::invalid_argument(option + " requires a value");
    }
    return argv[++index];
}

int parsePositiveInt(const std::string& text, const std::string& option) {
    std::size_t parsed = 0;
    const long value = std::stol(text, &parsed);
    if (parsed != text.size() || value <= 0 ||
        value > std::numeric_limits<int>::max()) {
        throw std::invalid_argument(option + " expects a positive integer");
    }
    return static_cast<int>(value);
}

std::uint64_t parseSeed(const std::string& text) {
    std::size_t parsed = 0;
    const unsigned long long value = std::stoull(text, &parsed);
    if (parsed != text.size() || value == 0 || value > 900000000ULL) {
        throw std::invalid_argument("--seed expects an integer in [1,900000000]");
    }
    return static_cast<std::uint64_t>(value);
}

Options parseOptions(int argc, char** argv) {
    Options options;
    bool outputSeen = false;
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument == "--help" || argument == "-h") {
            std::cout
                << "Usage: lowmc_pythia8_yield_convergence [output.csv] "
                   "[--events N] [--seed N]\n";
            std::exit(EXIT_SUCCESS);
        }
        if (argument == "--events") {
            options.events = parsePositiveInt(
                readOption(index, argc, argv, argument), argument);
        } else if (argument == "--seed") {
            options.seed = parseSeed(readOption(index, argc, argv, argument));
        } else if (!argument.empty() && argument[0] == '-') {
            throw std::invalid_argument("unknown option: " + argument);
        } else if (!outputSeen) {
            options.outputPath = argument;
            outputSeen = true;
        } else {
            throw std::invalid_argument("only one output path may be specified");
        }
    }
    return options;
}

void initializePythia(Pythia8::Pythia& pythia, const Options& options) {
    pythia.readString("Print:quiet = on");
    pythia.readString("Beams:idA = 2212");
    pythia.readString("Beams:idB = 2212");
    pythia.readString("Beams:eCM = 13000.");
    pythia.readString("Tune:pp = 14");
    pythia.readString("SoftQCD:nonDiffractive = on");
    pythia.readString("Random:setSeed = on");
    pythia.readString("Random:seed = " + std::to_string(options.seed));
    pythia.readString("TimeShower:QEDshowerByQ = off");
    pythia.readString("TimeShower:QEDshowerByL = off");
    pythia.readString("TimeShower:QEDshowerByOther = off");
    pythia.readString("TimeShower:QEDshowerByGamma = off");
    pythia.readString("SpaceShower:QEDshowerByQ = off");
    pythia.readString("SpaceShower:QEDshowerByL = off");
    pythia.readString("HadronLevel:QED = off");
    if (!pythia.init()) {
        throw std::runtime_error("Pythia initialization failed");
    }
}

lowmc::EventYieldEstimate estimate(
    const std::vector<lowmc::ChargedParticle>& particles,
    int level) {
    lowmc::GridAdaptiveEtaPhiProposalConfig grid =
        lowmc::tunedMassGridAdaptiveEtaPhiConfig();
    grid.etaBins = 32;
    grid.phiBins = 32;
    grid.refinementLevels = level;
    grid.refineFraction = 0.1;
    const lowmc::AdaptiveEventYieldIntegrationConfig config{16, 2, 2, grid};
    const lowmc::PhotonPhaseSpace phaseSpace{
        0.002, 0.2, -5.0, 5.0, 0.0, 2.0 * std::acos(-1.0)};
    return lowmc::estimateFactorizedAdaptiveEventYield(
        particles, phaseSpace, config);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Options options = parseOptions(argc, argv);
        Pythia8::Pythia pythia(LOWMC_PYTHIA8_XMLDOC, false);
        initializePythia(pythia, options);
        std::ofstream output(options.outputPath);
        if (!output) {
            throw std::runtime_error("could not open output CSV: " + options.outputPath);
        }
        output << std::setprecision(17);
        output << "seed,event_index,attempted_events,selected_charged,"
                  "signed_charge_sum,lambda_l7,delta_l7,lambda_l9,delta_l9,"
                  "lambda_l7_over_l9\n";

        int accepted = 0;
        int attempted = 0;
        while (accepted < options.events) {
            ++attempted;
            if (!pythia.next()) {
                if (attempted > 20 * options.events) {
                    throw std::runtime_error("too many Pythia generation failures");
                }
                continue;
            }
            const lowmc::AdapterParticleSelectionResult selection =
                lowmc::pythia8::selectChargedParticles(pythia.event);
            if (std::abs(selection.diagnostic.signedChargeSum) > 1.0e-12) {
                continue;
            }
            const lowmc::EventYieldEstimate level7 =
                estimate(selection.selectedParticles, 7);
            const lowmc::EventYieldEstimate level9 =
                estimate(selection.selectedParticles, 9);
            output << options.seed << ',' << accepted << ',' << attempted << ','
                   << selection.selectedParticles.size() << ','
                   << selection.diagnostic.signedChargeSum << ','
                   << level7.lambdaEvent << ',' << level7.relativeDifference << ','
                   << level9.lambdaEvent << ',' << level9.relativeDifference << ','
                   << level7.lambdaEvent / level9.lambdaEvent << '\n';
            ++accepted;
        }
        std::cout << "pythia8_yield_convergence=ok"
                  << " seed=" << options.seed
                  << " accepted_events=" << accepted
                  << " attempted_events=" << attempted
                  << " output=" << options.outputPath << '\n';
    } catch (const std::exception& error) {
        std::cerr << "FAILED: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
