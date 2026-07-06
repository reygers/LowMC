#include "lowmc/hepmc3/HepMC3Adapter.h"

#include <HepMC3/GenEvent.h>
#include <HepMC3/ReaderAscii.h>
#include <HepMC3/Units.h>

#include <cmath>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {

struct Options {
    std::string inputPath;
    std::string outputPath;
    int eventIndex = 0;
};

void printUsage(std::ostream& out) {
    out << "Usage: lowmc_hepmc3_selected_particle_text [--event-index N] "
        << "<input.hepmc3> <output.lowmc>\n"
        << "\n"
        << "Reads one HepMC3 event, runs the LowMC HepMC3 charged-particle selection,\n"
        << "and writes selected charged particles as diagnostic text columns:\n"
        << "charge inOutSign mass p eta phi.\n"
        << "No photons are generated or written to the HepMC3 event record.\n";
}

int parseNonNegativeInt(const std::string& value, const std::string& optionName) {
    std::size_t parsed = 0;
    const long parsedValue = std::stol(value, &parsed);
    if (parsed != value.size() ||
        parsedValue < 0 ||
        parsedValue > std::numeric_limits<int>::max()) {
        throw std::invalid_argument(optionName + " expects a non-negative integer");
    }
    return static_cast<int>(parsedValue);
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
            printUsage(std::cout);
            std::exit(EXIT_SUCCESS);
        }
        if (argument == "--event-index") {
            options.eventIndex =
                parseNonNegativeInt(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (!argument.empty() && argument.front() == '-') {
            throw std::invalid_argument("unknown option: " + argument);
        }
        if (options.inputPath.empty()) {
            options.inputPath = argument;
            continue;
        }
        if (options.outputPath.empty()) {
            options.outputPath = argument;
            continue;
        }
        throw std::invalid_argument("expected exactly one input and one output file path");
    }

    if (options.inputPath.empty() || options.outputPath.empty()) {
        throw std::invalid_argument("missing input or output file path");
    }
    return options;
}

[[nodiscard]] double momentumMagnitude(const lowmc::FourVector& momentum) {
    return std::sqrt(
        momentum.px * momentum.px +
        momentum.py * momentum.py +
        momentum.pz * momentum.pz);
}

[[nodiscard]] double invariantMass(const lowmc::FourVector& momentum) {
    const double mass2 =
        momentum.E * momentum.E -
        momentum.px * momentum.px -
        momentum.py * momentum.py -
        momentum.pz * momentum.pz;
    return std::sqrt(std::max(0.0, mass2));
}

[[nodiscard]] double pseudorapidity(const lowmc::FourVector& momentum) {
    const double pT = std::hypot(momentum.px, momentum.py);
    if (pT == 0.0) {
        throw std::runtime_error(
            "selected particle has zero transverse momentum and cannot be "
            "represented in the diagnostic p eta phi text format");
    }
    return std::asinh(momentum.pz / pT);
}

[[nodiscard]] const char* signName(lowmc::InOutSign sign) {
    return sign == lowmc::InOutSign::Incoming ? "incoming" : "outgoing";
}

void writeSelectedParticles(
    const std::string& path,
    const lowmc::AdapterParticleSelectionResult& result) {
    std::ofstream output(path);
    if (!output) {
        throw std::runtime_error("could not open output file: " + path);
    }

    output << "# LowMC diagnostic particle input\n"
           << "# Generated from HepMC3 adapter selection.\n"
           << "# columns: charge inOutSign mass p eta phi\n";
    output << std::setprecision(17);
    for (const lowmc::ChargedParticle& particle : result.selectedParticles) {
        const double phi = std::atan2(particle.momentum.py, particle.momentum.px);
        output << particle.charge << ' '
               << signName(particle.inOutSign) << ' '
               << invariantMass(particle.momentum) << ' '
               << momentumMagnitude(particle.momentum) << ' '
               << pseudorapidity(particle.momentum) << ' '
               << phi << '\n';
    }

    if (!output) {
        throw std::runtime_error("failed while writing output file: " + path);
    }
}

void printDiagnostic(
    int eventIndex,
    HepMC3::Units::MomentumUnit momentumUnit,
    const lowmc::AdapterParticleSelectionResult& result) {
    const lowmc::AdapterSelectionDiagnostic& diagnostic = result.diagnostic;
    std::cout << "event_index=" << eventIndex
              << " momentum_unit=" << HepMC3::Units::name(momentumUnit)
              << " selected_incoming_charged=" << diagnostic.selectedIncomingCharged
              << " selected_outgoing_charged=" << diagnostic.selectedOutgoingCharged
              << " ignored_neutral_candidates=" << diagnostic.ignoredNeutralCandidates
              << " rejected_charged_candidates=" << diagnostic.rejectedChargedCount()
              << " signed_charge_sum=" << diagnostic.signedChargeSum
              << " ordinary_charge_sum=" << diagnostic.ordinaryChargeSum
              << " output=" << result.selectedParticles.size()
              << '\n';
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Options options = parseOptions(argc, argv);
        HepMC3::ReaderAscii reader(options.inputPath);
        if (reader.failed()) {
            std::cerr << "FAILED: could not open HepMC3 input file: "
                      << options.inputPath << '\n';
            return EXIT_FAILURE;
        }

        int eventsRead = 0;
        HepMC3::GenEvent event;
        while (eventsRead <= options.eventIndex && !reader.failed()) {
            if (!reader.read_event(event)) {
                break;
            }
            if (eventsRead == options.eventIndex) {
                const lowmc::AdapterParticleSelectionResult result =
                    lowmc::hepmc3::selectChargedParticles(event);
                writeSelectedParticles(options.outputPath, result);
                printDiagnostic(eventsRead, event.momentum_unit(), result);
                return EXIT_SUCCESS;
            }
            ++eventsRead;
        }

        std::cerr << "FAILED: event index " << options.eventIndex
                  << " was not found in HepMC3 input file: " << options.inputPath << '\n';
        return EXIT_FAILURE;
    } catch (const std::exception& error) {
        std::cerr << "FAILED: " << error.what() << '\n';
        printUsage(std::cerr);
        return EXIT_FAILURE;
    }
}
