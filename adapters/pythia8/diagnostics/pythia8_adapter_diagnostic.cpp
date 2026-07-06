#include "lowmc/pythia8/Pythia8Adapter.h"

#include <Pythia8/Event.h>
#include <Pythia8/Pythia.h>

#include <cstdlib>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>

#ifndef LOWMC_PYTHIA8_XMLDOC
#define LOWMC_PYTHIA8_XMLDOC "../share/Pythia8/xmldoc"
#endif

namespace {

struct Options {
    std::string mode = "toy";
    int events = 1;
    int seed = 12345;
    double eCM = 13000.0;
};

void printUsage(std::ostream& os) {
    os << "Usage: lowmc_pythia8_adapter_diagnostic [options]\n"
       << "\n"
       << "Options:\n"
       << "  --mode toy|softqcd   Select deterministic toy event or generated SoftQCD event\n"
       << "  --events N           Number of accepted events to inspect in softqcd mode\n"
       << "  --seed N             Pythia random seed for softqcd mode\n"
       << "  --ecm GEV            Pythia center-of-mass energy for softqcd mode\n"
       << "  --help               Show this help text\n";
}

int parsePositiveInt(const std::string& text, const std::string& name) {
    const int value = std::stoi(text);
    if (value <= 0) {
        throw std::invalid_argument(name + " must be positive");
    }
    return value;
}

double parsePositiveDouble(const std::string& text, const std::string& name) {
    const double value = std::stod(text);
    if (value <= 0.0) {
        throw std::invalid_argument(name + " must be positive");
    }
    return value;
}

Options parseOptions(int argc, char** argv) {
    Options options;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            printUsage(std::cout);
            std::exit(EXIT_SUCCESS);
        }
        if (arg == "--mode") {
            if (++i >= argc) {
                throw std::invalid_argument("--mode requires a value");
            }
            options.mode = argv[i];
        } else if (arg == "--events") {
            if (++i >= argc) {
                throw std::invalid_argument("--events requires a value");
            }
            options.events = parsePositiveInt(argv[i], "--events");
        } else if (arg == "--seed") {
            if (++i >= argc) {
                throw std::invalid_argument("--seed requires a value");
            }
            options.seed = parsePositiveInt(argv[i], "--seed");
        } else if (arg == "--ecm") {
            if (++i >= argc) {
                throw std::invalid_argument("--ecm requires a value");
            }
            options.eCM = parsePositiveDouble(argv[i], "--ecm");
        } else {
            throw std::invalid_argument("unknown argument: " + arg);
        }
    }

    if (options.mode != "toy" && options.mode != "softqcd") {
        throw std::invalid_argument("--mode must be toy or softqcd");
    }
    return options;
}

Pythia8::Event makeToyEvent() {
    Pythia8::Event event;
    event.reset();
    event.append(2212, -12, 0, 0, 0, 0, 0, 0, 0.0, 0.0, 6500.0, 6500.0, 0.938);
    event.append(2212, -12, 0, 0, 0, 0, 0, 0, 0.0, 0.0, -6500.0, 6500.0, 0.938);
    event.append(2212, 1, 0, 0, 0, 0, 0, 0, 1.0, 0.0, 0.0, 1.37, 0.938);
    event.append(2212, 1, 0, 0, 0, 0, 0, 0, -1.0, 0.0, 0.0, 1.37, 0.938);
    event.append(22, 1, 0, 0, 0, 0, 0, 0, 0.0, 0.5, 0.0, 0.5, 0.0);
    event.append(321, -91, 0, 0, 0, 0, 0, 0, 0.0, -0.5, 0.0, 0.6, 0.494);
    return event;
}

void printDiagnosticLine(
    int eventIndex,
    const std::string& mode,
    const lowmc::AdapterParticleSelectionResult& result) {
    const lowmc::AdapterSelectionDiagnostic& diagnostic = result.diagnostic;
    std::cout << "event_index=" << eventIndex
              << " mode=" << mode
              << " selected_incoming_charged=" << diagnostic.selectedIncomingCharged
              << " selected_outgoing_charged=" << diagnostic.selectedOutgoingCharged
              << " ignored_neutral_candidates=" << diagnostic.ignoredNeutralCandidates
              << " rejected_charged_candidates=" << diagnostic.rejectedChargedCount()
              << " signed_charge_sum=" << diagnostic.signedChargeSum
              << " ordinary_charge_sum=" << diagnostic.ordinaryChargeSum
              << '\n';
}

void initializeSoftQCDPythia(Pythia8::Pythia& pythia, const Options& options) {
    pythia.readString("Print:quiet = on");
    pythia.readString("Beams:idA = 2212");
    pythia.readString("Beams:idB = 2212");
    pythia.readString("Beams:eCM = " + std::to_string(options.eCM));
    pythia.readString("SoftQCD:nonDiffractive = on");
    pythia.readString("Random:setSeed = on");
    pythia.readString("Random:seed = " + std::to_string(options.seed));
    if (!pythia.init()) {
        throw std::runtime_error("Pythia initialization failed");
    }
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Options options = parseOptions(argc, argv);

        if (options.mode == "toy") {
            const Pythia8::Event event = makeToyEvent();
            printDiagnosticLine(0, options.mode, lowmc::pythia8::selectChargedParticles(event));
            return EXIT_SUCCESS;
        }

        Pythia8::Pythia pythia(LOWMC_PYTHIA8_XMLDOC, false);
        initializeSoftQCDPythia(pythia, options);
        int acceptedEvents = 0;
        int attemptedEvents = 0;
        while (acceptedEvents < options.events) {
            ++attemptedEvents;
            if (!pythia.next()) {
                if (attemptedEvents > options.events * 20) {
                    throw std::runtime_error("too many Pythia event-generation failures");
                }
                continue;
            }

            printDiagnosticLine(
                acceptedEvents,
                options.mode,
                lowmc::pythia8::selectChargedParticles(pythia.event));
            ++acceptedEvents;
        }
    } catch (const std::exception& error) {
        std::cerr << "FAILED: " << error.what() << '\n';
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
