#include "lowmc/LowMC.h"
#include "lowmc/pythia8/Pythia8Adapter.h"

#include <Pythia8/Event.h>
#include <Pythia8/Pythia.h>

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

#ifndef LOWMC_PYTHIA8_XMLDOC
#define LOWMC_PYTHIA8_XMLDOC "../share/Pythia8/xmldoc"
#endif

namespace {

struct Options {
    std::string mode = "toy";
    int events = 1;
    int samples = 16;
    int oversamplingFactor = 1;
    std::uint64_t seed = 12345;
    double eCM = 13000.0;
    double kTMin = 0.002;
    double kTMax = 0.2;
    double etaMin = -5.0;
    double etaMax = 5.0;
    double phiMin = 0.0;
    double phiMax = 2.0 * std::acos(-1.0);
};

void printUsage(std::ostream& os) {
    os << "Usage: lowmc_pythia8_weighted_photon_diagnostic [options]\n"
       << "\n"
       << "Runs the direct Pythia8 in-memory adapter, selects charged particles,\n"
       << "and prints compact weighted Low-photon diagnostics. No HepMC3 file is\n"
       << "read or written.\n"
       << "\n"
       << "Options:\n"
       << "  --mode toy|toy-nonconserved|softqcd\n"
       << "                  Select deterministic toy input or generated SoftQCD events\n"
       << "  --events N      Number of accepted events to inspect in softqcd mode\n"
       << "  --samples N     Base weighted photon samples per accepted event\n"
       << "  --oversampling-factor N\n"
       << "                  Multiply samples and divide each photon weight by N\n"
       << "  --seed N        Random seed for Pythia and weighted photons\n"
       << "  --ecm GEV       Pythia center-of-mass energy in softqcd mode\n"
       << "  --kt-min GEV    Minimum photon k_T\n"
       << "  --kt-max GEV    Maximum photon k_T\n"
       << "  --eta-min X     Minimum photon eta\n"
       << "  --eta-max X     Maximum photon eta\n"
       << "  --phi-min X     Minimum photon phi\n"
       << "  --phi-max X     Maximum photon phi\n"
       << "  --help          Show this help text\n";
}

std::string readStringOption(int& index, int argc, char** argv, const std::string& optionName) {
    if (index + 1 >= argc) {
        throw std::invalid_argument(optionName + " requires a value");
    }
    return argv[++index];
}

int parsePositiveInt(const std::string& value, const std::string& optionName) {
    std::size_t parsed = 0;
    const long parsedValue = std::stol(value, &parsed);
    if (parsed != value.size() ||
        parsedValue <= 0 ||
        parsedValue > std::numeric_limits<int>::max()) {
        throw std::invalid_argument(optionName + " expects a positive integer");
    }
    return static_cast<int>(parsedValue);
}

std::uint64_t parseUInt64(const std::string& value, const std::string& optionName) {
    std::size_t parsed = 0;
    const unsigned long long parsedValue = std::stoull(value, &parsed);
    if (parsed != value.size()) {
        throw std::invalid_argument(optionName + " expects an unsigned integer");
    }
    return static_cast<std::uint64_t>(parsedValue);
}

double parsePositiveDouble(const std::string& value, const std::string& optionName) {
    std::size_t parsed = 0;
    const double parsedValue = std::stod(value, &parsed);
    if (parsed != value.size() || !std::isfinite(parsedValue) || parsedValue <= 0.0) {
        throw std::invalid_argument(optionName + " expects a positive finite number");
    }
    return parsedValue;
}

double parseFiniteDouble(const std::string& value, const std::string& optionName) {
    std::size_t parsed = 0;
    const double parsedValue = std::stod(value, &parsed);
    if (parsed != value.size() || !std::isfinite(parsedValue)) {
        throw std::invalid_argument(optionName + " expects a finite number");
    }
    return parsedValue;
}

Options parseOptions(int argc, char** argv) {
    Options options;
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument == "--help" || argument == "-h") {
            printUsage(std::cout);
            std::exit(EXIT_SUCCESS);
        }
        if (argument == "--mode") {
            options.mode = readStringOption(index, argc, argv, argument);
        } else if (argument == "--events") {
            options.events = parsePositiveInt(readStringOption(index, argc, argv, argument), argument);
        } else if (argument == "--samples") {
            options.samples = parsePositiveInt(readStringOption(index, argc, argv, argument), argument);
        } else if (argument == "--oversampling-factor") {
            options.oversamplingFactor =
                parsePositiveInt(readStringOption(index, argc, argv, argument), argument);
        } else if (argument == "--seed") {
            options.seed = parseUInt64(readStringOption(index, argc, argv, argument), argument);
        } else if (argument == "--ecm") {
            options.eCM = parsePositiveDouble(readStringOption(index, argc, argv, argument), argument);
        } else if (argument == "--kt-min") {
            options.kTMin = parsePositiveDouble(readStringOption(index, argc, argv, argument), argument);
        } else if (argument == "--kt-max") {
            options.kTMax = parsePositiveDouble(readStringOption(index, argc, argv, argument), argument);
        } else if (argument == "--eta-min") {
            options.etaMin = parseFiniteDouble(readStringOption(index, argc, argv, argument), argument);
        } else if (argument == "--eta-max") {
            options.etaMax = parseFiniteDouble(readStringOption(index, argc, argv, argument), argument);
        } else if (argument == "--phi-min") {
            options.phiMin = parseFiniteDouble(readStringOption(index, argc, argv, argument), argument);
        } else if (argument == "--phi-max") {
            options.phiMax = parseFiniteDouble(readStringOption(index, argc, argv, argument), argument);
        } else {
            throw std::invalid_argument("unknown argument: " + argument);
        }
    }

    if (options.mode != "toy" &&
        options.mode != "toy-nonconserved" &&
        options.mode != "softqcd") {
        throw std::invalid_argument("--mode must be toy, toy-nonconserved, or softqcd");
    }

    const lowmc::PhotonPhaseSpace phaseSpace{
        options.kTMin,
        options.kTMax,
        options.etaMin,
        options.etaMax,
        options.phiMin,
        options.phiMax,
    };
    if (!lowmc::validatePhotonPhaseSpace(phaseSpace).valid()) {
        throw std::invalid_argument("invalid photon phase space");
    }

    return options;
}

Pythia8::Event makeToyEvent(bool chargeConserving) {
    Pythia8::Event event;
    event.reset();
    event.append(2212, -12, 0, 0, 0, 0, 0, 0, 0.0, 0.0, 6500.0, 6500.0, 0.938);
    event.append(2212, -12, 0, 0, 0, 0, 0, 0, 0.0, 0.0, -6500.0, 6500.0, 0.938);
    if (chargeConserving) {
        event.append(2212, 1, 0, 0, 0, 0, 0, 0, 1.0, 0.0, 0.0, 1.37, 0.938);
        event.append(2212, 1, 0, 0, 0, 0, 0, 0, -1.0, 0.0, 0.0, 1.37, 0.938);
    } else {
        event.append(211, 1, 0, 0, 0, 0, 0, 0, 1.0, 0.0, 0.0, 1.01, 0.139);
        event.append(-211, 1, 0, 0, 0, 0, 0, 0, -1.0, 0.0, 0.0, 1.01, 0.139);
    }
    event.append(22, 1, 0, 0, 0, 0, 0, 0, 0.0, 0.5, 0.0, 0.5, 0.0);
    event.append(321, -91, 0, 0, 0, 0, 0, 0, 0.0, -0.5, 0.0, 0.6, 0.494);
    return event;
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

lowmc::WeightedEventGeneratorConfig makeGeneratorConfig(const Options& options) {
    return lowmc::WeightedEventGeneratorConfig{
        lowmc::PhotonPhaseSpace{
            options.kTMin,
            options.kTMax,
            options.etaMin,
            options.etaMax,
            options.phiMin,
            options.phiMax,
        },
        options.samples,
        options.seed,
        lowmc::tunedMassGridAdaptiveEtaPhiConfig(),
        options.oversamplingFactor,
    };
}

void printAdapterDiagnosticPrefix(
    int eventIndex,
    const std::string& mode,
    const lowmc::AdapterSelectionDiagnostic& diagnostic) {
    std::cout << "event_index=" << eventIndex
              << " mode=" << mode
              << " selected_incoming_charged=" << diagnostic.selectedIncomingCharged
              << " selected_outgoing_charged=" << diagnostic.selectedOutgoingCharged
              << " ignored_neutral_candidates=" << diagnostic.ignoredNeutralCandidates
              << " rejected_charged_candidates=" << diagnostic.rejectedChargedCount()
              << " signed_charge_sum=" << diagnostic.signedChargeSum
              << " ordinary_charge_sum=" << diagnostic.ordinaryChargeSum;
}

void printSkippedPhotonDiagnostic(
    int eventIndex,
    const std::string& mode,
    const lowmc::AdapterParticleSelectionResult& selection) {
    printAdapterDiagnosticPrefix(eventIndex, mode, selection.diagnostic);
    std::cout << " photon_status=skipped_nonconserved_charge\n";
}

void printGeneratedPhotonDiagnostic(
    int eventIndex,
    const std::string& mode,
    const lowmc::AdapterParticleSelectionResult& selection,
    const lowmc::WeightedEventGenerationResult& result,
    int baseSamples,
    int oversamplingFactor) {
    const lowmc::WeightSummary& summary = result.weightSummary;
    const double effectiveSampleFraction =
        summary.count > 0
            ? summary.effectiveSampleSize / static_cast<double>(summary.count)
            : 0.0;

    printAdapterDiagnosticPrefix(eventIndex, mode, selection.diagnostic);
    std::cout << " photon_status=generated"
              << " samples=" << summary.count
              << " base_samples=" << baseSamples
              << " oversampling_factor=" << oversamplingFactor
              << " sum_weight=" << summary.sumWeight
              << " mean_weight=" << summary.meanWeight
              << " weight_rms=" << std::sqrt(summary.varianceWeight)
              << " effective_sample_size=" << summary.effectiveSampleSize
              << " effective_sample_fraction=" << effectiveSampleFraction
              << " grid_build_kt=" << result.gridBuildKT
              << '\n';
}

void runSelection(
    int eventIndex,
    const std::string& mode,
    const Pythia8::Event& event,
    const lowmc::WeightedEventGeneratorConfig& config,
    int baseSamples,
    int oversamplingFactor) {
    const lowmc::AdapterParticleSelectionResult selection =
        lowmc::pythia8::selectChargedParticles(event);

    if (std::abs(selection.diagnostic.signedChargeSum) > 1.0e-12) {
        printSkippedPhotonDiagnostic(eventIndex, mode, selection);
        return;
    }

    const lowmc::WeightedEventGenerationResult result =
        lowmc::generateWeightedPhotons(selection.selectedParticles, config);
    printGeneratedPhotonDiagnostic(
        eventIndex,
        mode,
        selection,
        result,
        baseSamples,
        oversamplingFactor);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Options options = parseOptions(argc, argv);
        const lowmc::WeightedEventGeneratorConfig config = makeGeneratorConfig(options);

        if (options.mode == "toy" || options.mode == "toy-nonconserved") {
            const bool chargeConserving = options.mode == "toy";
            const Pythia8::Event event = makeToyEvent(chargeConserving);
            runSelection(
                0,
                options.mode,
                event,
                config,
                options.samples,
                options.oversamplingFactor);
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

            runSelection(
                acceptedEvents,
                options.mode,
                pythia.event,
                config,
                options.samples,
                options.oversamplingFactor);
            ++acceptedEvents;
        }
    } catch (const std::exception& error) {
        std::cerr << "FAILED: " << error.what() << '\n';
        printUsage(std::cerr);
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

