#include "lowmc/LowMC.h"
#include "lowmc/hepmc3/HepMC3Adapter.h"

#include <HepMC3/GenEvent.h>
#include <HepMC3/ReaderAscii.h>
#include <HepMC3/Units.h>

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {

struct Options {
    std::string inputPath;
    int maxEvents = 1;
    int samples = 16;
    int oversamplingFactor = 1;
    std::uint64_t seed = 12345;
    double kTMin = 0.002;
    double kTMax = 0.2;
    double etaMin = -5.0;
    double etaMax = 5.0;
    double phiMin = 0.0;
    double phiMax = 2.0 * std::acos(-1.0);
};

void printUsage(std::ostream& out) {
    out << "Usage: lowmc_hepmc3_weighted_photon_diagnostic [options] <input.hepmc3>\n"
        << "\n"
        << "Reads HepMC3 ASCII events, selects charged particles with the LowMC HepMC3\n"
        << "adapter, and prints compact weighted-photon diagnostics. No photons are\n"
        << "written back into the HepMC3 event record.\n"
        << "\n"
        << "Options:\n"
        << "  --max-events N   Number of input events to inspect (default: 1)\n"
        << "  --samples N      Base weighted photon samples per accepted event (default: 16)\n"
        << "  --oversampling-factor N\n"
        << "                  Multiply samples and divide each photon weight by N (default: 1)\n"
        << "  --seed N         Random seed (default: 12345)\n"
        << "  --kt-min GeV     Minimum photon k_T (default: 0.002)\n"
        << "  --kt-max GeV     Maximum photon k_T (default: 0.2)\n"
        << "  --eta-min X      Minimum photon eta (default: -5)\n"
        << "  --eta-max X      Maximum photon eta (default: 5)\n"
        << "  --phi-min X      Minimum photon phi (default: 0)\n"
        << "  --phi-max X      Maximum photon phi (default: 2 pi)\n";
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

double parseDouble(const std::string& value, const std::string& optionName) {
    std::size_t parsed = 0;
    const double parsedValue = std::stod(value, &parsed);
    if (parsed != value.size() || !std::isfinite(parsedValue)) {
        throw std::invalid_argument(optionName + " expects a finite number");
    }
    return parsedValue;
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
        if (argument == "--max-events") {
            options.maxEvents =
                parsePositiveInt(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--samples") {
            options.samples =
                parsePositiveInt(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--oversampling-factor") {
            options.oversamplingFactor =
                parsePositiveInt(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--seed") {
            options.seed = parseUInt64(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--kt-min") {
            options.kTMin = parseDouble(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--kt-max") {
            options.kTMax = parseDouble(readStringOption(index, argc, argv, argument), argument);
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

        if (!argument.empty() && argument.front() == '-') {
            throw std::invalid_argument("unknown option: " + argument);
        }
        if (!options.inputPath.empty()) {
            throw std::invalid_argument("expected exactly one input file path");
        }
        options.inputPath = argument;
    }

    if (options.inputPath.empty()) {
        throw std::invalid_argument("missing input file path");
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

void printAdapterDiagnosticPrefix(
    int eventIndex,
    HepMC3::Units::MomentumUnit momentumUnit,
    const lowmc::AdapterSelectionDiagnostic& diagnostic) {
    std::cout << "event_index=" << eventIndex
              << " momentum_unit=" << HepMC3::Units::name(momentumUnit)
              << " selected_incoming_charged=" << diagnostic.selectedIncomingCharged
              << " selected_outgoing_charged=" << diagnostic.selectedOutgoingCharged
              << " ignored_neutral_candidates=" << diagnostic.ignoredNeutralCandidates
              << " rejected_charged_candidates=" << diagnostic.rejectedChargedCount()
              << " signed_charge_sum=" << diagnostic.signedChargeSum
              << " ordinary_charge_sum=" << diagnostic.ordinaryChargeSum;
}

void printSkippedPhotonDiagnostic(
    int eventIndex,
    HepMC3::Units::MomentumUnit momentumUnit,
    const lowmc::AdapterParticleSelectionResult& selection) {
    printAdapterDiagnosticPrefix(eventIndex, momentumUnit, selection.diagnostic);
    std::cout << " photon_status=skipped_nonconserved_charge\n";
}

void printGeneratedPhotonDiagnostic(
    int eventIndex,
    HepMC3::Units::MomentumUnit momentumUnit,
    const lowmc::AdapterParticleSelectionResult& selection,
    const lowmc::WeightedEventGenerationResult& result,
    int baseSamples,
    int oversamplingFactor) {
    const lowmc::WeightSummary& summary = result.weightSummary;
    const double effectiveSampleFraction =
        summary.count > 0
            ? summary.effectiveSampleSize / static_cast<double>(summary.count)
            : 0.0;

    printAdapterDiagnosticPrefix(eventIndex, momentumUnit, selection.diagnostic);
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

        const lowmc::WeightedEventGeneratorConfig config{
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

        int eventsRead = 0;
        while (eventsRead < options.maxEvents && !reader.failed()) {
            HepMC3::GenEvent event;
            if (!reader.read_event(event)) {
                break;
            }

            const lowmc::AdapterParticleSelectionResult selection =
                lowmc::hepmc3::selectChargedParticles(event);

            if (std::abs(selection.diagnostic.signedChargeSum) > 1.0e-12) {
                printSkippedPhotonDiagnostic(eventsRead, event.momentum_unit(), selection);
            } else {
                const lowmc::WeightedEventGenerationResult result =
                    lowmc::generateWeightedPhotons(selection.selectedParticles, config);
                printGeneratedPhotonDiagnostic(
                    eventsRead,
                    event.momentum_unit(),
                    selection,
                    result,
                    options.samples,
                    options.oversamplingFactor);
            }

            ++eventsRead;
        }

        if (eventsRead == 0) {
            std::cerr << "FAILED: no events read from HepMC3 input file: "
                      << options.inputPath << '\n';
            return EXIT_FAILURE;
        }
    } catch (const std::exception& error) {
        std::cerr << "FAILED: " << error.what() << '\n';
        printUsage(std::cerr);
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
