#include "lowmc/hepmc3/HepMC3Adapter.h"
#include "lowmc/LowMC.h"
#include "lowmc/pythia8/Pythia8Adapter.h"

#include <HepMC3/GenEvent.h>
#include <Pythia8/Pythia.h>
#include <Pythia8Plugins/HepMC3.h>

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
    int events = 1;
    int samples = 8;
    int oversamplingFactor = 1;
    int seed = 12345;
    std::uint64_t photonSeed = 12345;
    double eCM = 13000.0;
    double kTMin = 0.002;
    double kTMax = 0.2;
    double etaMin = -5.0;
    double etaMax = 5.0;
    double phiMin = 0.0;
    double phiMax = 2.0 * std::acos(-1.0);
};

void printUsage(std::ostream& os) {
    os << "Usage: lowmc_pythia8_hepmc3_crosscheck [options]\n"
       << "\n"
       << "Generates deterministic Pythia8 events, runs the direct Pythia8\n"
       << "selected-particle path, converts the same event to HepMC3 in memory,\n"
       << "and compares the HepMC3 selected-particle path. This is a diagnostic\n"
       << "cross-check only; it writes no event files and generates no photons.\n"
       << "\n"
       << "Options:\n"
       << "  --events N   Number of accepted Pythia events to compare\n"
       << "  --samples N  Weighted photon samples per selected event\n"
       << "  --oversampling-factor N\n"
       << "               Multiply samples and divide each photon weight by N\n"
       << "  --seed N     Pythia random seed\n"
       << "  --photon-seed N\n"
       << "               Weighted photon random seed\n"
       << "  --ecm GEV    Pythia center-of-mass energy\n"
       << "  --kt-min GEV Minimum photon k_T\n"
       << "  --kt-max GEV Maximum photon k_T\n"
       << "  --eta-min X  Minimum photon eta\n"
       << "  --eta-max X  Maximum photon eta\n"
       << "  --phi-min X  Minimum photon phi\n"
       << "  --phi-max X  Maximum photon phi\n"
       << "  --help       Show this help text\n";
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

std::uint64_t parseUInt64(const std::string& value, const std::string& optionName) {
    std::size_t parsed = 0;
    const unsigned long long parsedValue = std::stoull(value, &parsed);
    if (parsed != value.size()) {
        throw std::invalid_argument(optionName + " expects an unsigned integer");
    }
    return static_cast<std::uint64_t>(parsedValue);
}

Options parseOptions(int argc, char** argv) {
    Options options;
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument == "--help" || argument == "-h") {
            printUsage(std::cout);
            std::exit(EXIT_SUCCESS);
        }
        if (argument == "--events") {
            options.events = parsePositiveInt(readStringOption(index, argc, argv, argument), argument);
        } else if (argument == "--samples") {
            options.samples = parsePositiveInt(readStringOption(index, argc, argv, argument), argument);
        } else if (argument == "--oversampling-factor") {
            options.oversamplingFactor =
                parsePositiveInt(readStringOption(index, argc, argv, argument), argument);
        } else if (argument == "--seed") {
            options.seed = parsePositiveInt(readStringOption(index, argc, argv, argument), argument);
        } else if (argument == "--photon-seed") {
            options.photonSeed = parseUInt64(readStringOption(index, argc, argv, argument), argument);
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

[[nodiscard]] bool nearlyEqual(double first, double second) {
    return std::abs(first - second) < 1.0e-12;
}

[[nodiscard]] bool nearlyEqualWeighted(double first, double second) {
    const double scale = std::max({1.0, std::abs(first), std::abs(second)});
    return std::abs(first - second) <= 1.0e-9 * scale;
}

[[nodiscard]] bool diagnosticsMatch(
    const lowmc::AdapterSelectionDiagnostic& direct,
    const lowmc::AdapterSelectionDiagnostic& hepmc3) {
    return direct.selectedIncomingCharged == hepmc3.selectedIncomingCharged &&
           direct.selectedOutgoingCharged == hepmc3.selectedOutgoingCharged &&
           direct.ignoredNeutralCandidates == hepmc3.ignoredNeutralCandidates &&
           direct.rejectedChargedCount() == hepmc3.rejectedChargedCount() &&
           nearlyEqual(direct.signedChargeSum, hepmc3.signedChargeSum) &&
           nearlyEqual(direct.ordinaryChargeSum, hepmc3.ordinaryChargeSum);
}

[[nodiscard]] double weightRms(const lowmc::WeightSummary& summary) {
    return std::sqrt(summary.varianceWeight);
}

[[nodiscard]] double effectiveSampleFraction(const lowmc::WeightSummary& summary) {
    return summary.count > 0
        ? summary.effectiveSampleSize / static_cast<double>(summary.count)
        : 0.0;
}

[[nodiscard]] bool weightedSummariesMatch(
    const lowmc::WeightedEventGenerationResult& direct,
    const lowmc::WeightedEventGenerationResult& hepmc3) {
    const lowmc::WeightSummary& directSummary = direct.weightSummary;
    const lowmc::WeightSummary& hepmc3Summary = hepmc3.weightSummary;
    return directSummary.count == hepmc3Summary.count &&
           nearlyEqualWeighted(directSummary.sumWeight, hepmc3Summary.sumWeight) &&
           nearlyEqualWeighted(directSummary.meanWeight, hepmc3Summary.meanWeight) &&
           nearlyEqualWeighted(weightRms(directSummary), weightRms(hepmc3Summary)) &&
           nearlyEqualWeighted(
               directSummary.effectiveSampleSize,
               hepmc3Summary.effectiveSampleSize) &&
           nearlyEqualWeighted(
               effectiveSampleFraction(directSummary),
               effectiveSampleFraction(hepmc3Summary)) &&
           nearlyEqualWeighted(direct.gridBuildKT, hepmc3.gridBuildKT);
}

[[nodiscard]] lowmc::WeightedEventGeneratorConfig makeGeneratorConfig(const Options& options) {
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
        options.photonSeed,
        lowmc::tunedMassGridAdaptiveEtaPhiConfig(),
        options.oversamplingFactor,
    };
}

void printComparisonLine(
    int eventIndex,
    const lowmc::AdapterSelectionDiagnostic& direct,
    const lowmc::AdapterSelectionDiagnostic& hepmc3,
    const lowmc::WeightedEventGenerationResult& directWeighted,
    const lowmc::WeightedEventGenerationResult& hepmc3Weighted,
    bool selectionMatch,
    bool weightedMatch) {
    const lowmc::WeightSummary& directSummary = directWeighted.weightSummary;
    const lowmc::WeightSummary& hepmc3Summary = hepmc3Weighted.weightSummary;
    std::cout << "event_index=" << eventIndex
              << " selection_status=" << (selectionMatch ? "match" : "mismatch")
              << " direct_selected_incoming_charged=" << direct.selectedIncomingCharged
              << " hepmc3_selected_incoming_charged=" << hepmc3.selectedIncomingCharged
              << " direct_selected_outgoing_charged=" << direct.selectedOutgoingCharged
              << " hepmc3_selected_outgoing_charged=" << hepmc3.selectedOutgoingCharged
              << " direct_ignored_neutral_candidates=" << direct.ignoredNeutralCandidates
              << " hepmc3_ignored_neutral_candidates=" << hepmc3.ignoredNeutralCandidates
              << " direct_rejected_charged_candidates=" << direct.rejectedChargedCount()
              << " hepmc3_rejected_charged_candidates=" << hepmc3.rejectedChargedCount()
              << " direct_signed_charge_sum=" << direct.signedChargeSum
              << " hepmc3_signed_charge_sum=" << hepmc3.signedChargeSum
              << " direct_ordinary_charge_sum=" << direct.ordinaryChargeSum
              << " hepmc3_ordinary_charge_sum=" << hepmc3.ordinaryChargeSum
              << " weighted_status=" << (weightedMatch ? "match" : "mismatch")
              << " direct_samples=" << directSummary.count
              << " hepmc3_samples=" << hepmc3Summary.count
              << " direct_sum_weight=" << directSummary.sumWeight
              << " hepmc3_sum_weight=" << hepmc3Summary.sumWeight
              << " direct_mean_weight=" << directSummary.meanWeight
              << " hepmc3_mean_weight=" << hepmc3Summary.meanWeight
              << " direct_weight_rms=" << weightRms(directSummary)
              << " hepmc3_weight_rms=" << weightRms(hepmc3Summary)
              << " direct_effective_sample_size=" << directSummary.effectiveSampleSize
              << " hepmc3_effective_sample_size=" << hepmc3Summary.effectiveSampleSize
              << " direct_effective_sample_fraction=" << effectiveSampleFraction(directSummary)
              << " hepmc3_effective_sample_fraction=" << effectiveSampleFraction(hepmc3Summary)
              << " direct_grid_build_kt=" << directWeighted.gridBuildKT
              << " hepmc3_grid_build_kt=" << hepmc3Weighted.gridBuildKT
              << '\n';
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Options options = parseOptions(argc, argv);
        const lowmc::WeightedEventGeneratorConfig generatorConfig = makeGeneratorConfig(options);
        Pythia8::Pythia pythia(LOWMC_PYTHIA8_XMLDOC, false);
        initializeSoftQCDPythia(pythia, options);
        HepMC3::Pythia8ToHepMC3 converter;

        int acceptedEvents = 0;
        int attemptedEvents = 0;
        bool allMatched = true;
        while (acceptedEvents < options.events) {
            ++attemptedEvents;
            if (!pythia.next()) {
                if (attemptedEvents > options.events * 20) {
                    throw std::runtime_error("too many Pythia event-generation failures");
                }
                continue;
            }

            const lowmc::AdapterParticleSelectionResult directSelection =
                lowmc::pythia8::selectChargedParticles(pythia.event);

            HepMC3::GenEvent hepmc3Event;
            if (!converter.fill_next_event(pythia, &hepmc3Event, acceptedEvents)) {
                throw std::runtime_error("Pythia8ToHepMC3 conversion failed");
            }

            const lowmc::AdapterParticleSelectionResult hepmc3Selection =
                lowmc::hepmc3::selectChargedParticles(hepmc3Event);

            const bool selectionMatch =
                diagnosticsMatch(directSelection.diagnostic, hepmc3Selection.diagnostic);
            if (!selectionMatch) {
                lowmc::WeightedEventGenerationResult emptyDirect;
                lowmc::WeightedEventGenerationResult emptyHepMC3;
                printComparisonLine(
                    acceptedEvents,
                    directSelection.diagnostic,
                    hepmc3Selection.diagnostic,
                    emptyDirect,
                    emptyHepMC3,
                    false,
                    false);
                allMatched = false;
                ++acceptedEvents;
                continue;
            }

            const lowmc::WeightedEventGenerationResult directWeighted =
                lowmc::generateWeightedPhotons(
                    directSelection.selectedParticles,
                    generatorConfig);
            const lowmc::WeightedEventGenerationResult hepmc3Weighted =
                lowmc::generateWeightedPhotons(
                    hepmc3Selection.selectedParticles,
                    generatorConfig);
            const bool weightedMatch = weightedSummariesMatch(directWeighted, hepmc3Weighted);
            printComparisonLine(
                acceptedEvents,
                directSelection.diagnostic,
                hepmc3Selection.diagnostic,
                directWeighted,
                hepmc3Weighted,
                selectionMatch,
                weightedMatch);
            allMatched = allMatched && selectionMatch && weightedMatch;
            ++acceptedEvents;
        }

        return allMatched ? EXIT_SUCCESS : EXIT_FAILURE;
    } catch (const std::exception& error) {
        std::cerr << "FAILED: " << error.what() << '\n';
        printUsage(std::cerr);
        return EXIT_FAILURE;
    }
}
