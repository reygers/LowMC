#include "lowmc/LowMC.h"
#include "lowmc/pythia8/Pythia8Adapter.h"

#include <Pythia8/Pythia.h>

#include <algorithm>
#include <chrono>
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

#ifndef LOWMC_PYTHIA8_XMLDOC
#define LOWMC_PYTHIA8_XMLDOC "../share/Pythia8/xmldoc"
#endif

namespace {

using Clock = std::chrono::steady_clock;

struct Options {
    std::string outputPath = "pythia8_leg_proposal_tail.csv";
    std::string mode = "ee-hadrons";
    int events = 10;
    int firstEvent = 0;
    int samples = 10000;
    std::uint64_t seed = 12345;
    std::uint64_t samplingSeed = 0;
    double eCM = 0.0;
    double kT = 0.02;
    double etaMin = -5.0;
    double etaMax = 5.0;
    double phiMin = 0.0;
    double phiMax = 2.0 * std::acos(-1.0);
    double legUniformFraction = 0.1;
};

struct TailSummary {
    std::size_t count = 0;
    double mean = 0.0;
    double standardError = 0.0;
    double effectiveSampleFraction = 0.0;
    double maximum = 0.0;
    double q99 = 0.0;
    double q999 = 0.0;
    double q9999 = 0.0;
    double insideFraction = 0.0;
};

struct ProposalRun {
    TailSummary summary;
    double buildMilliseconds = 0.0;
    double sampleMilliseconds = 0.0;
    std::size_t proposalElements = 0;
};

void printUsage(std::ostream& output) {
    output
        << "Usage: lowmc_pythia8_leg_proposal_tail_diagnostic [output.csv] [options]\n\n"
        << "Compare the diagnostic analytic leg-centred and tuned-grid proposals\n"
        << "on generated Pythia events at one fixed kT. No envelope is built.\n\n"
        << "Options:\n"
        << "  --mode ee-hadrons|softqcd\n"
        << "  --events N\n"
        << "  --first-event N         Skip this many valid selected events\n"
        << "  --samples N\n"
        << "  --seed N\n"
        << "  --sampling-seed N       Independent proposal RNG seed\n"
        << "  --ecm GEV              Defaults to mZ or 13000 according to mode\n"
        << "  --kt GEV\n"
        << "  --eta-min X --eta-max X --phi-min X --phi-max X\n"
        << "  --leg-uniform-fraction F\n"
        << "  --help\n";
}

std::string readStringOption(
    int& index, int argc, char** argv, const std::string& optionName) {
    if (index + 1 >= argc) {
        throw std::invalid_argument(optionName + " requires a value");
    }
    return argv[++index];
}

int parsePositiveInt(const std::string& text, const std::string& optionName) {
    std::size_t parsed = 0;
    const long value = std::stol(text, &parsed);
    if (parsed != text.size() || value <= 0 ||
        value > std::numeric_limits<int>::max()) {
        throw std::invalid_argument(optionName + " expects a positive integer");
    }
    return static_cast<int>(value);
}

int parseNonNegativeInt(const std::string& text, const std::string& optionName) {
    std::size_t parsed = 0;
    const long value = std::stol(text, &parsed);
    if (parsed != text.size() || value < 0 ||
        value > std::numeric_limits<int>::max()) {
        throw std::invalid_argument(optionName + " expects a non-negative integer");
    }
    return static_cast<int>(value);
}

std::uint64_t parseUInt64(const std::string& text, const std::string& optionName) {
    std::size_t parsed = 0;
    const unsigned long long value = std::stoull(text, &parsed);
    if (parsed != text.size()) {
        throw std::invalid_argument(optionName + " expects an unsigned integer");
    }
    return static_cast<std::uint64_t>(value);
}

double parseFiniteDouble(const std::string& text, const std::string& optionName) {
    std::size_t parsed = 0;
    const double value = std::stod(text, &parsed);
    if (parsed != text.size() || !std::isfinite(value)) {
        throw std::invalid_argument(optionName + " expects a finite number");
    }
    return value;
}

double parsePositiveDouble(const std::string& text, const std::string& optionName) {
    const double value = parseFiniteDouble(text, optionName);
    if (value <= 0.0) {
        throw std::invalid_argument(optionName + " expects a positive number");
    }
    return value;
}

Options parseOptions(int argc, char** argv) {
    Options options;
    bool outputSpecified = false;
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument == "--help" || argument == "-h") {
            printUsage(std::cout);
            std::exit(EXIT_SUCCESS);
        }
        if (argument == "--mode") {
            options.mode = readStringOption(index, argc, argv, argument);
        } else if (argument == "--events") {
            options.events = parsePositiveInt(
                readStringOption(index, argc, argv, argument), argument);
        } else if (argument == "--first-event") {
            options.firstEvent = parseNonNegativeInt(
                readStringOption(index, argc, argv, argument), argument);
        } else if (argument == "--samples") {
            options.samples = parsePositiveInt(
                readStringOption(index, argc, argv, argument), argument);
        } else if (argument == "--seed") {
            options.seed = parseUInt64(
                readStringOption(index, argc, argv, argument), argument);
        } else if (argument == "--sampling-seed") {
            options.samplingSeed = parseUInt64(
                readStringOption(index, argc, argv, argument), argument);
        } else if (argument == "--ecm") {
            options.eCM = parsePositiveDouble(
                readStringOption(index, argc, argv, argument), argument);
        } else if (argument == "--kt") {
            options.kT = parsePositiveDouble(
                readStringOption(index, argc, argv, argument), argument);
        } else if (argument == "--eta-min") {
            options.etaMin = parseFiniteDouble(
                readStringOption(index, argc, argv, argument), argument);
        } else if (argument == "--eta-max") {
            options.etaMax = parseFiniteDouble(
                readStringOption(index, argc, argv, argument), argument);
        } else if (argument == "--phi-min") {
            options.phiMin = parseFiniteDouble(
                readStringOption(index, argc, argv, argument), argument);
        } else if (argument == "--phi-max") {
            options.phiMax = parseFiniteDouble(
                readStringOption(index, argc, argv, argument), argument);
        } else if (argument == "--leg-uniform-fraction") {
            options.legUniformFraction = parseFiniteDouble(
                readStringOption(index, argc, argv, argument), argument);
        } else if (!argument.empty() && argument.front() == '-') {
            throw std::invalid_argument("unknown argument: " + argument);
        } else if (outputSpecified) {
            throw std::invalid_argument("only one output CSV path can be specified");
        } else {
            options.outputPath = argument;
            outputSpecified = true;
        }
    }

    if (options.mode != "ee-hadrons" && options.mode != "softqcd") {
        throw std::invalid_argument("--mode must be ee-hadrons or softqcd");
    }
    if (!(options.legUniformFraction > 0.0 && options.legUniformFraction < 1.0)) {
        throw std::invalid_argument(
            "--leg-uniform-fraction must satisfy 0 < f < 1");
    }
    const lowmc::PhotonPhaseSpace phaseSpace{
        options.kT,
        std::nextafter(options.kT, std::numeric_limits<double>::infinity()),
        options.etaMin,
        options.etaMax,
        options.phiMin,
        options.phiMax,
    };
    if (!lowmc::validatePhotonPhaseSpace(phaseSpace).valid()) {
        throw std::invalid_argument("invalid fixed-kT photon phase space");
    }
    return options;
}

void initializePythia(Pythia8::Pythia& pythia, const Options& options) {
    pythia.readString("Print:quiet = on");
    pythia.readString("Random:setSeed = on");
    pythia.readString("Random:seed = " + std::to_string(options.seed));
    if (options.mode == "softqcd") {
        pythia.readString("Beams:idA = 2212");
        pythia.readString("Beams:idB = 2212");
        pythia.readString(
            "Beams:eCM = " + std::to_string(options.eCM > 0.0 ? options.eCM : 13000.0));
        pythia.readString("SoftQCD:nonDiffractive = on");
    } else {
        pythia.readString("PDF:lepton = off");
        pythia.readString("WeakSingleBoson:ffbar2gmZ = on");
        pythia.readString("23:onMode = off");
        pythia.readString("23:onIfAny = 1 2 3 4 5");
        pythia.readString("Beams:idA = 11");
        pythia.readString("Beams:idB = -11");
        const double energy = options.eCM > 0.0 ? options.eCM : pythia.particleData.m0(23);
        pythia.readString("Beams:eCM = " + std::to_string(energy));
    }
    if (!pythia.init()) {
        throw std::runtime_error("Pythia initialization failed");
    }
}

double percentile(const std::vector<double>& sorted, double probability) {
    if (sorted.empty()) {
        return 0.0;
    }
    const double position = probability * static_cast<double>(sorted.size() - 1);
    const std::size_t lower = static_cast<std::size_t>(position);
    const std::size_t upper = std::min(lower + 1, sorted.size() - 1);
    const double fraction = position - static_cast<double>(lower);
    return sorted[lower] + fraction * (sorted[upper] - sorted[lower]);
}

TailSummary summarize(std::vector<double> weights, std::size_t insideCount) {
    TailSummary summary;
    summary.count = weights.size();
    if (weights.empty()) {
        return summary;
    }
    double sum = 0.0;
    double sum2 = 0.0;
    for (double weight : weights) {
        sum += weight;
        sum2 += weight * weight;
    }
    summary.mean = sum / static_cast<double>(weights.size());
    const double variance = std::max(
        0.0, sum2 / static_cast<double>(weights.size()) - summary.mean * summary.mean);
    summary.standardError = std::sqrt(variance / static_cast<double>(weights.size()));
    summary.effectiveSampleFraction =
        sum2 > 0.0 ? sum * sum / (static_cast<double>(weights.size()) * sum2) : 0.0;
    std::sort(weights.begin(), weights.end());
    summary.maximum = weights.back();
    summary.q99 = percentile(weights, 0.99);
    summary.q999 = percentile(weights, 0.999);
    summary.q9999 = percentile(weights, 0.9999);
    summary.insideFraction =
        static_cast<double>(insideCount) / static_cast<double>(weights.size());
    return summary;
}

ProposalRun runLegProposal(
    const lowmc::LowCurrent& current,
    const lowmc::PhotonPhaseSpace& phaseSpace,
    double kT,
    int sampleCount,
    std::uint64_t seed,
    double uniformFraction) {
    const auto buildStart = Clock::now();
    const lowmc::LegCenteredEtaPhiProposal proposal =
        lowmc::makeLegCenteredEtaPhiProposal(
            phaseSpace,
            current.particles(),
            lowmc::LegCenteredEtaPhiProposalConfig{uniformFraction, 1.0e-12});
    const auto buildStop = Clock::now();

    std::mt19937_64 rng(seed);
    std::vector<double> weights;
    weights.reserve(static_cast<std::size_t>(sampleCount));
    std::size_t insideCount = 0;
    const auto sampleStart = Clock::now();
    for (int index = 0; index < sampleCount; ++index) {
        const lowmc::FixedKTLegCenteredWeightedSample sample =
            lowmc::makeFixedKTLegCenteredEtaPhiWeightedSample(
                current,
                kT,
                proposal,
                std::generate_canonical<double, 53>(rng),
                std::generate_canonical<double, 53>(rng),
                std::generate_canonical<double, 53>(rng));
        weights.push_back(sample.weightedPhoton.photon.weight);
        insideCount += sample.etaPhiProposal.insideAcceptance ? 1 : 0;
    }
    const auto sampleStop = Clock::now();
    return ProposalRun{
        summarize(std::move(weights), insideCount),
        std::chrono::duration<double, std::milli>(buildStop - buildStart).count(),
        std::chrono::duration<double, std::milli>(sampleStop - sampleStart).count(),
        proposal.components.size(),
    };
}

ProposalRun runGridProposal(
    const lowmc::LowCurrent& current,
    const lowmc::PhotonPhaseSpace& phaseSpace,
    double kT,
    int sampleCount,
    std::uint64_t seed) {
    const auto buildStart = Clock::now();
    const lowmc::GridAdaptiveEtaPhiProposal proposal =
        lowmc::makeGridAdaptiveEtaPhiProposal(
            current,
            phaseSpace,
            kT,
            lowmc::tunedMassGridAdaptiveEtaPhiConfig());
    const auto buildStop = Clock::now();

    std::mt19937_64 rng(seed);
    std::vector<double> weights;
    weights.reserve(static_cast<std::size_t>(sampleCount));
    const auto sampleStart = Clock::now();
    for (int index = 0; index < sampleCount; ++index) {
        const lowmc::GridAdaptiveEtaPhiProposalPoint point =
            lowmc::sampleGridAdaptiveEtaPhi(
                proposal,
                std::generate_canonical<double, 53>(rng),
                std::generate_canonical<double, 53>(rng),
                std::generate_canonical<double, 53>(rng),
                std::generate_canonical<double, 53>(rng));
        const lowmc::FourVector momentum =
            lowmc::photonMomentumFromKTEtaPhi(kT, point.eta, point.phi);
        weights.push_back(
            lowmc::makeWeightedPhotonSample(
                current,
                phaseSpace,
                lowmc::WeightedPhotonProposal{momentum, point.density})
                .photon.weight);
    }
    const auto sampleStop = Clock::now();
    return ProposalRun{
        summarize(std::move(weights), static_cast<std::size_t>(sampleCount)),
        std::chrono::duration<double, std::milli>(buildStop - buildStart).count(),
        std::chrono::duration<double, std::milli>(sampleStop - sampleStart).count(),
        proposal.cells.size(),
    };
}

double particleBeta(const lowmc::ChargedParticle& particle) {
    const double momentum = std::sqrt(
        particle.momentum.px * particle.momentum.px +
        particle.momentum.py * particle.momentum.py +
        particle.momentum.pz * particle.momentum.pz);
    return particle.momentum.E != 0.0
               ? momentum / std::abs(particle.momentum.E)
               : 0.0;
}

double particleAbsEta(const lowmc::ChargedParticle& particle) {
    const double pt = particle.momentum.pt();
    if (pt <= 0.0) {
        return std::numeric_limits<double>::infinity();
    }
    return std::abs(std::asinh(particle.momentum.pz / pt));
}

void writeRow(
    std::ofstream& output,
    const Options& options,
    int eventIndex,
    std::uint64_t samplingSeed,
    const lowmc::AdapterParticleSelectionResult& selection,
    double maximumOutgoingBeta,
    double maximumOutgoingAbsEta,
    const char* proposalName,
    const ProposalRun& run) {
    const TailSummary& summary = run.summary;
    output << options.mode << ','
           << options.seed << ','
           << eventIndex << ','
           << samplingSeed << ','
           << selection.selectedParticles.size() << ','
           << selection.diagnostic.selectedIncomingCharged << ','
           << selection.diagnostic.selectedOutgoingCharged << ','
           << maximumOutgoingBeta << ','
           << maximumOutgoingAbsEta << ','
           << proposalName << ','
           << summary.count << ','
           << summary.mean << ','
           << summary.standardError << ','
           << summary.effectiveSampleFraction << ','
           << summary.maximum << ','
           << summary.q99 << ','
           << summary.q999 << ','
           << summary.q9999 << ','
           << summary.insideFraction << ','
           << run.buildMilliseconds << ','
           << run.sampleMilliseconds << ','
           << run.proposalElements << '\n';
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Options options = parseOptions(argc, argv);
        const lowmc::PhotonPhaseSpace phaseSpace{
            options.kT,
            std::nextafter(options.kT, std::numeric_limits<double>::infinity()),
            options.etaMin,
            options.etaMax,
            options.phiMin,
            options.phiMax,
        };

        std::ofstream output(options.outputPath);
        if (!output) {
            throw std::runtime_error("could not open output CSV: " + options.outputPath);
        }
        output << "mode,seed,event_index,sampling_seed,selected_charged,selected_incoming,"
                  "selected_outgoing,max_outgoing_beta,max_outgoing_abs_eta,proposal,"
                  "samples,mean_weight,"
                  "standard_error,effective_sample_fraction,max_weight,q99_weight,"
                  "q999_weight,q9999_weight,inside_fraction,build_milliseconds,"
                  "sample_milliseconds,proposal_elements\n";
        output << std::setprecision(17);

        Pythia8::Pythia pythia(LOWMC_PYTHIA8_XMLDOC, false);
        initializePythia(pythia, options);
        int validEvents = 0;
        int analyzedEvents = 0;
        int attemptedEvents = 0;
        while (analyzedEvents < options.events) {
            ++attemptedEvents;
            if (!pythia.next()) {
                if (attemptedEvents > (options.firstEvent + options.events) * 20) {
                    throw std::runtime_error("too many Pythia event-generation failures");
                }
                continue;
            }
            const lowmc::AdapterParticleSelectionResult selection =
                lowmc::pythia8::selectChargedParticles(pythia.event);
            if (std::abs(selection.diagnostic.signedChargeSum) > 1.0e-12 ||
                selection.selectedParticles.size() < 2) {
                continue;
            }
            const int eventIndex = validEvents++;
            if (eventIndex < options.firstEvent) {
                continue;
            }

            const lowmc::LowCurrent current(selection.selectedParticles);
            double maximumOutgoingBeta = 0.0;
            double maximumOutgoingAbsEta = 0.0;
            for (const lowmc::ChargedParticle& particle : selection.selectedParticles) {
                if (particle.inOutSign != lowmc::InOutSign::Outgoing) {
                    continue;
                }
                maximumOutgoingBeta =
                    std::max(maximumOutgoingBeta, particleBeta(particle));
                maximumOutgoingAbsEta =
                    std::max(maximumOutgoingAbsEta, particleAbsEta(particle));
            }

            const std::uint64_t eventSeed =
                (options.samplingSeed != 0 ? options.samplingSeed : options.seed) +
                104729ULL * static_cast<std::uint64_t>(eventIndex + 1);
            const ProposalRun leg = runLegProposal(
                current,
                phaseSpace,
                options.kT,
                options.samples,
                eventSeed,
                options.legUniformFraction);
            const ProposalRun grid = runGridProposal(
                current,
                phaseSpace,
                options.kT,
                options.samples,
                eventSeed + 1ULL);
            writeRow(
                output,
                options,
                eventIndex,
                eventSeed,
                selection,
                maximumOutgoingBeta,
                maximumOutgoingAbsEta,
                "leg",
                leg);
            writeRow(
                output,
                options,
                eventIndex,
                eventSeed + 1ULL,
                selection,
                maximumOutgoingBeta,
                maximumOutgoingAbsEta,
                "grid",
                grid);
            ++analyzedEvents;
        }
        if (!output) {
            throw std::runtime_error("failed while writing output CSV");
        }
    } catch (const std::exception& error) {
        std::cerr << "FAILED: " << error.what() << '\n';
        printUsage(std::cerr);
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
