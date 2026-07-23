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
    std::string outputPath = "pythia8_leg_envelope.csv";
    std::string mode = "ee-hadrons";
    int events = 10;
    int firstEvent = 0;
    int calibrationSamples = 20000;
    int trialSamples = 100000;
    std::uint64_t seed = 12345;
    std::uint64_t samplingSeed = 0;
    double eCM = 0.0;
    double kT = 0.02;
    double etaMin = -5.0;
    double etaMax = 5.0;
    double phiMin = 0.0;
    double phiMax = 2.0 * std::acos(-1.0);
    double legUniformFraction = 0.1;
    std::vector<double> safetyFactors = {1.2, 1.5, 2.0};
};

struct SampleBlock {
    std::vector<double> weights;
    double mean = 0.0;
    double maximum = 0.0;
    double q9999 = 0.0;
    double insideFraction = 0.0;
    double milliseconds = 0.0;
};

void printUsage(std::ostream& output) {
    output
        << "Usage: lowmc_pythia8_leg_envelope_diagnostic [output.csv] [options]\n\n"
        << "Calibrate and independently test an empirical global envelope over the\n"
        << "analytic leg-centred proposal. Any trial overshoot is a strict failure.\n\n"
        << "Options:\n"
        << "  --mode ee-hadrons|softqcd\n"
        << "  --events N --first-event N\n"
        << "  --calibration-samples N --trial-samples N\n"
        << "  --seed N --sampling-seed N\n"
        << "  --ecm GEV --kt GEV\n"
        << "  --eta-min X --eta-max X --phi-min X --phi-max X\n"
        << "  --leg-uniform-fraction F\n"
        << "  --safety-factor F       Repeat to replace defaults 1.2, 1.5, 2.0\n"
        << "  --help\n";
}

std::string readStringOption(
    int& index, int argc, char** argv, const std::string& optionName) {
    if (index + 1 >= argc) {
        throw std::invalid_argument(optionName + " requires a value");
    }
    return argv[++index];
}

int parseInteger(
    const std::string& text, const std::string& optionName, bool allowZero) {
    std::size_t parsed = 0;
    const long value = std::stol(text, &parsed);
    if (parsed != text.size() || value < (allowZero ? 0 : 1) ||
        value > std::numeric_limits<int>::max()) {
        throw std::invalid_argument(
            optionName + (allowZero ? " expects a non-negative integer"
                                    : " expects a positive integer"));
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
    bool customSafetyFactors = false;
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument == "--help" || argument == "-h") {
            printUsage(std::cout);
            std::exit(EXIT_SUCCESS);
        }
        if (argument == "--mode") {
            options.mode = readStringOption(index, argc, argv, argument);
        } else if (argument == "--events") {
            options.events = parseInteger(
                readStringOption(index, argc, argv, argument), argument, false);
        } else if (argument == "--first-event") {
            options.firstEvent = parseInteger(
                readStringOption(index, argc, argv, argument), argument, true);
        } else if (argument == "--calibration-samples") {
            options.calibrationSamples = parseInteger(
                readStringOption(index, argc, argv, argument), argument, false);
        } else if (argument == "--trial-samples") {
            options.trialSamples = parseInteger(
                readStringOption(index, argc, argv, argument), argument, false);
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
        } else if (argument == "--safety-factor") {
            if (!customSafetyFactors) {
                options.safetyFactors.clear();
                customSafetyFactors = true;
            }
            options.safetyFactors.push_back(parsePositiveDouble(
                readStringOption(index, argc, argv, argument), argument));
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
    if (std::any_of(
            options.safetyFactors.begin(),
            options.safetyFactors.end(),
            [](double factor) { return factor < 1.0; })) {
        throw std::invalid_argument("--safety-factor must be at least 1");
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
    const double position = probability * static_cast<double>(sorted.size() - 1);
    const std::size_t lower = static_cast<std::size_t>(position);
    const std::size_t upper = std::min(lower + 1, sorted.size() - 1);
    const double fraction = position - static_cast<double>(lower);
    return sorted[lower] + fraction * (sorted[upper] - sorted[lower]);
}

SampleBlock generateWeights(
    const lowmc::LowCurrent& current,
    const lowmc::LegCenteredEtaPhiProposal& proposal,
    double kT,
    int count,
    std::uint64_t seed) {
    std::mt19937_64 rng(seed);
    SampleBlock result;
    result.weights.reserve(static_cast<std::size_t>(count));
    std::size_t inside = 0;
    double sum = 0.0;
    const auto start = Clock::now();
    for (int index = 0; index < count; ++index) {
        const lowmc::FixedKTLegCenteredWeightedSample sample =
            lowmc::makeFixedKTLegCenteredEtaPhiWeightedSample(
                current,
                kT,
                proposal,
                std::generate_canonical<double, 53>(rng),
                std::generate_canonical<double, 53>(rng),
                std::generate_canonical<double, 53>(rng));
        const double weight = sample.weightedPhoton.photon.weight;
        result.weights.push_back(weight);
        sum += weight;
        result.maximum = std::max(result.maximum, weight);
        inside += sample.etaPhiProposal.insideAcceptance ? 1 : 0;
    }
    const auto stop = Clock::now();
    result.mean = sum / static_cast<double>(count);
    result.insideFraction = static_cast<double>(inside) / static_cast<double>(count);
    result.milliseconds =
        std::chrono::duration<double, std::milli>(stop - start).count();
    std::vector<double> sorted = result.weights;
    std::sort(sorted.begin(), sorted.end());
    result.q9999 = percentile(sorted, 0.9999);
    return result;
}

void writeEnvelopeRows(
    std::ofstream& output,
    const Options& options,
    int eventIndex,
    std::uint64_t eventSeed,
    std::size_t selectedCharged,
    double proposalBuildMilliseconds,
    const SampleBlock& calibration,
    const SampleBlock& trials) {
    for (std::size_t factorIndex = 0;
         factorIndex < options.safetyFactors.size(); ++factorIndex) {
        const double safetyFactor = options.safetyFactors[factorIndex];
        const double bound = safetyFactor * calibration.maximum;
        std::size_t overshoots = 0;
        double maximumRatio = 0.0;
        std::size_t accepted = 0;
        std::mt19937_64 acceptanceRng(
            eventSeed + 1000003ULL + static_cast<std::uint64_t>(factorIndex));
        for (double weight : trials.weights) {
            const double ratio = bound > 0.0 ? weight / bound : 0.0;
            maximumRatio = std::max(maximumRatio, ratio);
            overshoots += ratio > 1.0 ? 1 : 0;
            const double acceptanceProbability = std::min(1.0, ratio);
            accepted += std::generate_canonical<double, 53>(acceptanceRng) <
                                acceptanceProbability
                            ? 1
                            : 0;
        }
        const bool envelopeOk = overshoots == 0;
        output << options.mode << ','
               << options.seed << ','
               << eventIndex << ','
               << eventSeed << ','
               << selectedCharged << ','
               << options.legUniformFraction << ','
               << options.calibrationSamples << ','
               << options.trialSamples << ','
               << safetyFactor << ','
               << calibration.mean << ','
               << calibration.maximum << ','
               << calibration.q9999 << ','
               << calibration.insideFraction << ','
               << trials.mean << ','
               << trials.maximum << ','
               << trials.q9999 << ','
               << trials.insideFraction << ','
               << bound << ','
               << maximumRatio << ','
               << overshoots << ','
               << (envelopeOk ? 1 : 0) << ','
               << (bound > 0.0 ? trials.mean / bound : 0.0) << ','
               << static_cast<double>(accepted) /
                      static_cast<double>(options.trialSamples) << ','
               << accepted << ','
               << proposalBuildMilliseconds << ','
               << calibration.milliseconds << ','
               << trials.milliseconds << '\n';
    }
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
        output << "mode,seed,event_index,sampling_seed,selected_charged,"
                  "leg_uniform_fraction,calibration_samples,trial_samples,safety_factor,"
                  "calibration_mean,calibration_max,calibration_q9999,"
                  "calibration_inside_fraction,trial_mean,trial_max,trial_q9999,"
                  "trial_inside_fraction,envelope_bound,max_trial_bound_ratio,overshoots,"
                  "envelope_ok,expected_acceptance,observed_acceptance,accepted_trials,"
                  "proposal_build_milliseconds,calibration_milliseconds,trial_milliseconds\n";
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
            const auto buildStart = Clock::now();
            const lowmc::LegCenteredEtaPhiProposal proposal =
                lowmc::makeLegCenteredEtaPhiProposal(
                    phaseSpace,
                    current.particles(),
                    lowmc::LegCenteredEtaPhiProposalConfig{
                        options.legUniformFraction,
                        1.0e-12,
                    });
            const auto buildStop = Clock::now();
            const double buildMilliseconds =
                std::chrono::duration<double, std::milli>(buildStop - buildStart).count();
            const std::uint64_t eventSeed =
                (options.samplingSeed != 0 ? options.samplingSeed : options.seed) +
                104729ULL * static_cast<std::uint64_t>(eventIndex + 1);
            const SampleBlock calibration = generateWeights(
                current,
                proposal,
                options.kT,
                options.calibrationSamples,
                eventSeed);
            const SampleBlock trials = generateWeights(
                current,
                proposal,
                options.kT,
                options.trialSamples,
                eventSeed + 1ULL);
            writeEnvelopeRows(
                output,
                options,
                eventIndex,
                eventSeed,
                selection.selectedParticles.size(),
                buildMilliseconds,
                calibration,
                trials);
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
