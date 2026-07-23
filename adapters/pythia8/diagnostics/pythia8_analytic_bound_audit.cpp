#include "lowmc/LowMC.h"
#include "lowmc/pythia8/Pythia8Adapter.h"

#include <Pythia8/Pythia.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
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
    std::string outputPath = "pythia8_analytic_bound_audit.csv";
    std::string mode = "ee-hadrons";
    int events = 2;
    int firstEvent = 0;
    int trials = 10000;
    int calibrationTrials = 1000;
    std::uint64_t seed = 12345;
    std::uint64_t shapeSeed = 0;
    double eCM = 0.0;
};

void printUsage(std::ostream& output) {
    output
        << "Usage: lowmc_pythia8_analytic_bound_audit [output.csv] [options]\n\n"
        << "Audit analytic global and cell-local Low-current weight bounds\n"
        << "alongside the unchanged shape-unweighting envelope on selected\n"
        << "Pythia events. The analytic bounds are never used for sampling.\n\n"
        << "Options:\n"
        << "  --mode ee-hadrons|softqcd\n"
        << "  --events N --first-event N\n"
        << "  --trials N --calibration-trials N\n"
        << "  --seed N --shape-seed N\n"
        << "  --ecm GEV\n"
        << "  --help\n";
}

std::string readOption(
    int& index,
    int argc,
    char** argv,
    const std::string& name) {
    if (index + 1 >= argc) {
        throw std::invalid_argument(name + " requires a value");
    }
    return argv[++index];
}

int parseInteger(
    const std::string& text,
    const std::string& name,
    bool allowZero) {
    std::size_t parsed = 0;
    const long value = std::stol(text, &parsed);
    if (parsed != text.size() ||
        value < (allowZero ? 0 : 1) ||
        value > std::numeric_limits<int>::max()) {
        throw std::invalid_argument(
            name + (allowZero
                        ? " expects a non-negative integer"
                        : " expects a positive integer"));
    }
    return static_cast<int>(value);
}

std::uint64_t parseUInt64(
    const std::string& text,
    const std::string& name) {
    std::size_t parsed = 0;
    const unsigned long long value = std::stoull(text, &parsed);
    if (parsed != text.size()) {
        throw std::invalid_argument(
            name + " expects an unsigned integer");
    }
    return static_cast<std::uint64_t>(value);
}

double parsePositiveDouble(
    const std::string& text,
    const std::string& name) {
    std::size_t parsed = 0;
    const double value = std::stod(text, &parsed);
    if (parsed != text.size() || !std::isfinite(value) ||
        !(value > 0.0)) {
        throw std::invalid_argument(
            name + " expects a positive finite number");
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
            options.mode =
                readOption(index, argc, argv, argument);
        } else if (argument == "--events") {
            options.events = parseInteger(
                readOption(index, argc, argv, argument),
                argument,
                false);
        } else if (argument == "--first-event") {
            options.firstEvent = parseInteger(
                readOption(index, argc, argv, argument),
                argument,
                true);
        } else if (argument == "--trials") {
            options.trials = parseInteger(
                readOption(index, argc, argv, argument),
                argument,
                false);
        } else if (argument == "--calibration-trials") {
            options.calibrationTrials = parseInteger(
                readOption(index, argc, argv, argument),
                argument,
                true);
        } else if (argument == "--seed") {
            options.seed = parseUInt64(
                readOption(index, argc, argv, argument),
                argument);
        } else if (argument == "--shape-seed") {
            options.shapeSeed = parseUInt64(
                readOption(index, argc, argv, argument),
                argument);
        } else if (argument == "--ecm") {
            options.eCM = parsePositiveDouble(
                readOption(index, argc, argv, argument),
                argument);
        } else if (!argument.empty() &&
                   argument.front() == '-') {
            throw std::invalid_argument(
                "unknown argument: " + argument);
        } else if (outputSpecified) {
            throw std::invalid_argument(
                "only one output CSV path can be specified");
        } else {
            options.outputPath = argument;
            outputSpecified = true;
        }
    }
    if (options.mode != "ee-hadrons" &&
        options.mode != "softqcd") {
        throw std::invalid_argument(
            "--mode must be ee-hadrons or softqcd");
    }
    return options;
}

void initializePythia(
    Pythia8::Pythia& pythia,
    const Options& options) {
    pythia.readString("Print:quiet = on");
    pythia.readString("Random:setSeed = on");
    pythia.readString(
        "Random:seed = " + std::to_string(options.seed));
    if (options.mode == "softqcd") {
        pythia.readString("Beams:idA = 2212");
        pythia.readString("Beams:idB = 2212");
        pythia.readString(
            "Beams:eCM = " +
            std::to_string(
                options.eCM > 0.0 ? options.eCM : 13000.0));
        pythia.readString("SoftQCD:nonDiffractive = on");
    } else {
        pythia.readString("PDF:lepton = off");
        pythia.readString("WeakSingleBoson:ffbar2gmZ = on");
        pythia.readString("23:onMode = off");
        pythia.readString("23:onIfAny = 1 2 3 4 5");
        pythia.readString("Beams:idA = 11");
        pythia.readString("Beams:idB = -11");
        pythia.readString(
            "Beams:eCM = " +
            std::to_string(
                options.eCM > 0.0
                    ? options.eCM
                    : pythia.particleData.m0(23)));
    }
    if (!pythia.init()) {
        throw std::runtime_error(
            "Pythia initialization failed");
    }
}

const char* statusName(lowmc::ShapeUnweightingStatus status) {
    switch (status) {
        case lowmc::ShapeUnweightingStatus::Ok:
            return "ok";
        case lowmc::ShapeUnweightingStatus::EnvelopeFailed:
            return "envelope_failed";
        case lowmc::ShapeUnweightingStatus::InvalidSelection:
            return "invalid_selection";
        case lowmc::ShapeUnweightingStatus::NoTrials:
            return "no_trials";
    }
    return "unknown";
}

double maximumOutgoingBeta(
    const std::vector<lowmc::ChargedParticle>& particles) {
    double maximum = 0.0;
    for (const lowmc::ChargedParticle& particle : particles) {
        if (particle.inOutSign != lowmc::InOutSign::Outgoing) {
            continue;
        }
        const double momentum = std::sqrt(
            particle.momentum.px * particle.momentum.px +
            particle.momentum.py * particle.momentum.py +
            particle.momentum.pz * particle.momentum.pz);
        if (particle.momentum.E > 0.0) {
            maximum = std::max(
                maximum,
                momentum / particle.momentum.E);
        }
    }
    return maximum;
}

lowmc::ShapeUnweightingConfig makeShapeConfig(
    const Options& options,
    std::uint64_t eventSeed) {
    lowmc::ShapeUnweightingConfig config;
    config.phaseSpace = lowmc::PhotonPhaseSpace{
        0.002,
        0.2,
        -5.0,
        5.0,
        0.0,
        2.0 * std::acos(-1.0),
    };
    config.trials = options.trials;
    config.seed = eventSeed;
    config.calibrationTrials = options.calibrationTrials;
    config.probePolicy =
        lowmc::ShapeUnweightingProbePolicy::BoostedCells;
    config.failurePolicy =
        lowmc::EnvelopeFailurePolicy::Diagnostic;
    config.auditAnalyticBounds = true;
    return config;
}

void writeRow(
    std::ostream& output,
    const Options& options,
    int eventIndex,
    std::uint64_t eventSeed,
    const lowmc::AdapterParticleSelectionResult& selection,
    const lowmc::ShapeUnweightingResult& result) {
    const lowmc::ShapeUnweightingEnvelopeDiagnostics& envelope =
        result.envelope;
    output
        << options.mode << ','
        << options.seed << ','
        << eventIndex << ','
        << eventSeed << ','
        << selection.selectedParticles.size() << ','
        << selection.diagnostic.selectedIncomingCharged << ','
        << selection.diagnostic.selectedOutgoingCharged << ','
        << maximumOutgoingBeta(selection.selectedParticles) << ','
        << statusName(result.status) << ','
        << envelope.gridCellCount << ','
        << envelope.envelopeIntegral << ','
        << envelope.expectedAcceptanceEfficiency << ','
        << envelope.maxTrialBoundRatio << ','
        << envelope.failedEnvelopeCount << ','
        << (envelope.analyticBoundAvailable ? 1 : 0) << ','
        << envelope.analyticBoundInvalidParticleCount << ','
        << envelope.analyticBoundFallbackCellCount << ','
        << envelope.analyticBoundObservedViolationCount << ','
        << envelope.analyticBoundTrialViolationCount << ','
        << envelope.analyticGlobalWeightBound << ','
        << envelope.analyticCellEnvelopeIntegral << ','
        << envelope.analyticCellEnvelopeOverCurrent << ','
        << envelope.analyticExpectedAcceptanceEfficiency << ','
        << envelope.analyticMaxObservedOverBound << ','
        << envelope.analyticMaxTrialOverBound << ','
        << envelope.analyticBoundPairEvaluations << ','
        << envelope.analyticBoundSeconds << ','
        << envelope.gridSeconds << ','
        << envelope.baseProbeSeconds << ','
        << envelope.calibrationSeconds << ','
        << envelope.boostedProbeSeconds << ','
        << envelope.trialSeconds << '\n';
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Options options = parseOptions(argc, argv);
        std::ofstream output(options.outputPath);
        if (!output) {
            throw std::runtime_error(
                "could not open output CSV: " +
                options.outputPath);
        }
        output
            << "mode,seed,event_index,shape_seed,selected_charged,"
               "selected_incoming,selected_outgoing,max_outgoing_beta,"
               "status,grid_cells,current_envelope_integral,"
               "current_expected_acceptance,max_trial_bound_ratio,"
               "current_envelope_failures,analytic_bound_available,"
               "analytic_invalid_particles,analytic_fallback_cells,"
               "analytic_observed_violations,analytic_trial_violations,"
               "analytic_global_weight_bound,"
               "analytic_cell_envelope_integral,"
               "analytic_cell_envelope_over_current,"
               "analytic_expected_acceptance,"
               "analytic_max_observed_over_bound,"
               "analytic_max_trial_over_bound,"
               "analytic_pair_evaluations,analytic_seconds,"
               "grid_seconds,base_probe_seconds,calibration_seconds,"
               "boosted_probe_seconds,trial_seconds\n";
        output << std::setprecision(17);

        Pythia8::Pythia pythia(
            LOWMC_PYTHIA8_XMLDOC, false);
        initializePythia(pythia, options);
        int validEvents = 0;
        int analyzedEvents = 0;
        int attemptedEvents = 0;
        while (analyzedEvents < options.events) {
            ++attemptedEvents;
            if (!pythia.next()) {
                if (attemptedEvents >
                    (options.firstEvent + options.events) * 20) {
                    throw std::runtime_error(
                        "too many Pythia generation failures");
                }
                continue;
            }
            const lowmc::AdapterParticleSelectionResult selection =
                lowmc::pythia8::selectChargedParticles(
                    pythia.event);
            if (selection.selectedParticles.size() < 2 ||
                std::abs(
                    selection.diagnostic.signedChargeSum) >
                    1.0e-12) {
                continue;
            }
            const int eventIndex = validEvents++;
            if (eventIndex < options.firstEvent) {
                continue;
            }
            const std::uint64_t eventSeed =
                (options.shapeSeed != 0
                     ? options.shapeSeed
                     : options.seed) +
                104729ULL *
                    static_cast<std::uint64_t>(eventIndex + 1);
            const lowmc::ShapeUnweightingResult result =
                lowmc::generateShapeUnweightedPhotonCandidates(
                    selection.selectedParticles,
                    makeShapeConfig(options, eventSeed));
            writeRow(
                output,
                options,
                eventIndex,
                eventSeed,
                selection,
                result);
            ++analyzedEvents;
        }
        if (!output) {
            throw std::runtime_error(
                "failed while writing analytic-bound audit CSV");
        }
        std::cout
            << "analytic_bound_audit_status=ok"
            << " mode=" << options.mode
            << " seed=" << options.seed
            << " events=" << analyzedEvents
            << " output=" << options.outputPath
            << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "FAILED: " << error.what() << "\n\n";
        printUsage(std::cerr);
        return EXIT_FAILURE;
    }
}
