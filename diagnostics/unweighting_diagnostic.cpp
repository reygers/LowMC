#include "lowmc/LowMC.h"

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

enum class EnvelopeMode {
    Global,
    Cell,
    ProbedCell,
};

enum class ProbeRule {
    Interior,
    Edge,
};

struct Options {
    int calibrationSamples = 1000;
    int trialSamples = 1000;
    std::uint64_t seed = 12345;
    double safetyFactor = 1.5;
    EnvelopeMode envelopeMode = EnvelopeMode::Global;
    int probeBins = 11;
    ProbeRule probeRule = ProbeRule::Interior;
    double momentum = 100.0;
    double particleEta = 0.0;
    double particlePhi = 0.5 * std::acos(-1.0);
    double kTMin = 0.002;
    double kTMax = 0.2;
    double etaMin = -5.0;
    double etaMax = 5.0;
    double phiMin = 0.0;
    double phiMax = 2.0 * std::acos(-1.0);
};

struct SampleRecord {
    lowmc::WeightedPhotonSample weightedPhoton;
    std::size_t cellIndex = 0;
};

void printUsage(std::ostream& out, const char* programName) {
    out << "Usage: " << programName << " [options]\n"
        << "\n"
        << "Probe accept/reject unweighting for weighted Low photons in a two-pion toy event.\n"
        << "This diagnostic does not write HepMC3 output events.\n"
        << "\n"
        << "Options:\n"
        << "  --calibration-samples N  Samples used to estimate max weight (default: 1000)\n"
        << "  --trial-samples N        Fresh samples used to test the envelope (default: 1000)\n"
        << "  --seed N                 Random seed (default: 12345)\n"
        << "  --safety-factor X        Envelope multiplier on calibration max (default: 1.5)\n"
        << "  --envelope-mode MODE     Envelope mode: global, cell, or probed-cell (default: global)\n"
        << "  --probe-bins N           Probe lattice bins per cell dimension (default: 11)\n"
        << "  --probe-rule RULE        Probe rule: interior or edge (default: interior)\n"
        << "  --momentum GeV           Charged-pion momentum (default: 100)\n"
        << "  --particle-eta X         First charged-pion pseudorapidity (default: 0)\n"
        << "  --particle-phi X         First charged-pion azimuth (default: pi/2)\n"
        << "  --kt-min GeV             Minimum photon k_T (default: 0.002)\n"
        << "  --kt-max GeV             Maximum photon k_T (default: 0.2)\n"
        << "  --eta-min X              Minimum photon eta (default: -5)\n"
        << "  --eta-max X              Maximum photon eta (default: 5)\n"
        << "  --phi-min X              Minimum photon phi (default: 0)\n"
        << "  --phi-max X              Maximum photon phi (default: 2 pi)\n";
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

EnvelopeMode parseEnvelopeMode(const std::string& text) {
    if (text == "global") {
        return EnvelopeMode::Global;
    }
    if (text == "cell") {
        return EnvelopeMode::Cell;
    }
    if (text == "probed-cell") {
        return EnvelopeMode::ProbedCell;
    }
    throw std::invalid_argument("--envelope-mode must be global, cell, or probed-cell");
}

const char* envelopeModeName(EnvelopeMode mode) {
    switch (mode) {
        case EnvelopeMode::Global:
            return "global";
        case EnvelopeMode::Cell:
            return "cell";
        case EnvelopeMode::ProbedCell:
            return "probed-cell";
    }
    return "unknown";
}

ProbeRule parseProbeRule(const std::string& text) {
    if (text == "interior") {
        return ProbeRule::Interior;
    }
    if (text == "edge") {
        return ProbeRule::Edge;
    }
    throw std::invalid_argument("--probe-rule must be interior or edge");
}

const char* probeRuleName(ProbeRule rule) {
    switch (rule) {
        case ProbeRule::Interior:
            return "interior";
        case ProbeRule::Edge:
            return "edge";
    }
    return "unknown";
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
            printUsage(std::cout, argv[0]);
            std::exit(EXIT_SUCCESS);
        }
        if (argument == "--calibration-samples") {
            options.calibrationSamples =
                parsePositiveInt(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--trial-samples") {
            options.trialSamples =
                parsePositiveInt(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--seed") {
            options.seed = parseUInt64(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--safety-factor") {
            options.safetyFactor =
                parsePositiveDouble(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--envelope-mode") {
            options.envelopeMode =
                parseEnvelopeMode(readStringOption(index, argc, argv, argument));
            continue;
        }
        if (argument == "--probe-bins") {
            options.probeBins =
                parsePositiveInt(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--probe-rule") {
            options.probeRule =
                parseProbeRule(readStringOption(index, argc, argv, argument));
            continue;
        }
        if (argument == "--momentum") {
            options.momentum =
                parsePositiveDouble(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--particle-eta") {
            options.particleEta =
                parseDouble(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--particle-phi") {
            options.particlePhi =
                parseDouble(readStringOption(index, argc, argv, argument), argument);
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
        throw std::invalid_argument("unknown option: " + argument);
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

std::vector<lowmc::ChargedParticle> makeBackToBackPions(
    double momentum,
    double particleEta,
    double particlePhi) {
    const double mass = 0.13957039;
    const double pT = momentum / std::cosh(particleEta);
    const double px = pT * std::cos(particlePhi);
    const double py = pT * std::sin(particlePhi);
    const double pz = pT * std::sinh(particleEta);
    const double energy = std::sqrt(momentum * momentum + mass * mass);
    return {
        {lowmc::FourVector{energy, px, py, pz}, +1.0, lowmc::InOutSign::Outgoing},
        {lowmc::FourVector{energy, -px, -py, -pz}, -1.0, lowmc::InOutSign::Outgoing},
    };
}

std::size_t findCellIndex(
    const lowmc::GridAdaptiveEtaPhiProposal& proposal,
    double eta,
    double phi) {
    for (std::size_t index = 0; index < proposal.cells.size(); ++index) {
        const lowmc::GridAdaptiveEtaPhiCell& cell = proposal.cells[index];
        if (eta >= cell.etaMin && eta <= cell.etaMax &&
            phi >= cell.phiMin && phi <= cell.phiMax) {
            return index;
        }
    }
    throw std::runtime_error("sampled eta-phi point is outside all grid cells");
}

std::vector<SampleRecord> generateToyRecords(
    const lowmc::LowCurrent& current,
    const lowmc::PhotonPhaseSpace& phaseSpace,
    const lowmc::GridAdaptiveEtaPhiProposal& gridProposal,
    int samples,
    std::uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::vector<SampleRecord> records;
    records.reserve(static_cast<std::size_t>(samples));

    for (int index = 0; index < samples; ++index) {
        const double unitKT = std::generate_canonical<double, 53>(rng);
        const double unitEta = std::generate_canonical<double, 53>(rng);
        const double unitPhi = std::generate_canonical<double, 53>(rng);
        const double unitComponent = std::generate_canonical<double, 53>(rng);
        const double unitCell = std::generate_canonical<double, 53>(rng);

        const lowmc::GridAdaptiveWeightedPhotonSample sample =
            lowmc::makeGridAdaptiveWeightedPhotonSample(
                current,
                phaseSpace,
                gridProposal,
                unitKT,
                unitComponent,
                unitCell,
                unitEta,
                unitPhi);

        records.push_back(SampleRecord{
            sample.weightedPhoton,
            findCellIndex(gridProposal, sample.etaPhiProposal.eta, sample.etaPhiProposal.phi),
        });
    }

    return records;
}

double probeCoordinate(double lower, double upper, int index, int bins, ProbeRule rule) {
    if (rule == ProbeRule::Interior) {
        return lower + (upper - lower) *
                           (static_cast<double>(index) + 0.5) /
                           static_cast<double>(bins);
    }

    if (bins == 1) {
        return 0.5 * (lower + upper);
    }
    const double fraction = static_cast<double>(index) / static_cast<double>(bins - 1);
    if (index == 0) {
        return std::nextafter(lower, upper);
    }
    if (index == bins - 1) {
        return std::nextafter(upper, lower);
    }
    return lower + (upper - lower) * fraction;
}

double exactWeightAtEtaPhi(
    const lowmc::LowCurrent& current,
    const lowmc::PhotonPhaseSpace& phaseSpace,
    const lowmc::GridAdaptiveEtaPhiProposal& gridProposal,
    double kT,
    double eta,
    double phi) {
    const double proposalDensity =
        lowmc::logUniformKTDensityPhysicalMeasure(phaseSpace, kT) *
        lowmc::gridAdaptiveEtaPhiDensity(gridProposal, eta, phi);
    const lowmc::WeightedPhotonSample sample =
        lowmc::makeWeightedPhotonSample(
            current,
            phaseSpace,
            lowmc::WeightedPhotonProposal{
                lowmc::photonMomentumFromKTEtaPhi(kT, eta, phi),
                proposalDensity,
            });
    return sample.photon.weight;
}

std::vector<double> probeCellMaxWeights(
    const lowmc::LowCurrent& current,
    const lowmc::PhotonPhaseSpace& phaseSpace,
    const lowmc::GridAdaptiveEtaPhiProposal& gridProposal,
    int probeBins,
    ProbeRule probeRule) {
    const double kT = std::sqrt(phaseSpace.kTMin * phaseSpace.kTMax);
    std::vector<double> cellMaxWeights(gridProposal.cells.size(), 0.0);
    for (std::size_t cellIndex = 0; cellIndex < gridProposal.cells.size(); ++cellIndex) {
        const lowmc::GridAdaptiveEtaPhiCell& cell = gridProposal.cells[cellIndex];
        double maxWeight = 0.0;
        for (int etaIndex = 0; etaIndex < probeBins; ++etaIndex) {
            const double eta =
                probeCoordinate(cell.etaMin, cell.etaMax, etaIndex, probeBins, probeRule);
            for (int phiIndex = 0; phiIndex < probeBins; ++phiIndex) {
                const double phi =
                    probeCoordinate(cell.phiMin, cell.phiMax, phiIndex, probeBins, probeRule);
                maxWeight = std::max(
                    maxWeight,
                    exactWeightAtEtaPhi(current, phaseSpace, gridProposal, kT, eta, phi));
            }
        }
        cellMaxWeights[cellIndex] = maxWeight;
    }
    return cellMaxWeights;
}

lowmc::WeightSummary summarizeRecords(const std::vector<SampleRecord>& records) {
    std::vector<lowmc::WeightedPhotonSample> samples;
    samples.reserve(records.size());
    for (const SampleRecord& record : records) {
        samples.push_back(record.weightedPhoton);
    }
    return lowmc::summarizeWeightedPhotonSamples(samples);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Options options = parseOptions(argc, argv);

        const lowmc::PhotonPhaseSpace phaseSpace{
            options.kTMin,
            options.kTMax,
            options.etaMin,
            options.etaMax,
            options.phiMin,
            options.phiMax,
        };
        const lowmc::LowCurrent current(
            makeBackToBackPions(options.momentum, options.particleEta, options.particlePhi));
        const lowmc::GridAdaptiveEtaPhiProposal gridProposal =
            lowmc::makeGridAdaptiveEtaPhiProposal(
                current,
                phaseSpace,
                std::sqrt(phaseSpace.kTMin * phaseSpace.kTMax),
                lowmc::tunedMassGridAdaptiveEtaPhiConfig());

        const std::vector<SampleRecord> calibration =
            generateToyRecords(
                current,
                phaseSpace,
                gridProposal,
                options.calibrationSamples,
                options.seed);
        const lowmc::WeightSummary calibrationSummary = summarizeRecords(calibration);
        const double calibrationMaxWeight = calibrationSummary.maxWeight;
        const double globalEnvelopeWeight = options.safetyFactor * calibrationMaxWeight;

        std::vector<double> calibrationCellMaxWeights(gridProposal.cells.size(), 0.0);
        std::vector<int> calibrationCellCounts(gridProposal.cells.size(), 0);
        for (const SampleRecord& record : calibration) {
            calibrationCellMaxWeights[record.cellIndex] =
                std::max(calibrationCellMaxWeights[record.cellIndex], record.weightedPhoton.photon.weight);
            ++calibrationCellCounts[record.cellIndex];
        }

        int calibrationEmptyCells = 0;
        for (const int count : calibrationCellCounts) {
            if (count == 0) {
                ++calibrationEmptyCells;
            }
        }

        const std::vector<double> probedCellMaxWeights =
            probeCellMaxWeights(
                current,
                phaseSpace,
                gridProposal,
                options.probeBins,
                options.probeRule);
        const std::size_t probeEvaluationCount =
            gridProposal.cells.size() *
            static_cast<std::size_t>(options.probeBins) *
            static_cast<std::size_t>(options.probeBins);
        int suspectCellCount = 0;
        for (std::size_t index = 0; index < gridProposal.cells.size(); ++index) {
            if (calibrationCellMaxWeights[index] > probedCellMaxWeights[index]) {
                ++suspectCellCount;
            }
        }

        const std::vector<SampleRecord> trial =
            generateToyRecords(
                current,
                phaseSpace,
                gridProposal,
                options.trialSamples,
                options.seed + 1);
        const lowmc::WeightSummary trialSummary = summarizeRecords(trial);

        int acceptedPhotons = 0;
        int failedEnvelopeCount = 0;
        int fallbackEnvelopeCount = 0;
        double maxObservedWeight = 0.0;
        double maxLocalRatio = 0.0;
        double cappedAcceptanceSum = 0.0;

        std::mt19937_64 rng(options.seed + 2);
        for (const SampleRecord& sample : trial) {
            const double weight = sample.weightedPhoton.photon.weight;
            maxObservedWeight = std::max(maxObservedWeight, weight);

            double envelopeWeight = globalEnvelopeWeight;
            if (options.envelopeMode == EnvelopeMode::Cell) {
                const double cellMax = calibrationCellMaxWeights[sample.cellIndex];
                if (cellMax > 0.0) {
                    envelopeWeight = options.safetyFactor * cellMax;
                    maxLocalRatio = std::max(maxLocalRatio, weight / cellMax);
                } else {
                    ++fallbackEnvelopeCount;
                }
            } else if (options.envelopeMode == EnvelopeMode::ProbedCell) {
                const double cellMax = std::max(
                    probedCellMaxWeights[sample.cellIndex],
                    calibrationCellMaxWeights[sample.cellIndex]);
                if (cellMax > 0.0) {
                    envelopeWeight = options.safetyFactor * cellMax;
                    maxLocalRatio = std::max(maxLocalRatio, weight / cellMax);
                } else {
                    ++fallbackEnvelopeCount;
                }
            }

            if (weight > envelopeWeight) {
                ++failedEnvelopeCount;
            }

            const double acceptanceProbability =
                envelopeWeight > 0.0 ? std::min(1.0, weight / envelopeWeight) : 0.0;
            cappedAcceptanceSum += acceptanceProbability;

            const double unit = std::generate_canonical<double, 53>(rng);
            if (unit < acceptanceProbability) {
                ++acceptedPhotons;
            }
        }

        const double trialCount = static_cast<double>(trial.size());
        const double observedAcceptance = trial.empty()
            ? 0.0
            : static_cast<double>(acceptedPhotons) / trialCount;
        const double expectedAcceptance = trial.empty() ? 0.0 : cappedAcceptanceSum / trialCount;

        std::cout << "calibration_samples=" << calibrationSummary.count
                  << " trial_samples=" << trialSummary.count
                  << " momentum=" << options.momentum
                  << " particle_eta=" << options.particleEta
                  << " particle_phi=" << options.particlePhi
                  << " safety_factor=" << options.safetyFactor
                  << " envelope_mode=" << envelopeModeName(options.envelopeMode)
                  << " probe_bins=" << options.probeBins
                  << " probe_rule=" << probeRuleName(options.probeRule)
                  << " grid_cells=" << gridProposal.cells.size()
                  << " probe_evaluations=" << probeEvaluationCount
                  << " calibration_empty_cells=" << calibrationEmptyCells
                  << " suspect_cell_count=" << suspectCellCount
                  << " calibration_max_weight=" << calibrationMaxWeight
                  << " global_envelope_weight=" << globalEnvelopeWeight
                  << " trial_max_weight=" << maxObservedWeight
                  << " trial_calibration_max_ratio="
                  << (calibrationMaxWeight > 0.0 ? maxObservedWeight / calibrationMaxWeight : 0.0)
                  << " max_local_trial_calibration_ratio=" << maxLocalRatio
                  << " fallback_envelope_count=" << fallbackEnvelopeCount
                  << " failed_envelope_count=" << failedEnvelopeCount
                  << " envelope_status=" << (failedEnvelopeCount == 0 ? "held" : "failed")
                  << " accepted_photons=" << acceptedPhotons
                  << " observed_acceptance_efficiency=" << observedAcceptance
                  << " expected_acceptance_efficiency=" << expectedAcceptance
                  << '\n';
    } catch (const std::exception& error) {
        std::cerr << "FAILED: " << error.what() << "\n\n";
        printUsage(std::cerr, argv[0]);
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
