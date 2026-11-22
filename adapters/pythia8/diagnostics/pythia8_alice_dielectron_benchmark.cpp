#include "lowmc/LowMC.h"
#include "lowmc/pythia8/Pythia8Adapter.h"

#include <Pythia8/Pythia.h>

#include <algorithm>
#include <array>
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

constexpr std::array<double, 7> MassEdges{
    0.15, 0.20, 0.25, 0.30, 0.40, 0.50, 0.60};
constexpr std::array<double, 9> QTEdges{
    0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};
constexpr std::array<double, 7> ElectronPTEdges{
    0.075, 0.10, 0.15, 0.20, 0.30, 0.50, 1.0};
constexpr std::array<double, 11> AgreementMassEdges{
    0.02, 0.03, 0.04, 0.05, 0.07, 0.10, 0.15, 0.20, 0.30, 0.40, 0.60};
constexpr std::array<double, 11> AgreementQTEdges{
    0.0, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0};
constexpr double ElectronEtaMax = 0.8;
constexpr double ElectronPTMin = 0.075;
constexpr double ElectronPTMax = 10.0;
// In detector-acceptance mode this is proposal support: requiring both
// leptons inside |eta_e| < ElectronEtaMax already confines the pair rapidity.
// In publication mode it is the explicit provisional midrapidity interval
// used for the acceptance-corrected cross section.
constexpr double SamplingRapidityMin = -ElectronEtaMax;
constexpr double SamplingRapidityMax = ElectronEtaMax;

enum class Method : std::size_t { Sampled = 0, Direct = 1 };
enum class Decay : std::size_t { Polarized = 0, Isotropic = 1 };
enum class Observable : std::size_t { Mass = 0, QT = 1 };
enum class ComparisonMode { DetectorAcceptance, PublishedCrossSection };

struct Options {
    std::string outputPath = "alice_dielectron_prediction.csv";
    int events = 100;
    int samples = 20000;
    int directSamples = 32768;
    std::uint64_t seed = 12345;
    double eCM = 13000.0;
    bool extendedAgreementMap = false;
    ComparisonMode comparisonMode = ComparisonMode::DetectorAcceptance;
    lowmc::VirtualCurrentModel virtualCurrentModel =
        lowmc::VirtualCurrentModel::LeadingLow;
};

struct LogQuantileHistogram {
    static constexpr int BinCount = 240;
    static constexpr double LogMin = -10.0;
    static constexpr double LogMax = 5.0;
    std::array<double, BinCount> weight{};
    double totalWeight = 0.0;

    void add(double value, double sampleWeight) {
        if (!(sampleWeight > 0.0) || !std::isfinite(sampleWeight) || value < 0.0) {
            return;
        }
        double logValue = LogMin;
        if (value > 0.0 && std::isfinite(value)) {
            logValue = std::log10(value);
        } else if (!std::isfinite(value)) {
            logValue = LogMax;
        }
        const double position =
            (std::clamp(logValue, LogMin, LogMax) - LogMin) /
            (LogMax - LogMin);
        const int index = std::min(
            BinCount - 1,
            static_cast<int>(position * static_cast<double>(BinCount)));
        weight[static_cast<std::size_t>(index)] += sampleWeight;
        totalWeight += sampleWeight;
    }

    [[nodiscard]] double quantile(double probability) const {
        if (!(totalWeight > 0.0)) {
            return 0.0;
        }
        const double target = probability * totalWeight;
        double cumulative = 0.0;
        for (int index = 0; index < BinCount; ++index) {
            cumulative += weight[static_cast<std::size_t>(index)];
            if (cumulative >= target) {
                const double center = LogMin +
                    (static_cast<double>(index) + 0.5) *
                    (LogMax - LogMin) / static_cast<double>(BinCount);
                return std::pow(10.0, center);
            }
        }
        return std::pow(10.0, LogMax);
    }
};

struct BinAccumulator {
    double sumEventYield = 0.0;
    double sumEventYield2 = 0.0;
    LogQuantileHistogram softEpsilon;
    LogQuantileHistogram energyFraction;
};

using MassBins = std::array<BinAccumulator, MassEdges.size() - 1>;
using QTBins = std::array<BinAccumulator, QTEdges.size() - 1>;
using ElectronPTBins =
    std::array<BinAccumulator, ElectronPTEdges.size() - 1>;
using AgreementQTBins =
    std::array<BinAccumulator, AgreementQTEdges.size() - 1>;
using AgreementBins =
    std::array<AgreementQTBins, AgreementMassEdges.size() - 1>;

struct Accumulators {
    std::array<std::array<MassBins, 2>, 2> mass;
    std::array<std::array<QTBins, 2>, 2> qT;
    std::array<std::array<ElectronPTBins, 2>, 2> electronPT;
    std::array<std::array<AgreementBins, 2>, 2> agreement;
};

struct EventYields {
    std::array<std::array<std::array<double, MassEdges.size() - 1>, 2>, 2> mass{};
    std::array<std::array<std::array<double, QTEdges.size() - 1>, 2>, 2> qT{};
    std::array<
        std::array<
            std::array<double, ElectronPTEdges.size() - 1>,
            2>,
        2> electronPT{};
    std::array<
        std::array<
            std::array<
                std::array<double, AgreementQTEdges.size() - 1>,
                AgreementMassEdges.size() - 1>,
            2>,
        2> agreement{};
};

void printUsage(std::ostream& output) {
    output
        << "Usage: lowmc_pythia8_alice_dielectron_benchmark [output.csv] [options]\n\n"
        << "Generate non-diffractive 13 TeV pp events and compare ordinary\n"
        << "weighted dielectron sampling with a fixed six-dimensional Sobol\n"
        << "integration on the same events. Detector mode applies the ALICE\n"
        << "single-electron cuts and writes pair and single-electron projections.\n"
        << "Publication mode writes the angle-integrated\n"
        << "pair cross section used for comparison with corrected Fig. 4 data.\n\n"
        << "Options:\n"
        << "  --events N          Accepted Pythia events (default: 100)\n"
        << "  --samples N         Random samples per event (default: 20000)\n"
        << "  --direct-samples N  Sobol samples per event (default: 32768)\n"
        << "  --seed N            Pythia and random-sampling seed\n"
        << "  --ecm GEV           Center-of-mass energy (default: 13000)\n"
        << "  --current MODEL     leading (default) or lichard\n"
        << "  --comparison MODE  detector (default) or publication\n"
        << "  --extended-agreement-map\n"
        << "                      Sample 0.02 < M < 0.60 GeV and pT < 2 GeV\n"
        << "                      for joint current-agreement rows\n"
        << "  --help\n";
}

ComparisonMode parseComparisonMode(const std::string& text) {
    if (text == "detector") {
        return ComparisonMode::DetectorAcceptance;
    }
    if (text == "publication") {
        return ComparisonMode::PublishedCrossSection;
    }
    throw std::invalid_argument(
        "--comparison expects detector or publication");
}

const char* comparisonModeName(ComparisonMode mode) {
    return mode == ComparisonMode::PublishedCrossSection
        ? "published_cross_section"
        : "detector_acceptance";
}

lowmc::VirtualCurrentModel parseVirtualCurrentModel(const std::string& text) {
    if (text == "leading") {
        return lowmc::VirtualCurrentModel::LeadingLow;
    }
    if (text == "lichard") {
        return lowmc::VirtualCurrentModel::LichardEq2_8;
    }
    throw std::invalid_argument("--current expects leading or lichard");
}

const char* virtualCurrentModelName(lowmc::VirtualCurrentModel model) {
    return model == lowmc::VirtualCurrentModel::LichardEq2_8
        ? "lichard_eq2_8"
        : "leading_low";
}

std::string readOption(
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

std::uint64_t parseUInt64(const std::string& text, const std::string& optionName) {
    std::size_t parsed = 0;
    const unsigned long long value = std::stoull(text, &parsed);
    if (parsed != text.size()) {
        throw std::invalid_argument(optionName + " expects an unsigned integer");
    }
    return static_cast<std::uint64_t>(value);
}

double parsePositiveDouble(const std::string& text, const std::string& optionName) {
    std::size_t parsed = 0;
    const double value = std::stod(text, &parsed);
    if (parsed != text.size() || !std::isfinite(value) || value <= 0.0) {
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
        if (argument == "--events") {
            options.events = parsePositiveInt(
                readOption(index, argc, argv, argument), argument);
        } else if (argument == "--samples") {
            options.samples = parsePositiveInt(
                readOption(index, argc, argv, argument), argument);
        } else if (argument == "--direct-samples") {
            options.directSamples = parsePositiveInt(
                readOption(index, argc, argv, argument), argument);
        } else if (argument == "--seed") {
            options.seed = parseUInt64(
                readOption(index, argc, argv, argument), argument);
        } else if (argument == "--ecm") {
            options.eCM = parsePositiveDouble(
                readOption(index, argc, argv, argument), argument);
        } else if (argument == "--current") {
            options.virtualCurrentModel = parseVirtualCurrentModel(
                readOption(index, argc, argv, argument));
        } else if (argument == "--comparison") {
            options.comparisonMode = parseComparisonMode(
                readOption(index, argc, argv, argument));
        } else if (argument == "--extended-agreement-map") {
            options.extendedAgreementMap = true;
        } else if (!argument.empty() && argument.front() == '-') {
            throw std::invalid_argument("unknown argument: " + argument);
        } else if (outputSpecified) {
            throw std::invalid_argument("only one output path may be specified");
        } else {
            options.outputPath = argument;
            outputSpecified = true;
        }
    }
    return options;
}

double pi() {
    return std::acos(-1.0);
}

void initializePythia(Pythia8::Pythia& pythia, const Options& options) {
    pythia.readString("Print:quiet = on");
    pythia.readString("Beams:idA = 2212");
    pythia.readString("Beams:idB = 2212");
    pythia.readString("Beams:eCM = " + std::to_string(options.eCM));
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

class Sobol6 {
public:
    Sobol6() {
        initializeDimension(0, 0, 0, {});
        initializeDimension(1, 1, 0, {1});
        initializeDimension(2, 2, 1, {1, 3});
        initializeDimension(3, 3, 1, {1, 3, 1});
        initializeDimension(4, 3, 2, {1, 1, 1});
        initializeDimension(5, 4, 1, {1, 3, 5, 13});
    }

    std::array<double, 6> next() {
        std::uint64_t value = index_++;
        int bit = 0;
        while ((value & 1U) != 0U) {
            value >>= 1U;
            ++bit;
        }
        if (bit >= 32) {
            throw std::overflow_error("Sobol sequence exhausted");
        }
        std::array<double, 6> point{};
        for (std::size_t dimension = 0; dimension < point.size(); ++dimension) {
            state_[dimension] ^= directions_[dimension][static_cast<std::size_t>(bit)];
            point[dimension] =
                (static_cast<double>(state_[dimension]) + 0.5) / 4294967296.0;
        }
        return point;
    }

private:
    void initializeDimension(
        std::size_t dimension,
        int degree,
        unsigned coefficient,
        std::initializer_list<unsigned> initial) {
        if (dimension == 0) {
            for (int bit = 0; bit < 32; ++bit) {
                directions_[dimension][static_cast<std::size_t>(bit)] =
                    std::uint32_t{1} << (31 - bit);
            }
            return;
        }
        int bit = 0;
        for (const unsigned value : initial) {
            directions_[dimension][static_cast<std::size_t>(bit)] =
                static_cast<std::uint32_t>(value) << (31 - bit);
            ++bit;
        }
        for (; bit < 32; ++bit) {
            std::uint32_t direction =
                directions_[dimension][static_cast<std::size_t>(bit - degree)];
            direction ^= direction >> degree;
            for (int offset = 1; offset < degree; ++offset) {
                if (((coefficient >> (degree - 1 - offset)) & 1U) != 0U) {
                    direction ^=
                        directions_[dimension][static_cast<std::size_t>(bit - offset)];
                }
            }
            directions_[dimension][static_cast<std::size_t>(bit)] = direction;
        }
    }

    std::array<std::array<std::uint32_t, 32>, 6> directions_{};
    std::array<std::uint32_t, 6> state_{};
    std::uint64_t index_ = 0;
};

template <std::size_t N>
int findBin(double value, const std::array<double, N>& edges) {
    if (value < edges.front() || value >= edges.back()) {
        return -1;
    }
    const auto upper = std::upper_bound(edges.begin(), edges.end(), value);
    return static_cast<int>(std::distance(edges.begin(), upper)) - 1;
}

bool passesElectronCuts(const lowmc::FourVector& lepton) {
    const double pT = lepton.pt();
    if (!(pT > ElectronPTMin && pT < ElectronPTMax)) {
        return false;
    }
    const double eta = std::asinh(lepton.pz / pT);
    return std::isfinite(eta) && std::abs(eta) < ElectronEtaMax;
}

std::pair<double, double> softDiagnostics(
    const std::vector<lowmc::ChargedParticle>& particles,
    const lowmc::WeightedDielectronSample& sample) {
    double maximumEpsilon = 0.0;
    double maximumEnergyFraction = 0.0;
    const double mass2 = sample.mass * sample.mass;
    for (const lowmc::ChargedParticle& particle : particles) {
        const double denominator =
            2.0 * std::abs(particle.momentum.dot(sample.virtualPhoton));
        const double epsilon = denominator > 0.0
            ? mass2 / denominator
            : std::numeric_limits<double>::infinity();
        const double energyFraction = particle.momentum.E > 0.0
            ? sample.virtualPhoton.E / particle.momentum.E
            : std::numeric_limits<double>::infinity();
        maximumEpsilon = std::max(maximumEpsilon, epsilon);
        maximumEnergyFraction = std::max(maximumEnergyFraction, energyFraction);
    }
    return {maximumEpsilon, maximumEnergyFraction};
}

void recordSample(
    Method method,
    int sampleCount,
    const std::vector<lowmc::ChargedParticle>& particles,
    const lowmc::WeightedDielectronSample& sample,
    const Options& options,
    EventYields& eventYields,
    Accumulators& accumulators) {
    const bool publicationMode =
        options.comparisonMode == ComparisonMode::PublishedCrossSection;
    if (!publicationMode &&
        (!passesElectronCuts(sample.dielectron.electron) ||
         !passesElectronCuts(sample.dielectron.positron))) {
        return;
    }
    const int massBin = findBin(sample.mass, MassEdges);
    const int qTBin = findBin(sample.qT, QTEdges);
    const int agreementMassBin = findBin(sample.mass, AgreementMassEdges);
    const int agreementQTBin = findBin(sample.qT, AgreementQTEdges);
    const auto [epsilon, energyFraction] = softDiagnostics(particles, sample);
    const std::size_t methodIndex = static_cast<std::size_t>(method);
    const std::array<double, 2> weights{
        sample.dielectron.weight,
        sample.productionDensity * sample.conversionDensityDMass2 /
            sample.proposalDensity,
    };
    // The second weight is angle integrated. In publication mode the sampled
    // decay angles are therefore dummy integration coordinates and only this
    // row is retained. Detector mode keeps both the polarization-aware and
    // isotropic acceptance-folded results.
    const std::size_t firstDecay = publicationMode
        ? static_cast<std::size_t>(Decay::Isotropic)
        : static_cast<std::size_t>(Decay::Polarized);
    for (std::size_t decay = firstDecay; decay < weights.size(); ++decay) {
        const double eventContribution = weights[decay] /
            static_cast<double>(sampleCount);
        if (massBin >= 0 && sample.qT < 0.4) {
            eventYields.mass[methodIndex][decay][static_cast<std::size_t>(massBin)] +=
                eventContribution;
            BinAccumulator& accumulator =
                accumulators.mass[methodIndex][decay][static_cast<std::size_t>(massBin)];
            accumulator.softEpsilon.add(epsilon, weights[decay]);
            accumulator.energyFraction.add(energyFraction, weights[decay]);
        }
        if (qTBin >= 0 && sample.mass >= MassEdges.front() &&
            sample.mass < MassEdges.back()) {
            eventYields.qT[methodIndex][decay][static_cast<std::size_t>(qTBin)] +=
                eventContribution;
            BinAccumulator& accumulator =
                accumulators.qT[methodIndex][decay][static_cast<std::size_t>(qTBin)];
            accumulator.softEpsilon.add(epsilon, weights[decay]);
            accumulator.energyFraction.add(energyFraction, weights[decay]);
        }
        for (const lowmc::FourVector& lepton : {
                 sample.dielectron.electron,
                 sample.dielectron.positron}) {
            const int electronPTBin = findBin(lepton.pt(), ElectronPTEdges);
            if (electronPTBin >= 0) {
                eventYields.electronPT[methodIndex][decay]
                    [static_cast<std::size_t>(electronPTBin)] += eventContribution;
                BinAccumulator& accumulator =
                    accumulators.electronPT[methodIndex][decay]
                        [static_cast<std::size_t>(electronPTBin)];
                accumulator.softEpsilon.add(epsilon, weights[decay]);
                accumulator.energyFraction.add(energyFraction, weights[decay]);
            }
        }
        if (agreementMassBin >= 0 && agreementQTBin >= 0) {
            const std::size_t massIndex =
                static_cast<std::size_t>(agreementMassBin);
            const std::size_t qTIndex =
                static_cast<std::size_t>(agreementQTBin);
            eventYields.agreement[methodIndex][decay][massIndex][qTIndex] +=
                eventContribution;
            BinAccumulator& accumulator =
                accumulators.agreement[methodIndex][decay][massIndex][qTIndex];
            accumulator.softEpsilon.add(epsilon, weights[decay]);
            accumulator.energyFraction.add(energyFraction, weights[decay]);
        }
    }
}

lowmc::DielectronPhaseSpace alicePhaseSpace(const Options& options) {
    return lowmc::DielectronPhaseSpace{
        options.extendedAgreementMap
            ? AgreementMassEdges.front()
            : MassEdges.front(),
        options.extendedAgreementMap
            ? AgreementMassEdges.back()
            : MassEdges.back(),
        options.extendedAgreementMap
            ? AgreementQTEdges.front()
            : QTEdges.front(),
        options.extendedAgreementMap
            ? AgreementQTEdges.back()
            : QTEdges.back(),
        SamplingRapidityMin,
        SamplingRapidityMax,
        0.0,
        2.0 * pi(),
    };
}

void sampleEvent(
    const std::vector<lowmc::ChargedParticle>& particles,
    const Options& options,
    int eventIndex,
    EventYields& eventYields,
    Accumulators& accumulators) {
    const lowmc::LowCurrent current(particles);
    const lowmc::DielectronPhaseSpace phaseSpace = alicePhaseSpace(options);
    std::mt19937_64 engine(
        options.seed + 0x9e3779b97f4a7c15ULL *
            static_cast<std::uint64_t>(eventIndex + 1));
    auto randomUnit = [&engine]() {
        return std::generate_canonical<double, 53>(engine);
    };
    for (int sampleIndex = 0; sampleIndex < options.samples; ++sampleIndex) {
        const lowmc::WeightedDielectronSample sample =
            lowmc::makeUniformMassSquaredQT2DielectronSample(
                current,
                phaseSpace,
                randomUnit(),
                randomUnit(),
                randomUnit(),
                randomUnit(),
                randomUnit(),
                randomUnit(),
                lowmc::ElectronMassGeV,
                lowmc::FineStructureConstant,
                options.virtualCurrentModel);
        recordSample(
            Method::Sampled,
            options.samples,
            particles,
            sample,
            options,
            eventYields,
            accumulators);
    }

    Sobol6 sobol;
    std::mt19937_64 shiftEngine(
        (options.seed ^ 0xd1b54a32d192ed03ULL) +
        0x94d049bb133111ebULL * static_cast<std::uint64_t>(eventIndex + 1));
    std::array<double, 6> sobolShift{};
    for (double& shift : sobolShift) {
        shift = std::generate_canonical<double, 53>(shiftEngine);
    }
    for (int sampleIndex = 0; sampleIndex < options.directSamples; ++sampleIndex) {
        std::array<double, 6> unit = sobol.next();
        for (std::size_t dimension = 0; dimension < unit.size(); ++dimension) {
            unit[dimension] += sobolShift[dimension];
            unit[dimension] -= std::floor(unit[dimension]);
        }
        const lowmc::WeightedDielectronSample sample =
            lowmc::makeUniformMassSquaredQT2DielectronSample(
                current,
                phaseSpace,
                unit[0],
                unit[1],
                unit[2],
                unit[3],
                unit[4],
                unit[5],
                lowmc::ElectronMassGeV,
                lowmc::FineStructureConstant,
                options.virtualCurrentModel);
        recordSample(
            Method::Direct,
            options.directSamples,
            particles,
            sample,
            options,
            eventYields,
            accumulators);
    }
}

void finishEvent(const EventYields& eventYields, Accumulators& accumulators) {
    for (std::size_t method = 0; method < 2; ++method) {
        for (std::size_t decay = 0; decay < 2; ++decay) {
            for (std::size_t bin = 0; bin + 1 < MassEdges.size(); ++bin) {
                const double value = eventYields.mass[method][decay][bin];
                BinAccumulator& accumulator = accumulators.mass[method][decay][bin];
                accumulator.sumEventYield += value;
                accumulator.sumEventYield2 += value * value;
            }
            for (std::size_t bin = 0; bin + 1 < QTEdges.size(); ++bin) {
                const double value = eventYields.qT[method][decay][bin];
                BinAccumulator& accumulator = accumulators.qT[method][decay][bin];
                accumulator.sumEventYield += value;
                accumulator.sumEventYield2 += value * value;
            }
            for (std::size_t bin = 0;
                 bin + 1 < ElectronPTEdges.size(); ++bin) {
                const double value = eventYields.electronPT[method][decay][bin];
                BinAccumulator& accumulator =
                    accumulators.electronPT[method][decay][bin];
                accumulator.sumEventYield += value;
                accumulator.sumEventYield2 += value * value;
            }
            for (std::size_t mass = 0;
                 mass + 1 < AgreementMassEdges.size(); ++mass) {
                for (std::size_t qT = 0;
                     qT + 1 < AgreementQTEdges.size(); ++qT) {
                    const double value =
                        eventYields.agreement[method][decay][mass][qT];
                    BinAccumulator& accumulator =
                        accumulators.agreement[method][decay][mass][qT];
                    accumulator.sumEventYield += value;
                    accumulator.sumEventYield2 += value * value;
                }
            }
        }
    }
}

double standardError(const BinAccumulator& accumulator, int eventCount) {
    if (eventCount < 2) {
        return 0.0;
    }
    const double count = static_cast<double>(eventCount);
    const double centered = std::max(
        0.0,
        accumulator.sumEventYield2 -
            accumulator.sumEventYield * accumulator.sumEventYield / count);
    return std::sqrt(centered / (count * (count - 1.0)));
}

const char* methodName(std::size_t method) {
    return method == static_cast<std::size_t>(Method::Sampled)
        ? "sampled"
        : "direct_sobol";
}

const char* decayName(std::size_t decay, const Options& options) {
    if (options.comparisonMode == ComparisonMode::PublishedCrossSection) {
        return "angle_integrated";
    }
    return decay == static_cast<std::size_t>(Decay::Polarized)
        ? "polarized"
        : "isotropic";
}

std::size_t firstDecayIndex(const Options& options) {
    return options.comparisonMode == ComparisonMode::PublishedCrossSection
        ? static_cast<std::size_t>(Decay::Isotropic)
        : static_cast<std::size_t>(Decay::Polarized);
}

template <std::size_t N>
void writeObservable(
    std::ofstream& output,
    const char* observable,
    const std::array<double, N>& edges,
    const std::array<std::array<std::array<BinAccumulator, N - 1>, 2>, 2>& bins,
    int eventCount,
    int attemptedEvents,
    double sigmaGen,
    const Options& options) {
    const double rapidityWidth =
        SamplingRapidityMax - SamplingRapidityMin;
    for (std::size_t method = 0; method < 2; ++method) {
        for (std::size_t decay = firstDecayIndex(options); decay < 2; ++decay) {
            for (std::size_t bin = 0; bin + 1 < N; ++bin) {
                const BinAccumulator& accumulator = bins[method][decay][bin];
                const double binWidth = edges[bin + 1] - edges[bin];
                const double scale = sigmaGen /
                    (static_cast<double>(eventCount) * rapidityWidth * binWidth);
                output
                    << observable << ','
                    << methodName(method) << ','
                    << decayName(decay, options) << ','
                    << edges[bin] << ','
                    << edges[bin + 1] << ','
                    << 0.5 * (edges[bin] + edges[bin + 1]) << ','
                    << accumulator.sumEventYield * scale << ','
                    << standardError(accumulator, eventCount) *
                        sigmaGen / (rapidityWidth * binWidth) << ','
                    << accumulator.softEpsilon.quantile(0.5) << ','
                    << accumulator.softEpsilon.quantile(0.9) << ','
                    << accumulator.softEpsilon.quantile(0.99) << ','
                    << accumulator.energyFraction.quantile(0.5) << ','
                    << accumulator.energyFraction.quantile(0.9) << ','
                    << accumulator.energyFraction.quantile(0.99) << ','
                    << eventCount << ','
                    << attemptedEvents << ','
                    << (method == 0 ? options.samples : options.directSamples) << ','
                    << sigmaGen << ','
                    << options.seed << ','
                    << virtualCurrentModelName(options.virtualCurrentModel) << ','
                    << comparisonModeName(options.comparisonMode)
                    << ",,,,,,,mb_per_gev_per_y\n";
            }
        }
    }
}

void writeAgreementObservable(
    std::ofstream& output,
    const std::array<std::array<AgreementBins, 2>, 2>& bins,
    int eventCount,
    int attemptedEvents,
    double sigmaGen,
    const Options& options) {
    const double rapidityWidth =
        SamplingRapidityMax - SamplingRapidityMin;
    for (std::size_t method = 0; method < 2; ++method) {
        for (std::size_t decay = firstDecayIndex(options); decay < 2; ++decay) {
            for (std::size_t mass = 0;
                 mass + 1 < AgreementMassEdges.size(); ++mass) {
                const double massWidth =
                    AgreementMassEdges[mass + 1] - AgreementMassEdges[mass];
                for (std::size_t qT = 0;
                     qT + 1 < AgreementQTEdges.size(); ++qT) {
                    const BinAccumulator& accumulator =
                        bins[method][decay][mass][qT];
                    const double qTWidth =
                        AgreementQTEdges[qT + 1] - AgreementQTEdges[qT];
                    const double binArea = massWidth * qTWidth;
                    const double scale = sigmaGen /
                        (static_cast<double>(eventCount) * rapidityWidth * binArea);
                    output
                        << "mass_qT" << ','
                        << methodName(method) << ','
                        << decayName(decay, options) << ','
                        << AgreementMassEdges[mass] << ','
                        << AgreementMassEdges[mass + 1] << ','
                        << 0.5 * (AgreementMassEdges[mass] +
                            AgreementMassEdges[mass + 1]) << ','
                        << accumulator.sumEventYield * scale << ','
                        << standardError(accumulator, eventCount) *
                            sigmaGen / (rapidityWidth * binArea) << ','
                        << accumulator.softEpsilon.quantile(0.5) << ','
                        << accumulator.softEpsilon.quantile(0.9) << ','
                        << accumulator.softEpsilon.quantile(0.99) << ','
                        << accumulator.energyFraction.quantile(0.5) << ','
                        << accumulator.energyFraction.quantile(0.9) << ','
                        << accumulator.energyFraction.quantile(0.99) << ','
                        << eventCount << ','
                        << attemptedEvents << ','
                        << (method == 0 ? options.samples : options.directSamples) << ','
                        << sigmaGen << ','
                        << options.seed << ','
                        << virtualCurrentModelName(options.virtualCurrentModel) << ','
                        << comparisonModeName(options.comparisonMode) << ','
                        << AgreementMassEdges[mass] << ','
                        << AgreementMassEdges[mass + 1] << ','
                        << 0.5 * (AgreementMassEdges[mass] +
                            AgreementMassEdges[mass + 1]) << ','
                        << AgreementQTEdges[qT] << ','
                        << AgreementQTEdges[qT + 1] << ','
                        << 0.5 * (AgreementQTEdges[qT] +
                            AgreementQTEdges[qT + 1]) << ','
                        << "mb_per_gev2_per_y" << '\n';
                }
            }
        }
    }
}

void writeOutput(
    const Options& options,
    const Accumulators& accumulators,
    int eventCount,
    int attemptedEvents,
    double sigmaGen) {
    std::ofstream output(options.outputPath);
    if (!output) {
        throw std::runtime_error("could not open output CSV: " + options.outputPath);
    }
    output << std::setprecision(17);
    output
        << "observable,method,decay,bin_low,bin_high,bin_center,"
        << "cross_section_mb_per_gev_per_y,standard_error_mb_per_gev_per_y,"
        << "max_soft_epsilon_p50,max_soft_epsilon_p90,max_soft_epsilon_p99,"
        << "max_energy_fraction_p50,max_energy_fraction_p90,max_energy_fraction_p99,"
        << "accepted_events,attempted_events,samples_per_event,sigma_gen_mb,seed,"
        << "virtual_current_model,comparison_mode,"
        << "mass_bin_low,mass_bin_high,mass_bin_center,"
        << "qT_bin_low,qT_bin_high,qT_bin_center,density_units\n";
    writeObservable(
        output,
        "mass",
        MassEdges,
        accumulators.mass,
        eventCount,
        attemptedEvents,
        sigmaGen,
        options);
    writeObservable(
        output,
        "qT",
        QTEdges,
        accumulators.qT,
        eventCount,
        attemptedEvents,
        sigmaGen,
        options);
    writeObservable(
        output,
        "electron_pT",
        ElectronPTEdges,
        accumulators.electronPT,
        eventCount,
        attemptedEvents,
        sigmaGen,
        options);
    writeAgreementObservable(
        output,
        accumulators.agreement,
        eventCount,
        attemptedEvents,
        sigmaGen,
        options);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Options options = parseOptions(argc, argv);
        Pythia8::Pythia pythia(LOWMC_PYTHIA8_XMLDOC, false);
        initializePythia(pythia, options);

        Accumulators accumulators;
        int acceptedEvents = 0;
        int attemptedEvents = 0;
        while (acceptedEvents < options.events) {
            ++attemptedEvents;
            if (!pythia.next()) {
                if (attemptedEvents > 20 * options.events) {
                    throw std::runtime_error("too many Pythia event-generation failures");
                }
                continue;
            }
            const lowmc::AdapterParticleSelectionResult selection =
                lowmc::pythia8::selectChargedParticles(pythia.event);
            if (std::abs(selection.diagnostic.signedChargeSum) > 1.0e-12) {
                continue;
            }
            EventYields eventYields;
            sampleEvent(
                selection.selectedParticles,
                options,
                acceptedEvents,
                eventYields,
                accumulators);
            finishEvent(eventYields, accumulators);
            ++acceptedEvents;
        }

        const double sigmaGen = pythia.info.sigmaGen();
        writeOutput(
            options,
            accumulators,
            acceptedEvents,
            attemptedEvents,
            sigmaGen);
        std::cout << "aggregate_status=ok"
                  << " process=SoftQCD:nonDiffractive"
                  << " accepted_events=" << acceptedEvents
                  << " attempted_events=" << attemptedEvents
                  << " samples_per_event=" << options.samples
                  << " direct_samples_per_event=" << options.directSamples
                  << " virtual_current_model="
                  << virtualCurrentModelName(options.virtualCurrentModel)
                  << " comparison_mode="
                  << comparisonModeName(options.comparisonMode)
                  << " sigma_gen_mb=" << sigmaGen
                  << " output=" << options.outputPath << '\n';
    } catch (const std::exception& error) {
        std::cerr << "FAILED: " << error.what() << '\n';
        printUsage(std::cerr);
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
