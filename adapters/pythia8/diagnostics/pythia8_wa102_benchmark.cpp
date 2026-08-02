#include "lowmc/LowMC.h"
#include "lowmc/pythia8/Pythia8Adapter.h"

#include <Pythia8/Pythia.h>

#include <algorithm>
#include <array>
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

constexpr int BinCount = 10;
constexpr int StageCount = 4;
constexpr int TruthStage = 0;
constexpr int IsolatedStage = 1;
constexpr int TruthSmearedStage = 2;
constexpr int SmearedStage = 3;
constexpr double BeamMomentumGeV = 450.0;
constexpr double EnergyMin = 0.2;
constexpr double EnergyMax = 1.0;
constexpr double ThetaMax = 0.020;
constexpr double ThetaMinForIntegration = 1.0e-5;
constexpr double PtBinWidth = 0.002;
constexpr double PtMax = BinCount * PtBinWidth;
constexpr double MagneticConversion = 0.299792458;
constexpr int ComparisonFirstBin = 1;
constexpr int ComparisonLastBinExclusive = 9;

using StageBins = std::array<std::array<double, BinCount>, StageCount>;

struct Options {
    std::string outputPath = "pythia8_wa102_benchmark.csv";
    int events = 2;
    int firstEvent = 0;
    int mcSamples = 200000;
    int directLogEBins = 24;
    int directEtaBins = 32;
    int directPhiBins = 32;
    int directRefinementLevels = 2;
    int directSubBins = 2;
    std::uint64_t seed = 12345;
    std::uint64_t samplingSeed = 0;
    double legUniformFraction = 0.1;
    double leadDistanceM = 0.66;
    double isolationDistanceM = 0.003;
    double magneticFieldT = 1.2;
    double sheetHalfSizeM = 0.25;
    double ptResolutionGeV = 0.0029;
};

struct SheetHit {
    double x = 0.0;
    double y = 0.0;
};

struct EventSelection {
    bool accepted = false;
    int chargedMultiplicity = 0;
    std::vector<SheetHit> trackHits;
};

struct Estimate {
    StageBins yield{};
    StageBins standardError{};
    std::array<double, StageCount> totalStandardError{};
    std::array<double, StageCount> comparisonStandardError{};
    int evaluations = 0;
    int gridCells = 0;
    double milliseconds = 0.0;
};

double pi() {
    return std::acos(-1.0);
}

void printUsage(std::ostream& output) {
    output
        << "Usage: lowmc_pythia8_wa102_benchmark [output.csv] [options]\n\n"
        << "Compare direct Low integration and weighted LowMC sampling for\n"
        << "WA102 fixed-target pp events. The diagnostic records truth,\n"
        << "approximate 66 cm / 3 mm isolation, and isolated-plus-smeared\n"
        << "spectra. It does not emit photons.\n\n"
        << "Published photon cuts:\n"
        << "  beam momentum = 450 GeV/c, 0.2 < E_gamma < 1 GeV,\n"
        << "  theta_gamma < 20 mrad, ten 2 MeV/c pT bins.\n\n"
        << "Options:\n"
        << "  --events N --first-event N --mc-samples N\n"
        << "  --direct-loge-bins N\n"
        << "  --direct-eta-bins N --direct-phi-bins N\n"
        << "  --direct-refinement-levels N --direct-subbins N\n"
        << "  --seed N --sampling-seed N\n"
        << "  --leg-uniform-fraction F\n"
        << "  --lead-distance-m L       (default 0.66)\n"
        << "  --isolation-distance-m D  (default 0.003)\n"
        << "  --magnetic-field-t B      (default 1.2)\n"
        << "  --sheet-half-size-m H     (default 0.25)\n"
        << "  --pt-resolution-gev S     (default 0.0029)\n"
        << "  --help\n";
}

std::string readOption(
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

std::uint64_t parseUInt64(
    const std::string& text, const std::string& optionName) {
    std::size_t parsed = 0;
    const unsigned long long value = std::stoull(text, &parsed);
    if (parsed != text.size()) {
        throw std::invalid_argument(
            optionName + " expects an unsigned integer");
    }
    return static_cast<std::uint64_t>(value);
}

double parseFinite(
    const std::string& text,
    const std::string& optionName,
    bool allowZero) {
    std::size_t parsed = 0;
    const double value = std::stod(text, &parsed);
    if (parsed != text.size() || !std::isfinite(value) ||
        value < 0.0 || (!allowZero && value == 0.0)) {
        throw std::invalid_argument(
            optionName + (allowZero ? " expects a non-negative number"
                                    : " expects a positive number"));
    }
    return value;
}

Options parseOptions(int argc, char** argv) {
    Options options;
    bool outputSet = false;
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument == "--help") {
            printUsage(std::cout);
            std::exit(EXIT_SUCCESS);
        }
        if (argument == "--events") {
            options.events = parseInteger(
                readOption(index, argc, argv, argument), argument, false);
        } else if (argument == "--first-event") {
            options.firstEvent = parseInteger(
                readOption(index, argc, argv, argument), argument, true);
        } else if (argument == "--mc-samples") {
            options.mcSamples = parseInteger(
                readOption(index, argc, argv, argument), argument, false);
        } else if (argument == "--direct-loge-bins") {
            options.directLogEBins = parseInteger(
                readOption(index, argc, argv, argument), argument, false);
        } else if (argument == "--direct-eta-bins") {
            options.directEtaBins = parseInteger(
                readOption(index, argc, argv, argument), argument, false);
        } else if (argument == "--direct-phi-bins") {
            options.directPhiBins = parseInteger(
                readOption(index, argc, argv, argument), argument, false);
        } else if (argument == "--direct-refinement-levels") {
            options.directRefinementLevels = parseInteger(
                readOption(index, argc, argv, argument), argument, true);
        } else if (argument == "--direct-subbins") {
            options.directSubBins = parseInteger(
                readOption(index, argc, argv, argument), argument, false);
        } else if (argument == "--seed") {
            options.seed =
                parseUInt64(readOption(index, argc, argv, argument), argument);
        } else if (argument == "--sampling-seed") {
            options.samplingSeed =
                parseUInt64(readOption(index, argc, argv, argument), argument);
        } else if (argument == "--leg-uniform-fraction") {
            options.legUniformFraction = parseFinite(
                readOption(index, argc, argv, argument), argument, false);
            if (!(options.legUniformFraction < 1.0)) {
                throw std::invalid_argument(
                    argument + " expects a value below one");
            }
        } else if (argument == "--lead-distance-m") {
            options.leadDistanceM = parseFinite(
                readOption(index, argc, argv, argument), argument, false);
        } else if (argument == "--isolation-distance-m") {
            options.isolationDistanceM = parseFinite(
                readOption(index, argc, argv, argument), argument, true);
        } else if (argument == "--magnetic-field-t") {
            options.magneticFieldT = parseFinite(
                readOption(index, argc, argv, argument), argument, true);
        } else if (argument == "--sheet-half-size-m") {
            options.sheetHalfSizeM = parseFinite(
                readOption(index, argc, argv, argument), argument, false);
        } else if (argument == "--pt-resolution-gev") {
            options.ptResolutionGeV = parseFinite(
                readOption(index, argc, argv, argument), argument, true);
        } else if (!argument.empty() && argument.front() == '-') {
            throw std::invalid_argument("unknown option: " + argument);
        } else if (!outputSet) {
            options.outputPath = argument;
            outputSet = true;
        } else {
            throw std::invalid_argument(
                "unexpected positional argument: " + argument);
        }
    }
    return options;
}

double thetaToEta(double theta) {
    return -std::log(std::tan(0.5 * theta));
}

lowmc::PhotonPhaseSpace angularPhaseSpace() {
    return lowmc::PhotonPhaseSpace{
        1.0e-8,
        PtMax,
        thetaToEta(ThetaMax),
        thetaToEta(ThetaMinForIntegration),
        0.0,
        2.0 * pi(),
    };
}

lowmc::FourVector photonFromEnergyEtaPhi(
    double energy, double eta, double phi) {
    const double inverseCosh = 1.0 / std::cosh(eta);
    return lowmc::FourVector{
        energy,
        energy * inverseCosh * std::cos(phi),
        energy * inverseCosh * std::sin(phi),
        energy * std::tanh(eta),
    };
}

void initializePythia(Pythia8::Pythia& pythia, const Options& options) {
    const double protonMass = pythia.particleData.m0(2212);
    const double beamEnergy =
        std::sqrt(BeamMomentumGeV * BeamMomentumGeV + protonMass * protonMass);
    pythia.readString("Print:quiet = on");
    pythia.readString("Random:setSeed = on");
    pythia.readString("Random:seed = " + std::to_string(options.seed));
    pythia.readString("Beams:idA = 2212");
    pythia.readString("Beams:idB = 2212");
    pythia.readString("Beams:frameType = 2");
    pythia.readString("Beams:eA = " + std::to_string(beamEnergy));
    pythia.readString("Beams:eB = " + std::to_string(protonMass));
    pythia.readString("SoftQCD:inelastic = on");
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

double chargedTrackZ(
    double path,
    double px,
    double pz,
    double momentum,
    double omega) {
    if (std::abs(omega) < 1.0e-14) {
        return path * pz / momentum;
    }
    return (
        px * (1.0 - std::cos(omega * path)) +
        pz * std::sin(omega * path)) /
        (momentum * omega);
}

SheetHit chargedTrackPosition(
    double path,
    double px,
    double py,
    double pz,
    double momentum,
    double omega) {
    if (std::abs(omega) < 1.0e-14) {
        return SheetHit{
            path * px / momentum,
            path * py / momentum,
        };
    }
    return SheetHit{
        (
            px * std::sin(omega * path) +
            pz * (std::cos(omega * path) - 1.0)) /
            (momentum * omega),
        path * py / momentum,
    };
}

bool trackSheetHit(
    const Pythia8::Particle& particle,
    const Options& options,
    SheetHit& hit) {
    const double momentum = particle.pAbs();
    if (!(momentum > 0.0) || !std::isfinite(momentum)) {
        return false;
    }
    const double omega =
        particle.charge() * MagneticConversion * options.magneticFieldT /
        momentum;
    constexpr double pathStep = 0.001;
    constexpr double maxPath = 5.0;
    double lower = 0.0;
    double lowerZ = 0.0;
    for (double upper = pathStep; upper <= maxPath; upper += pathStep) {
        const double upperZ = chargedTrackZ(
            upper,
            particle.px(),
            particle.pz(),
            momentum,
            omega);
        if (lowerZ < options.leadDistanceM &&
            upperZ >= options.leadDistanceM) {
            double left = lower;
            double right = upper;
            for (int iteration = 0; iteration < 48; ++iteration) {
                const double middle = 0.5 * (left + right);
                const double middleZ = chargedTrackZ(
                    middle,
                    particle.px(),
                    particle.pz(),
                    momentum,
                    omega);
                if (middleZ < options.leadDistanceM) {
                    left = middle;
                } else {
                    right = middle;
                }
            }
            hit = chargedTrackPosition(
                0.5 * (left + right),
                particle.px(),
                particle.py(),
                particle.pz(),
                momentum,
                omega);
            return std::abs(hit.x) <= options.sheetHalfSizeM &&
                   std::abs(hit.y) <= options.sheetHalfSizeM;
        }
        lower = upper;
        lowerZ = upperZ;
    }
    return false;
}

EventSelection selectEvent(
    const Pythia8::Event& event, const Options& options) {
    EventSelection result;
    for (int index = 0; index < event.size(); ++index) {
        const Pythia8::Particle& particle = event[index];
        if (!particle.isFinal() || !particle.isCharged()) {
            continue;
        }
        ++result.chargedMultiplicity;
        SheetHit hit;
        if (trackSheetHit(particle, options, hit)) {
            result.trackHits.push_back(hit);
        }
    }
    result.accepted =
        result.chargedMultiplicity > 0 &&
        result.chargedMultiplicity < 8;
    return result;
}

int ptBin(double pt) {
    if (!(pt >= 0.0 && pt < PtMax)) {
        return -1;
    }
    return std::min(BinCount - 1, static_cast<int>(pt / PtBinWidth));
}

bool passesIsolation(
    double nx,
    double ny,
    double nz,
    const std::vector<SheetHit>& trackHits,
    const Options& options) {
    if (options.isolationDistanceM == 0.0) {
        return true;
    }
    if (!(nz > 0.0)) {
        return false;
    }
    const double photonX = options.leadDistanceM * nx / nz;
    const double photonY = options.leadDistanceM * ny / nz;
    if (std::abs(photonX) > options.sheetHalfSizeM ||
        std::abs(photonY) > options.sheetHalfSizeM) {
        return false;
    }
    const double isolation2 =
        options.isolationDistanceM * options.isolationDistanceM;
    for (const SheetHit& hit : trackHits) {
        const double dx = photonX - hit.x;
        const double dy = photonY - hit.y;
        if (dx * dx + dy * dy <= isolation2) {
            return false;
        }
    }
    return true;
}

double normalCdf(double value) {
    return 0.5 * std::erfc(-value / std::sqrt(2.0));
}

double foldedGaussianBinProbability(
    double truthPt, double low, double high, double sigma) {
    if (sigma == 0.0) {
        return truthPt >= low && truthPt < high ? 1.0 : 0.0;
    }
    return std::clamp(
        normalCdf((high - truthPt) / sigma) -
            normalCdf((low - truthPt) / sigma) +
            normalCdf((-low - truthPt) / sigma) -
            normalCdf((-high - truthPt) / sigma),
        0.0,
        1.0);
}

void addContribution(
    StageBins& bins,
    int truthBin,
    double truthPt,
    double weight,
    bool isolated,
    const Options& options) {
    bins[TruthStage][static_cast<std::size_t>(truthBin)] += weight;
    for (int bin = 0; bin < BinCount; ++bin) {
        const double smearedContribution =
            weight *
            foldedGaussianBinProbability(
                truthPt,
                static_cast<double>(bin) * PtBinWidth,
                static_cast<double>(bin + 1) * PtBinWidth,
                options.ptResolutionGeV);
        bins[TruthSmearedStage][static_cast<std::size_t>(bin)] +=
            smearedContribution;
        if (isolated) {
            bins[SmearedStage][static_cast<std::size_t>(bin)] +=
                smearedContribution;
        }
    }
    if (isolated) {
        bins[IsolatedStage][static_cast<std::size_t>(truthBin)] += weight;
    }
}

Estimate directIntegrate(
    const lowmc::LowCurrent& current,
    const std::vector<SheetHit>& trackHits,
    const Options& options,
    int logEBins,
    int etaBins,
    int phiBins,
    int refinementLevels,
    int subBins) {
    const auto start = Clock::now();
    const lowmc::PhotonPhaseSpace phaseSpace = angularPhaseSpace();
    lowmc::GridAdaptiveEtaPhiProposalConfig gridConfig =
        lowmc::tunedMassGridAdaptiveEtaPhiConfig();
    gridConfig.etaBins = etaBins;
    gridConfig.phiBins = phiBins;
    gridConfig.refinementLevels = refinementLevels;
    const lowmc::GridAdaptiveEtaPhiProposal grid =
        lowmc::makeGridAdaptiveEtaPhiProposal(
            current,
            phaseSpace,
            0.002,
            gridConfig);

    Estimate result;
    result.gridCells = static_cast<int>(grid.cells.size());
    const double logEnergyMin = std::log(EnergyMin);
    const double dLogEnergy =
        std::log(EnergyMax / EnergyMin) / static_cast<double>(logEBins);
    for (int energyIndex = 0; energyIndex < logEBins; ++energyIndex) {
        const double energy = std::exp(
            logEnergyMin +
            (static_cast<double>(energyIndex) + 0.5) * dLogEnergy);
        for (const lowmc::GridAdaptiveEtaPhiCell& cell : grid.cells) {
            const double dEta =
                (cell.etaMax - cell.etaMin) / static_cast<double>(subBins);
            const double dPhi =
                (cell.phiMax - cell.phiMin) / static_cast<double>(subBins);
            for (int etaIndex = 0; etaIndex < subBins; ++etaIndex) {
                const double eta =
                    cell.etaMin +
                    (static_cast<double>(etaIndex) + 0.5) * dEta;
                const double inverseCosh = 1.0 / std::cosh(eta);
                const double nz = std::tanh(eta);
                for (int phiIndex = 0; phiIndex < subBins; ++phiIndex) {
                    const double phi =
                        cell.phiMin +
                        (static_cast<double>(phiIndex) + 0.5) * dPhi;
                    const double nx = inverseCosh * std::cos(phi);
                    const double ny = inverseCosh * std::sin(phi);
                    const double pt = energy * inverseCosh;
                    const int bin = ptBin(pt);
                    if (bin >= 0) {
                        const lowmc::FourVector photon =
                            photonFromEnergyEtaPhi(energy, eta, phi);
                        const double physicalMeasure =
                            energy * energy * inverseCosh * inverseCosh *
                            dLogEnergy * dEta * dPhi;
                        const double contribution =
                            lowmc::lowPhotonDensity(current, photon) *
                            physicalMeasure;
                        addContribution(
                            result.yield,
                            bin,
                            pt,
                            contribution,
                            passesIsolation(
                                nx, ny, nz, trackHits, options),
                            options);
                    }
                    ++result.evaluations;
                }
            }
        }
    }
    result.milliseconds =
        std::chrono::duration<double, std::milli>(Clock::now() - start).count();
    return result;
}

Estimate sampleLowMC(
    const lowmc::LowCurrent& current,
    const std::vector<SheetHit>& trackHits,
    const Options& options,
    std::uint64_t seed) {
    const auto start = Clock::now();
    const lowmc::PhotonPhaseSpace phaseSpace = angularPhaseSpace();
    const lowmc::LegCenteredEtaPhiProposal proposal =
        lowmc::makeLegCenteredEtaPhiProposal(
            phaseSpace,
            current.particles(),
            lowmc::LegCenteredEtaPhiProposalConfig{
                options.legUniformFraction,
                1.0e-12,
            });
    StageBins sums{};
    StageBins sums2{};
    std::array<double, StageCount> totalSums{};
    std::array<double, StageCount> totalSums2{};
    std::array<double, StageCount> comparisonSums{};
    std::array<double, StageCount> comparisonSums2{};
    std::mt19937_64 rng(seed);
    const double logEnergyRange = std::log(EnergyMax / EnergyMin);
    for (int sampleIndex = 0; sampleIndex < options.mcSamples; ++sampleIndex) {
        const double energy =
            EnergyMin *
            std::exp(
                std::generate_canonical<double, 53>(rng) *
                logEnergyRange);
        const lowmc::LegCenteredEtaPhiProposalPoint point =
            lowmc::sampleLegCenteredEtaPhi(
                proposal,
                std::generate_canonical<double, 53>(rng),
                std::generate_canonical<double, 53>(rng),
                std::generate_canonical<double, 53>(rng));
        if (!point.insideAcceptance) {
            continue;
        }
        const double inverseCosh = 1.0 / std::cosh(point.eta);
        const double nx = inverseCosh * std::cos(point.phi);
        const double ny = inverseCosh * std::sin(point.phi);
        const double nz = std::tanh(point.eta);
        const double pt = energy * inverseCosh;
        const int bin = ptBin(pt);
        if (bin < 0) {
            continue;
        }
        const lowmc::FourVector photon =
            photonFromEnergyEtaPhi(energy, point.eta, point.phi);
        const double weight =
            lowmc::lowPhotonDensity(current, photon) *
            energy * energy * inverseCosh * inverseCosh *
            logEnergyRange / point.density;
        StageBins contribution{};
        addContribution(
            contribution,
            bin,
            pt,
            weight,
            passesIsolation(nx, ny, nz, trackHits, options),
            options);
        for (int stage = 0; stage < StageCount; ++stage) {
            double total = 0.0;
            double comparisonTotal = 0.0;
            for (int outputBin = 0; outputBin < BinCount; ++outputBin) {
                const double value =
                    contribution[static_cast<std::size_t>(stage)]
                                [static_cast<std::size_t>(outputBin)];
                sums[static_cast<std::size_t>(stage)]
                    [static_cast<std::size_t>(outputBin)] += value;
                sums2[static_cast<std::size_t>(stage)]
                     [static_cast<std::size_t>(outputBin)] += value * value;
                total += value;
                if (outputBin >= ComparisonFirstBin &&
                    outputBin < ComparisonLastBinExclusive) {
                    comparisonTotal += value;
                }
            }
            totalSums[static_cast<std::size_t>(stage)] += total;
            totalSums2[static_cast<std::size_t>(stage)] += total * total;
            comparisonSums[static_cast<std::size_t>(stage)] +=
                comparisonTotal;
            comparisonSums2[static_cast<std::size_t>(stage)] +=
                comparisonTotal * comparisonTotal;
        }
    }

    Estimate result;
    result.evaluations = options.mcSamples;
    const double sampleCount = static_cast<double>(options.mcSamples);
    for (int stage = 0; stage < StageCount; ++stage) {
        for (int bin = 0; bin < BinCount; ++bin) {
            const double mean =
                sums[static_cast<std::size_t>(stage)]
                    [static_cast<std::size_t>(bin)] /
                sampleCount;
            const double secondMoment =
                sums2[static_cast<std::size_t>(stage)]
                     [static_cast<std::size_t>(bin)] /
                sampleCount;
            result.yield[static_cast<std::size_t>(stage)]
                        [static_cast<std::size_t>(bin)] = mean;
            result.standardError[static_cast<std::size_t>(stage)]
                                [static_cast<std::size_t>(bin)] =
                std::sqrt(
                    std::max(0.0, secondMoment - mean * mean) /
                    sampleCount);
        }
        const double totalMean =
            totalSums[static_cast<std::size_t>(stage)] / sampleCount;
        const double totalSecondMoment =
            totalSums2[static_cast<std::size_t>(stage)] / sampleCount;
        result.totalStandardError[static_cast<std::size_t>(stage)] =
            std::sqrt(
                std::max(
                    0.0,
                    totalSecondMoment - totalMean * totalMean) /
                sampleCount);
        const double comparisonMean =
            comparisonSums[static_cast<std::size_t>(stage)] / sampleCount;
        const double comparisonSecondMoment =
            comparisonSums2[static_cast<std::size_t>(stage)] / sampleCount;
        result.comparisonStandardError[static_cast<std::size_t>(stage)] =
            std::sqrt(
                std::max(
                    0.0,
                    comparisonSecondMoment -
                        comparisonMean * comparisonMean) /
                sampleCount);
    }
    result.milliseconds =
        std::chrono::duration<double, std::milli>(Clock::now() - start).count();
    return result;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Options options = parseOptions(argc, argv);
        std::ofstream output(options.outputPath);
        if (!output) {
            throw std::runtime_error(
                "could not open output CSV: " + options.outputPath);
        }
        output
            << "seed,event_index,sampling_seed,attempt_index,"
               "charged_multiplicity,selected_charged,selected_incoming,"
               "selected_outgoing,signed_charge_sum,track_hits,"
               "beam_momentum_GeV,lead_distance_m,isolation_distance_m,"
               "magnetic_field_T,sheet_half_size_m,pt_resolution_GeV,"
               "bin_index,pt_min_GeV,pt_max_GeV,"
               "direct_truth_yield_event,direct_isolated_yield_event,"
               "direct_truth_smeared_yield_event,"
               "direct_isolated_smeared_yield_event,"
               "direct_coarse_truth_smeared_yield_event,"
               "direct_coarse_isolated_smeared_yield_event,"
               "lowmc_truth_yield_event,lowmc_truth_standard_error,"
               "lowmc_isolated_yield_event,"
               "lowmc_isolated_standard_error,"
               "lowmc_truth_smeared_yield_event,"
               "lowmc_truth_smeared_standard_error,"
               "lowmc_truth_smeared_total_standard_error_event,"
               "lowmc_truth_smeared_comparison_standard_error_event,"
               "lowmc_isolated_smeared_yield_event,"
               "lowmc_isolated_smeared_standard_error,"
               "lowmc_isolated_smeared_total_standard_error_event,"
               "mc_samples,direct_evaluations,direct_grid_cells,"
               "direct_milliseconds,direct_coarse_milliseconds,"
               "lowmc_mc_milliseconds\n";
        output << std::setprecision(17);

        Pythia8::Pythia pythia(LOWMC_PYTHIA8_XMLDOC, false);
        initializePythia(pythia, options);

        StageBins directTotals{};
        StageBins coarseTotals{};
        StageBins mcTotals{};
        StageBins mcVarianceTotals{};
        std::array<double, StageCount> mcTotalVariance{};
        int acceptedEvents = 0;
        int validEventIndex = 0;
        int attemptedEvents = 0;
        int totalTrackHits = 0;

        while (acceptedEvents < options.events) {
            ++attemptedEvents;
            if (attemptedEvents >
                (options.firstEvent + options.events) * 10000) {
                throw std::runtime_error(
                    "too many generated events failed the WA102 selection");
            }
            if (!pythia.next()) {
                continue;
            }
            const EventSelection selection =
                selectEvent(pythia.event, options);
            if (!selection.accepted) {
                continue;
            }
            const lowmc::AdapterParticleSelectionResult particles =
                lowmc::pythia8::selectChargedParticles(pythia.event);
            if (!particles.diagnostic.selectedChargedCount() ||
                std::abs(particles.diagnostic.signedChargeSum) > 1.0e-12) {
                continue;
            }
            const int eventIndex = validEventIndex++;
            if (eventIndex < options.firstEvent) {
                continue;
            }

            const lowmc::LowCurrent current(particles.selectedParticles);
            const Estimate direct = directIntegrate(
                current,
                selection.trackHits,
                options,
                options.directLogEBins,
                options.directEtaBins,
                options.directPhiBins,
                options.directRefinementLevels,
                options.directSubBins);
            const Estimate coarse = directIntegrate(
                current,
                selection.trackHits,
                options,
                std::max(1, options.directLogEBins / 2),
                std::max(2, options.directEtaBins / 2),
                std::max(2, options.directPhiBins / 2),
                std::max(0, options.directRefinementLevels - 1),
                options.directSubBins);
            const std::uint64_t eventSamplingSeed =
                (options.samplingSeed != 0
                     ? options.samplingSeed
                     : options.seed) +
                104729ULL * static_cast<std::uint64_t>(eventIndex + 1);
            const Estimate mc = sampleLowMC(
                current,
                selection.trackHits,
                options,
                eventSamplingSeed);

            for (int stage = 0; stage < StageCount; ++stage) {
                mcTotalVariance[static_cast<std::size_t>(stage)] +=
                    mc.totalStandardError[static_cast<std::size_t>(stage)] *
                    mc.totalStandardError[static_cast<std::size_t>(stage)];
                for (int bin = 0; bin < BinCount; ++bin) {
                    const std::size_t stageIndex =
                        static_cast<std::size_t>(stage);
                    const std::size_t binIndex =
                        static_cast<std::size_t>(bin);
                    directTotals[stageIndex][binIndex] +=
                        direct.yield[stageIndex][binIndex];
                    coarseTotals[stageIndex][binIndex] +=
                        coarse.yield[stageIndex][binIndex];
                    mcTotals[stageIndex][binIndex] +=
                        mc.yield[stageIndex][binIndex];
                    mcVarianceTotals[stageIndex][binIndex] +=
                        mc.standardError[stageIndex][binIndex] *
                        mc.standardError[stageIndex][binIndex];
                }
            }
            totalTrackHits += static_cast<int>(selection.trackHits.size());

            for (int bin = 0; bin < BinCount; ++bin) {
                const std::size_t index = static_cast<std::size_t>(bin);
                output
                    << options.seed << ','
                    << eventIndex << ','
                    << eventSamplingSeed << ','
                    << attemptedEvents << ','
                    << selection.chargedMultiplicity << ','
                    << particles.diagnostic.selectedChargedCount() << ','
                    << particles.diagnostic.selectedIncomingCharged << ','
                    << particles.diagnostic.selectedOutgoingCharged << ','
                    << particles.diagnostic.signedChargeSum << ','
                    << selection.trackHits.size() << ','
                    << BeamMomentumGeV << ','
                    << options.leadDistanceM << ','
                    << options.isolationDistanceM << ','
                    << options.magneticFieldT << ','
                    << options.sheetHalfSizeM << ','
                    << options.ptResolutionGeV << ','
                    << bin << ','
                    << static_cast<double>(bin) * PtBinWidth << ','
                    << static_cast<double>(bin + 1) * PtBinWidth << ','
                    << direct.yield[TruthStage][index] << ','
                    << direct.yield[IsolatedStage][index] << ','
                    << direct.yield[TruthSmearedStage][index] << ','
                    << direct.yield[SmearedStage][index] << ','
                    << coarse.yield[TruthSmearedStage][index] << ','
                    << coarse.yield[SmearedStage][index] << ','
                    << mc.yield[TruthStage][index] << ','
                    << mc.standardError[TruthStage][index] << ','
                    << mc.yield[IsolatedStage][index] << ','
                    << mc.standardError[IsolatedStage][index] << ','
                    << mc.yield[TruthSmearedStage][index] << ','
                    << mc.standardError[TruthSmearedStage][index] << ','
                    << mc.totalStandardError[TruthSmearedStage] << ','
                    << mc.comparisonStandardError[TruthSmearedStage] << ','
                    << mc.yield[SmearedStage][index] << ','
                    << mc.standardError[SmearedStage][index] << ','
                    << mc.totalStandardError[SmearedStage] << ','
                    << options.mcSamples << ','
                    << direct.evaluations << ','
                    << direct.gridCells << ','
                    << direct.milliseconds << ','
                    << coarse.milliseconds << ','
                    << mc.milliseconds << '\n';
            }
            ++acceptedEvents;
        }
        if (!output) {
            throw std::runtime_error("failed while writing WA102 CSV");
        }

        auto total = [](const StageBins& values, int stage) {
            double result = 0.0;
            for (double value :
                 values[static_cast<std::size_t>(stage)]) {
                result += value;
            }
            return result;
        };
        for (int bin = 0; bin < BinCount; ++bin) {
            const std::size_t index = static_cast<std::size_t>(bin);
            std::cout
                << "bin_index=" << bin
                << " pt_min_GeV=" << static_cast<double>(bin) * PtBinWidth
                << " pt_max_GeV="
                << static_cast<double>(bin + 1) * PtBinWidth
                << " direct_isolated_smeared_per_event="
                << directTotals[SmearedStage][index] /
                       static_cast<double>(acceptedEvents)
                << " direct_coarse_isolated_smeared_per_event="
                << coarseTotals[SmearedStage][index] /
                       static_cast<double>(acceptedEvents)
                << " lowmc_isolated_smeared_per_event="
                << mcTotals[SmearedStage][index] /
                       static_cast<double>(acceptedEvents)
                << " lowmc_standard_error="
                << std::sqrt(mcVarianceTotals[SmearedStage][index]) /
                       static_cast<double>(acceptedEvents)
                << '\n';
        }

        const double directTruth =
            total(directTotals, TruthStage) /
            static_cast<double>(acceptedEvents);
        const double directIsolated =
            total(directTotals, IsolatedStage) /
            static_cast<double>(acceptedEvents);
        const double directTruthSmeared =
            total(directTotals, TruthSmearedStage) /
            static_cast<double>(acceptedEvents);
        const double directSmeared =
            total(directTotals, SmearedStage) /
            static_cast<double>(acceptedEvents);
        const double coarseSmeared =
            total(coarseTotals, SmearedStage) /
            static_cast<double>(acceptedEvents);
        const double mcSmeared =
            total(mcTotals, SmearedStage) /
            static_cast<double>(acceptedEvents);
        std::cout
            << "aggregate_status=ok"
            << " accepted_events=" << acceptedEvents
            << " attempted_events=" << attemptedEvents
            << " mean_track_hits="
            << static_cast<double>(totalTrackHits) /
                   static_cast<double>(acceptedEvents)
            << " lead_distance_m=" << options.leadDistanceM
            << " isolation_distance_m=" << options.isolationDistanceM
            << " magnetic_field_T=" << options.magneticFieldT
            << " pt_resolution_GeV=" << options.ptResolutionGeV
            << " direct_truth_per_event=" << directTruth
            << " direct_isolated_per_event=" << directIsolated
            << " isolation_efficiency="
            << (directTruth != 0.0 ? directIsolated / directTruth : 0.0)
            << " direct_truth_smeared_per_event=" << directTruthSmeared
            << " direct_isolated_smeared_per_event=" << directSmeared
            << " smear_display_efficiency="
            << (directIsolated != 0.0 ? directSmeared / directIsolated : 0.0)
            << " direct_coarse_relative_difference="
            << (directSmeared != 0.0
                    ? std::abs(directSmeared - coarseSmeared) /
                          std::abs(directSmeared)
                    : 0.0)
            << " lowmc_isolated_smeared_per_event=" << mcSmeared
            << " lowmc_standard_error="
            << std::sqrt(mcTotalVariance[SmearedStage]) /
                   static_cast<double>(acceptedEvents)
            << " lowmc_over_direct="
            << (directSmeared != 0.0 ? mcSmeared / directSmeared : 0.0)
            << '\n';
    } catch (const std::exception& error) {
        std::cerr << "FAILED: " << error.what() << '\n';
        printUsage(std::cerr);
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
