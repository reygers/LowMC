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
constexpr int JointPtBinCount = 50;
constexpr int ThetaBinCount = 25;
constexpr int JointBinCount = JointPtBinCount * ThetaBinCount;
constexpr double EnergyMin = 0.2;
constexpr double EnergyMax = 1.0;
constexpr double PtBinWidth = 0.008;
constexpr double PtMax = BinCount * PtBinWidth;
constexpr double JointPtMax = JointPtBinCount * PtBinWidth;
constexpr double JetAngleMax = 0.4;
constexpr double ThetaBinWidth = JetAngleMax / ThetaBinCount;

struct Options {
    std::string outputPath = "pythia8_delphi_hadronic_benchmark.csv";
    std::string jointOutputPath;
    int events = 2;
    int firstEvent = 0;
    int mcSamples = 200000;
    int directLogEBins = 24;
    int directEtaBins = 24;
    int directPhiBins = 32;
    int directRefinementLevels = 2;
    int directSubBins = 2;
    std::uint64_t seed = 12345;
    std::uint64_t samplingSeed = 0;
    double legUniformFraction = 0.1;
    bool electronVeto = true;
};

struct JetAxis {
    int clusterIndex = -1;
    double nx = 0.0;
    double ny = 0.0;
    double nz = 0.0;
};

struct BinEstimate {
    std::array<double, BinCount> yield{};
    std::array<double, BinCount> standardError{};
    std::array<double, JointBinCount> jointYield{};
    std::array<double, JointBinCount> jointStandardError{};
    std::array<double, ThetaBinCount> thetaYield{};
    std::array<double, ThetaBinCount> thetaStandardError{};
    double totalStandardError = 0.0;
    int evaluations = 0;
    int gridCells = 0;
    double milliseconds = 0.0;
};

struct EventSelection {
    bool accepted = false;
    int chargedMultiplicity = 0;
    double visibleEnergy = 0.0;
    int totalJets = 0;
    std::vector<JetAxis> acceptedJets;
};

void printUsage(std::ostream& output) {
    output
        << "Usage: lowmc_pythia8_delphi_hadronic_benchmark [output.csv] [options]\n\n"
        << "Compare direct bin integration and weighted LowMC sampling for the\n"
        << "DELPHI Figure 4d inner-bremsstrahlung phase space on identical\n"
        << "Pythia hadronic-Z events. A joint pT-theta CSV also supports the\n"
        << "Figure 4b angular projection through pT = 0.4 GeV. This diagnostic\n"
        << "does not emit photons.\n\n"
        << "Primary fixed cuts:\n"
        << "  0.2 < E_gamma < 1 GeV, theta_beam in [20,160] deg,\n"
        << "  theta_gamma_jet < 0.4 rad, pT_gamma_jet < 80 MeV,\n"
        << "  ten 8 MeV pT bins, normalized per accepted jet.\n"
        << "  The joint diagnostic retains 50 bins through pT = 0.4 GeV.\n\n"
        << "Options:\n"
        << "  --events N --first-event N\n"
        << "  --joint-output PATH\n"
        << "  --mc-samples N\n"
        << "  --direct-loge-bins N\n"
        << "  --direct-eta-bins N --direct-phi-bins N\n"
        << "  --direct-refinement-levels N --direct-subbins N\n"
        << "  --seed N --sampling-seed N\n"
        << "  --leg-uniform-fraction F\n"
        << "  --no-electron-veto\n"
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

std::uint64_t parseUInt64(const std::string& text, const std::string& optionName) {
    std::size_t parsed = 0;
    const unsigned long long value = std::stoull(text, &parsed);
    if (parsed != text.size()) {
        throw std::invalid_argument(optionName + " expects an unsigned integer");
    }
    return static_cast<std::uint64_t>(value);
}

double parseFraction(const std::string& text, const std::string& optionName) {
    std::size_t parsed = 0;
    const double value = std::stod(text, &parsed);
    if (parsed != text.size() || !std::isfinite(value) ||
        !(value > 0.0 && value < 1.0)) {
        throw std::invalid_argument(optionName + " expects 0 < value < 1");
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
            options.events = parseInteger(
                readOption(index, argc, argv, argument), argument, false);
        } else if (argument == "--first-event") {
            options.firstEvent = parseInteger(
                readOption(index, argc, argv, argument), argument, true);
        } else if (argument == "--joint-output") {
            options.jointOutputPath =
                readOption(index, argc, argv, argument);
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
            options.seed = parseUInt64(
                readOption(index, argc, argv, argument), argument);
        } else if (argument == "--sampling-seed") {
            options.samplingSeed = parseUInt64(
                readOption(index, argc, argv, argument), argument);
        } else if (argument == "--leg-uniform-fraction") {
            options.legUniformFraction = parseFraction(
                readOption(index, argc, argv, argument), argument);
        } else if (argument == "--no-electron-veto") {
            options.electronVeto = false;
        } else if (!argument.empty() && argument.front() == '-') {
            throw std::invalid_argument("unknown argument: " + argument);
        } else if (outputSpecified) {
            throw std::invalid_argument("only one output CSV path can be specified");
        } else {
            options.outputPath = argument;
            outputSpecified = true;
        }
    }
    if (options.jointOutputPath.empty()) {
        const std::string suffix = ".csv";
        if (options.outputPath.size() >= suffix.size() &&
            options.outputPath.compare(
                options.outputPath.size() - suffix.size(),
                suffix.size(),
                suffix) == 0) {
            options.jointOutputPath =
                options.outputPath.substr(
                    0, options.outputPath.size() - suffix.size()) +
                "_joint.csv";
        } else {
            options.jointOutputPath = options.outputPath + "_joint.csv";
        }
    }
    return options;
}

double pi() {
    return std::acos(-1.0);
}

double thetaToEta(double theta) {
    return -std::log(std::tan(0.5 * theta));
}

lowmc::PhotonPhaseSpace angularPhaseSpace() {
    const double thetaMin = 20.0 * pi() / 180.0;
    const double thetaMax = 160.0 * pi() / 180.0;
    return lowmc::PhotonPhaseSpace{
        0.001,
        1.0,
        thetaToEta(thetaMax),
        thetaToEta(thetaMin),
        0.0,
        2.0 * pi(),
    };
}

void initializePythia(Pythia8::Pythia& pythia, const Options& options) {
    pythia.readString("Print:quiet = on");
    pythia.readString("Random:setSeed = on");
    pythia.readString("Random:seed = " + std::to_string(options.seed));
    pythia.readString("PDF:lepton = off");
    pythia.readString("WeakSingleBoson:ffbar2gmZ = on");
    pythia.readString("23:onMode = off");
    pythia.readString("23:onIfAny = 1 2 3 4 5");
    pythia.readString("Beams:idA = 11");
    pythia.readString("Beams:idB = -11");
    pythia.readString(
        "Beams:eCM = " + std::to_string(pythia.particleData.m0(23)));
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

bool jetContainsElectron(
    const Pythia8::Event& event,
    const Pythia8::ClusterJet& clustering,
    int jetIndex) {
    for (int index = 0; index < event.size(); ++index) {
        if (event[index].isFinal() && std::abs(event[index].id()) == 11 &&
            clustering.jetAssignment(index) == jetIndex) {
            return true;
        }
    }
    return false;
}

EventSelection selectEvent(
    const Pythia8::Event& event,
    Pythia8::Thrust& thrust,
    Pythia8::ClusterJet& clustering,
    bool electronVeto) {
    EventSelection result;
    for (int index = 0; index < event.size(); ++index) {
        const Pythia8::Particle& particle = event[index];
        if (!particle.isFinal()) {
            continue;
        }
        result.chargedMultiplicity += particle.isCharged() ? 1 : 0;
        result.visibleEnergy += particle.isVisible() ? particle.e() : 0.0;
    }
    if (result.chargedMultiplicity < 5 ||
        result.visibleEnergy < 0.2 * event[0].e()) {
        return result;
    }
    if (!thrust.analyze(event)) {
        return result;
    }
    const double thrustAbsCosine = std::abs(thrust.eventAxis(1).pz());
    if (thrustAbsCosine > std::cos(30.0 * pi() / 180.0)) {
        return result;
    }
    if (!clustering.analyze(event, 0.0, 3.0) || clustering.size() < 2) {
        return result;
    }
    result.totalJets = clustering.size();
    for (int jetIndex = 0; jetIndex < clustering.size(); ++jetIndex) {
        const Pythia8::Vec4 momentum = clustering.p(jetIndex);
        const double magnitude = momentum.pAbs();
        if (!(magnitude >= 5.0)) {
            continue;
        }
        const double cosine = momentum.pz() / magnitude;
        const double theta = std::acos(std::clamp(cosine, -1.0, 1.0));
        if (theta < 30.0 * pi() / 180.0 ||
            theta > 150.0 * pi() / 180.0) {
            continue;
        }
        if (electronVeto &&
            jetContainsElectron(event, clustering, jetIndex)) {
            continue;
        }
        result.acceptedJets.push_back(JetAxis{
            jetIndex,
            momentum.px() / magnitude,
            momentum.py() / magnitude,
            momentum.pz() / magnitude,
        });
    }
    result.accepted = !result.acceptedJets.empty();
    return result;
}

double closestJetCosine(
    double nx,
    double ny,
    double nz,
    const std::vector<JetAxis>& jets) {
    double result = -1.0;
    for (const JetAxis& jet : jets) {
        result = std::max(
            result,
            nx * jet.nx + ny * jet.ny + nz * jet.nz);
    }
    return std::clamp(result, -1.0, 1.0);
}

struct PhotonBins {
    int primaryPt = -1;
    int jointPt = -1;
    int theta = -1;
    int joint = -1;
};

PhotonBins photonBins(double energy, double closestCosine) {
    PhotonBins result;
    if (closestCosine <= std::cos(JetAngleMax)) {
        return result;
    }
    const double theta =
        std::acos(std::clamp(closestCosine, -1.0, 1.0));
    const double sine =
        std::sqrt(std::max(0.0, 1.0 - closestCosine * closestCosine));
    const double pt = energy * sine;
    if (!(pt >= 0.0 && pt < JointPtMax) ||
        !(theta >= 0.0 && theta < JetAngleMax)) {
        return result;
    }
    result.jointPt = std::min(
        JointPtBinCount - 1,
        static_cast<int>(pt / PtBinWidth));
    result.theta = std::min(
        ThetaBinCount - 1,
        static_cast<int>(theta / ThetaBinWidth));
    result.joint = result.jointPt * ThetaBinCount + result.theta;
    if (pt < PtMax) {
        result.primaryPt = std::min(
            BinCount - 1,
            static_cast<int>(pt / PtBinWidth));
    }
    return result;
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

BinEstimate directIntegrate(
    const lowmc::LowCurrent& current,
    const std::vector<JetAxis>& jets,
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
            0.02,
            gridConfig);

    BinEstimate result;
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
                    const PhotonBins bins =
                        photonBins(
                            energy,
                            closestJetCosine(nx, ny, nz, jets));
                    if (bins.joint >= 0) {
                        const lowmc::FourVector photon =
                            photonFromEnergyEtaPhi(energy, eta, phi);
                        const double physicalMeasure =
                            energy * energy * inverseCosh * inverseCosh *
                            dLogEnergy * dEta * dPhi;
                        const double contribution =
                            lowmc::lowPhotonDensity(current, photon) *
                            physicalMeasure;
                        result.jointYield[
                            static_cast<std::size_t>(bins.joint)] +=
                            contribution;
                        result.thetaYield[
                            static_cast<std::size_t>(bins.theta)] +=
                            contribution;
                        if (bins.primaryPt >= 0) {
                            result.yield[
                                static_cast<std::size_t>(bins.primaryPt)] +=
                                contribution;
                        }
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

BinEstimate sampleLowMC(
    const lowmc::LowCurrent& current,
    const std::vector<JetAxis>& jets,
    int samples,
    std::uint64_t seed,
    double uniformFraction) {
    const auto start = Clock::now();
    const lowmc::PhotonPhaseSpace phaseSpace = angularPhaseSpace();
    const lowmc::LegCenteredEtaPhiProposal proposal =
        lowmc::makeLegCenteredEtaPhiProposal(
            phaseSpace,
            current.particles(),
            lowmc::LegCenteredEtaPhiProposalConfig{
                uniformFraction,
                1.0e-12,
            });
    std::array<double, BinCount> sums{};
    std::array<double, BinCount> sums2{};
    std::array<double, JointBinCount> jointSums{};
    std::array<double, JointBinCount> jointSums2{};
    std::array<double, ThetaBinCount> thetaSums{};
    std::array<double, ThetaBinCount> thetaSums2{};
    double totalSum = 0.0;
    double totalSum2 = 0.0;
    std::mt19937_64 rng(seed);
    const double logEnergyRange = std::log(EnergyMax / EnergyMin);
    for (int sampleIndex = 0; sampleIndex < samples; ++sampleIndex) {
        const double unitEnergy =
            std::generate_canonical<double, 53>(rng);
        const double energy =
            EnergyMin * std::exp(unitEnergy * logEnergyRange);
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
        const PhotonBins bins =
            photonBins(
                energy,
                closestJetCosine(nx, ny, nz, jets));
        if (bins.joint < 0) {
            continue;
        }
        const lowmc::FourVector photon =
            photonFromEnergyEtaPhi(energy, point.eta, point.phi);
        const double weight =
            lowmc::lowPhotonDensity(current, photon) *
            energy * energy * inverseCosh * inverseCosh *
            logEnergyRange / point.density;
        jointSums[static_cast<std::size_t>(bins.joint)] += weight;
        jointSums2[static_cast<std::size_t>(bins.joint)] += weight * weight;
        thetaSums[static_cast<std::size_t>(bins.theta)] += weight;
        thetaSums2[static_cast<std::size_t>(bins.theta)] += weight * weight;
        if (bins.primaryPt >= 0) {
            sums[static_cast<std::size_t>(bins.primaryPt)] += weight;
            sums2[static_cast<std::size_t>(bins.primaryPt)] +=
                weight * weight;
            totalSum += weight;
            totalSum2 += weight * weight;
        }
    }

    BinEstimate result;
    result.evaluations = samples;
    for (int bin = 0; bin < BinCount; ++bin) {
        const double mean =
            sums[static_cast<std::size_t>(bin)] /
            static_cast<double>(samples);
        const double secondMoment =
            sums2[static_cast<std::size_t>(bin)] /
            static_cast<double>(samples);
        const double variance =
            std::max(0.0, secondMoment - mean * mean);
        result.yield[static_cast<std::size_t>(bin)] = mean;
        result.standardError[static_cast<std::size_t>(bin)] =
            std::sqrt(variance / static_cast<double>(samples));
    }
    for (int bin = 0; bin < JointBinCount; ++bin) {
        const double mean =
            jointSums[static_cast<std::size_t>(bin)] /
            static_cast<double>(samples);
        const double secondMoment =
            jointSums2[static_cast<std::size_t>(bin)] /
            static_cast<double>(samples);
        const double variance =
            std::max(0.0, secondMoment - mean * mean);
        result.jointYield[static_cast<std::size_t>(bin)] = mean;
        result.jointStandardError[static_cast<std::size_t>(bin)] =
            std::sqrt(variance / static_cast<double>(samples));
    }
    for (int bin = 0; bin < ThetaBinCount; ++bin) {
        const double mean =
            thetaSums[static_cast<std::size_t>(bin)] /
            static_cast<double>(samples);
        const double secondMoment =
            thetaSums2[static_cast<std::size_t>(bin)] /
            static_cast<double>(samples);
        const double variance =
            std::max(0.0, secondMoment - mean * mean);
        result.thetaYield[static_cast<std::size_t>(bin)] = mean;
        result.thetaStandardError[static_cast<std::size_t>(bin)] =
            std::sqrt(variance / static_cast<double>(samples));
    }
    const double totalMean = totalSum / static_cast<double>(samples);
    const double totalSecondMoment =
        totalSum2 / static_cast<double>(samples);
    result.totalStandardError = std::sqrt(
        std::max(0.0, totalSecondMoment - totalMean * totalMean) /
        static_cast<double>(samples));
    result.milliseconds =
        std::chrono::duration<double, std::milli>(Clock::now() - start).count();
    return result;
}

void validateJointProjections(
    const BinEstimate& estimate,
    const std::string& label) {
    constexpr double tolerance = 1.0e-12;
    for (int ptBinIndex = 0; ptBinIndex < BinCount; ++ptBinIndex) {
        double projection = 0.0;
        for (int thetaBin = 0; thetaBin < ThetaBinCount; ++thetaBin) {
            projection += estimate.jointYield[
                static_cast<std::size_t>(
                    ptBinIndex * ThetaBinCount + thetaBin)];
        }
        if (std::abs(
                projection -
                estimate.yield[static_cast<std::size_t>(ptBinIndex)]) >
            tolerance *
                std::max(1.0, std::abs(projection))) {
            throw std::runtime_error(
                label + " joint-to-pT projection is inconsistent");
        }
    }
    for (int thetaBin = 0; thetaBin < ThetaBinCount; ++thetaBin) {
        double projection = 0.0;
        for (int ptBinIndex = 0;
             ptBinIndex < JointPtBinCount;
             ++ptBinIndex) {
            projection += estimate.jointYield[
                static_cast<std::size_t>(
                    ptBinIndex * ThetaBinCount + thetaBin)];
        }
        if (std::abs(
                projection -
                estimate.thetaYield[static_cast<std::size_t>(thetaBin)]) >
            tolerance *
                std::max(1.0, std::abs(projection))) {
            throw std::runtime_error(
                label + " joint-to-theta projection is inconsistent");
        }
    }
}

int hardFlavor(const Pythia8::Event& event) {
    for (int index = 0; index < event.size(); ++index) {
        const int idAbs = std::abs(event[index].id());
        if (idAbs >= 1 && idAbs <= 5 &&
            std::abs(event[index].status()) == 23) {
            return idAbs;
        }
    }
    return 0;
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
        std::ofstream jointOutput(options.jointOutputPath);
        if (!jointOutput) {
            throw std::runtime_error(
                "could not open joint output CSV: " +
                options.jointOutputPath);
        }
        output
            << "seed,event_index,sampling_seed,hard_flavor,charged_multiplicity,"
               "visible_energy,total_jets,accepted_jets,selected_charged,"
               "selected_outgoing,signed_charge_sum,electron_veto,"
               "generator_qed_radiation,theta_gamma_jet_max,bin_index,"
               "pt_min_GeV,pt_max_GeV,direct_yield_event,"
               "direct_coarse_yield_event,direct_relative_difference,"
               "lowmc_mc_yield_event,lowmc_mc_standard_error,"
               "lowmc_mc_total_standard_error_event,mc_samples,"
               "direct_evaluations,direct_grid_cells,direct_milliseconds,"
               "direct_coarse_milliseconds,lowmc_mc_milliseconds\n";
        output << std::setprecision(17);
        jointOutput
            << "seed,event_index,sampling_seed,accepted_jets,"
               "pt_bin_index,pt_min_GeV,pt_max_GeV,"
               "theta_bin_index,theta_min_rad,theta_max_rad,"
               "direct_yield_event,direct_coarse_yield_event,"
               "lowmc_mc_yield_event,lowmc_mc_joint_standard_error,"
               "lowmc_mc_theta_standard_error_event,mc_samples\n";
        jointOutput << std::setprecision(17);

        Pythia8::Pythia pythia(LOWMC_PYTHIA8_XMLDOC, false);
        initializePythia(pythia, options);
        Pythia8::Thrust thrust;
        Pythia8::ClusterJet clustering("Lund");

        std::array<double, BinCount> directTotal{};
        std::array<double, BinCount> coarseTotal{};
        std::array<double, BinCount> mcTotal{};
        std::array<double, BinCount> mcVarianceTotal{};
        std::array<double, ThetaBinCount> directThetaTotal{};
        std::array<double, ThetaBinCount> coarseThetaTotal{};
        std::array<double, ThetaBinCount> mcThetaTotal{};
        std::array<double, ThetaBinCount> mcThetaVarianceTotal{};
        double mcTotalVariance = 0.0;
        int acceptedEvents = 0;
        int validEventIndex = 0;
        int attemptedEvents = 0;
        int totalAcceptedJets = 0;

        while (acceptedEvents < options.events) {
            ++attemptedEvents;
            if (attemptedEvents > (options.firstEvent + options.events) * 100) {
                throw std::runtime_error(
                    "too many generated events failed the DELPHI selection");
            }
            if (!pythia.next()) {
                continue;
            }
            const EventSelection eventSelection =
                selectEvent(
                    pythia.event,
                    thrust,
                    clustering,
                    options.electronVeto);
            if (!eventSelection.accepted) {
                continue;
            }
            const lowmc::AdapterParticleSelectionResult particleSelection =
                lowmc::pythia8::selectChargedParticles(pythia.event);
            if (!particleSelection.diagnostic.selectedChargedCount() ||
                std::abs(particleSelection.diagnostic.signedChargeSum) >
                    1.0e-12) {
                continue;
            }
            const int eventIndex = validEventIndex++;
            if (eventIndex < options.firstEvent) {
                continue;
            }

            const lowmc::LowCurrent current(
                particleSelection.selectedParticles);
            const BinEstimate direct = directIntegrate(
                current,
                eventSelection.acceptedJets,
                options.directLogEBins,
                options.directEtaBins,
                options.directPhiBins,
                options.directRefinementLevels,
                options.directSubBins);
            const BinEstimate coarse = directIntegrate(
                current,
                eventSelection.acceptedJets,
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
            const BinEstimate mc = sampleLowMC(
                current,
                eventSelection.acceptedJets,
                options.mcSamples,
                eventSamplingSeed,
                options.legUniformFraction);
            validateJointProjections(direct, "direct");
            validateJointProjections(coarse, "direct coarse");
            validateJointProjections(mc, "LowMC");
            mcTotalVariance +=
                mc.totalStandardError * mc.totalStandardError;

            totalAcceptedJets +=
                static_cast<int>(eventSelection.acceptedJets.size());
            for (int bin = 0; bin < BinCount; ++bin) {
                const std::size_t index = static_cast<std::size_t>(bin);
                directTotal[index] += direct.yield[index];
                coarseTotal[index] += coarse.yield[index];
                mcTotal[index] += mc.yield[index];
                mcVarianceTotal[index] +=
                    mc.standardError[index] * mc.standardError[index];
                const double relativeDifference =
                    direct.yield[index] != 0.0
                        ? std::abs(
                              direct.yield[index] -
                              coarse.yield[index]) /
                              std::abs(direct.yield[index])
                        : 0.0;
                output
                    << options.seed << ','
                    << eventIndex << ','
                    << eventSamplingSeed << ','
                    << hardFlavor(pythia.event) << ','
                    << eventSelection.chargedMultiplicity << ','
                    << eventSelection.visibleEnergy << ','
                    << eventSelection.totalJets << ','
                    << eventSelection.acceptedJets.size() << ','
                    << particleSelection.diagnostic.selectedChargedCount()
                    << ','
                    << particleSelection.diagnostic.selectedOutgoingCharged
                    << ','
                    << particleSelection.diagnostic.signedChargeSum << ','
                    << (options.electronVeto ? 1 : 0) << ','
                    << 0 << ','
                    << JetAngleMax << ','
                    << bin << ','
                    << static_cast<double>(bin) * PtBinWidth << ','
                    << static_cast<double>(bin + 1) * PtBinWidth << ','
                    << direct.yield[index] << ','
                    << coarse.yield[index] << ','
                    << relativeDifference << ','
                    << mc.yield[index] << ','
                    << mc.standardError[index] << ','
                    << mc.totalStandardError << ','
                    << options.mcSamples << ','
                    << direct.evaluations << ','
                    << direct.gridCells << ','
                    << direct.milliseconds << ','
                    << coarse.milliseconds << ','
                    << mc.milliseconds << '\n';
            }
            for (int thetaBin = 0; thetaBin < ThetaBinCount; ++thetaBin) {
                const std::size_t thetaIndex =
                    static_cast<std::size_t>(thetaBin);
                directThetaTotal[thetaIndex] +=
                    direct.thetaYield[thetaIndex];
                coarseThetaTotal[thetaIndex] +=
                    coarse.thetaYield[thetaIndex];
                mcThetaTotal[thetaIndex] += mc.thetaYield[thetaIndex];
                mcThetaVarianceTotal[thetaIndex] +=
                    mc.thetaStandardError[thetaIndex] *
                    mc.thetaStandardError[thetaIndex];
                for (int ptBinIndex = 0;
                     ptBinIndex < JointPtBinCount;
                     ++ptBinIndex) {
                    const int jointBin =
                        ptBinIndex * ThetaBinCount + thetaBin;
                    const std::size_t jointIndex =
                        static_cast<std::size_t>(jointBin);
                    jointOutput
                        << options.seed << ','
                        << eventIndex << ','
                        << eventSamplingSeed << ','
                        << eventSelection.acceptedJets.size() << ','
                        << ptBinIndex << ','
                        << static_cast<double>(ptBinIndex) * PtBinWidth
                        << ','
                        << static_cast<double>(ptBinIndex + 1) *
                               PtBinWidth
                        << ','
                        << thetaBin << ','
                        << static_cast<double>(thetaBin) * ThetaBinWidth
                        << ','
                        << static_cast<double>(thetaBin + 1) *
                               ThetaBinWidth
                        << ','
                        << direct.jointYield[jointIndex] << ','
                        << coarse.jointYield[jointIndex] << ','
                        << mc.jointYield[jointIndex] << ','
                        << mc.jointStandardError[jointIndex] << ','
                        << mc.thetaStandardError[thetaIndex] << ','
                        << options.mcSamples << '\n';
                }
            }
            ++acceptedEvents;
        }
        if (!output || !jointOutput) {
            throw std::runtime_error("failed while writing benchmark CSV");
        }

        double directIntegral = 0.0;
        double coarseIntegral = 0.0;
        double mcIntegral = 0.0;
        for (int bin = 0; bin < BinCount; ++bin) {
            const std::size_t index = static_cast<std::size_t>(bin);
            directIntegral += directTotal[index];
            coarseIntegral += coarseTotal[index];
            mcIntegral += mcTotal[index];
            std::cout
                << "bin_index=" << bin
                << " pt_min_GeV=" << static_cast<double>(bin) * PtBinWidth
                << " pt_max_GeV=" << static_cast<double>(bin + 1) * PtBinWidth
                << " direct_per_jet="
                << directTotal[index] /
                       static_cast<double>(totalAcceptedJets)
                << " direct_coarse_per_jet="
                << coarseTotal[index] /
                       static_cast<double>(totalAcceptedJets)
                << " lowmc_mc_per_jet="
                << mcTotal[index] /
                       static_cast<double>(totalAcceptedJets)
                << " lowmc_mc_standard_error="
                << std::sqrt(mcVarianceTotal[index]) /
                       static_cast<double>(totalAcceptedJets)
                << '\n';
        }
        for (int bin = 0; bin < ThetaBinCount; ++bin) {
            const std::size_t index = static_cast<std::size_t>(bin);
            std::cout
                << "theta_bin_index=" << bin
                << " theta_min_rad="
                << static_cast<double>(bin) * ThetaBinWidth
                << " theta_max_rad="
                << static_cast<double>(bin + 1) * ThetaBinWidth
                << " direct_per_jet="
                << directThetaTotal[index] /
                       static_cast<double>(totalAcceptedJets)
                << " direct_coarse_per_jet="
                << coarseThetaTotal[index] /
                       static_cast<double>(totalAcceptedJets)
                << " lowmc_mc_per_jet="
                << mcThetaTotal[index] /
                       static_cast<double>(totalAcceptedJets)
                << " lowmc_mc_standard_error="
                << std::sqrt(mcThetaVarianceTotal[index]) /
                       static_cast<double>(totalAcceptedJets)
                << '\n';
        }
        const double directPerJet =
            directIntegral / static_cast<double>(totalAcceptedJets);
        const double coarsePerJet =
            coarseIntegral / static_cast<double>(totalAcceptedJets);
        const double mcPerJet =
            mcIntegral / static_cast<double>(totalAcceptedJets);
        const double mcStandardError =
            std::sqrt(mcTotalVariance) /
            static_cast<double>(totalAcceptedJets);
        std::cout
            << "aggregate_status=ok"
            << " accepted_events=" << acceptedEvents
            << " accepted_jets=" << totalAcceptedJets
            << " theta_gamma_jet_max=0.4"
            << " direct_per_jet=" << directPerJet
            << " direct_coarse_per_jet=" << coarsePerJet
            << " direct_relative_difference="
            << (directPerJet != 0.0
                    ? std::abs(directPerJet - coarsePerJet) /
                          std::abs(directPerJet)
                    : 0.0)
            << " lowmc_mc_per_jet=" << mcPerJet
            << " lowmc_mc_standard_error=" << mcStandardError
            << " lowmc_mc_over_direct="
            << (directPerJet != 0.0 ? mcPerJet / directPerJet : 0.0)
            << '\n';
    } catch (const std::exception& error) {
        std::cerr << "FAILED: " << error.what() << '\n';
        printUsage(std::cerr);
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
