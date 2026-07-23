#include "ParticleInput.h"

#include "lowmc/LowMC.h"

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct Options {
    int trials = 1000;
    int calibrationTrials = 1000;
    std::uint64_t seed = 12345;
    std::string topology = "two-pion";
    std::string particlesPath;
    double momentum = 100.0;
    double particleEta = 0.0;
    double particlePhi = 0.5 * std::acos(-1.0);
    std::optional<double> secondMomentum;
    std::optional<double> secondParticleEta;
    std::optional<double> secondParticlePhi;
    double safetyFactor = 3.0;
    int probeBins = 11;
    lowmc::ShapeUnweightingProbeRule probeRule = lowmc::ShapeUnweightingProbeRule::Interior;
    lowmc::ShapeUnweightingProbePolicy probePolicy = lowmc::ShapeUnweightingProbePolicy::Fixed;
    int boostedProbeBins = 21;
    double boostedCellSafetyFactor = 20.0;
    double boostedCellNeighborScale = 1.0;
    int targetedRefinementLevels = 0;
    double targetedRefinementNeighborScale = 0.0;
    lowmc::EnvelopeFailurePolicy failurePolicy = lowmc::EnvelopeFailurePolicy::Diagnostic;
    double kTMin = 0.002;
    double kTMax = 0.2;
    double etaMin = -5.0;
    double etaMax = 5.0;
    double phiMin = 0.0;
    double phiMax = 2.0 * std::acos(-1.0);
    lowmc::GridAdaptiveEtaPhiProposalConfig gridConfig = lowmc::tunedMassGridAdaptiveEtaPhiConfig();
    std::string samplesCsv;
};

void printUsage(std::ostream& out, const char* programName) {
    out << "Usage: " << programName << " [options]\n"
        << "\n"
        << "Run the generator-independent shape-unweighting API on a two-pion toy event.\n"
        << "This diagnostic does not write HepMC3 output events and does not model photon multiplicity.\n"
        << "\n"
        << "Options:\n"
        << "  --trials N               Trial photons for accept/reject (default: 1000)\n"
        << "  --calibration-trials N   Random calibration trials for hybrid bound (default: 1000)\n"
        << "  --seed N                 Random seed (default: 12345)\n"
        << "  --topology NAME          Topology: two-pion, event4, forward4, or multi8 (default: two-pion)\n"
        << "  --particles PATH         Read particles with columns: charge inOutSign mass p eta phi\n"
        << "  --momentum GeV           Charged-pion momentum (default: 100)\n"
        << "  --particle-eta X         First charged-pion pseudorapidity (default: 0)\n"
        << "  --particle-phi X         First charged-pion azimuth (default: pi/2)\n"
        << "  --second-momentum GeV    Optional second charged-pion momentum\n"
        << "  --second-particle-eta X  Optional second charged-pion pseudorapidity\n"
        << "  --second-particle-phi X  Optional second charged-pion azimuth\n"
        << "  --safety-factor X        Envelope safety factor (default: 3)\n"
        << "  --probe-bins N           Probe lattice bins per cell dimension (default: 11)\n"
        << "  --probe-rule RULE        Probe rule: interior or edge (default: interior)\n"
        << "  --probe-policy POLICY    Probe policy: fixed or boosted-cells (default: fixed)\n"
        << "  --boosted-probe-bins N   Dense probe bins for promoted cells (default: 21)\n"
        << "  --boosted-safety-factor X  Safety factor for promoted cells (default: 20)\n"
        << "  --boosted-neighbor-scale X  Charged-direction cell-neighborhood scale (default: 1)\n"
        << "  --targeted-refinement-levels N  Extra charged-direction split levels (default: 0)\n"
        << "  --targeted-refinement-neighbor-scale X  Extra split neighborhood scale (default: 0)\n"
        << "  --failure-policy POLICY  Failure policy: diagnostic or strict (default: diagnostic)\n"
        << "  --kt-min GeV             Minimum photon k_T (default: 0.002)\n"
        << "  --kt-max GeV             Maximum photon k_T (default: 0.2)\n"
        << "  --eta-min X              Minimum photon eta (default: -5)\n"
        << "  --eta-max X              Maximum photon eta (default: 5)\n"
        << "  --phi-min X              Minimum photon phi (default: 0)\n"
        << "  --phi-max X              Maximum photon phi (default: 2 pi)\n"
        << "  --grid-eta-bins N        Initial eta bins for grid proposal (default: 32)\n"
        << "  --grid-phi-bins N        Initial phi bins for grid proposal (default: 32)\n"
        << "  --grid-refinement-levels N  Grid refinement levels (default: 3)\n"
        << "  --grid-refine-fraction X    Fraction of cells refined per level (default: 0.1)\n"
        << "  --grid-uniform-fraction X   Uniform proposal floor fraction (default: 0.1)\n"
        << "  --samples-csv PATH       Optional CSV with accepted candidates\n";
}

int parseNonNegativeInt(const std::string& text, const std::string& optionName) {
    std::size_t parsed = 0;
    const int value = std::stoi(text, &parsed);
    if (parsed != text.size() || value < 0) {
        throw std::invalid_argument(optionName + " must be a non-negative integer");
    }
    return value;
}

int parsePositiveInt(const std::string& text, const std::string& optionName) {
    const int value = parseNonNegativeInt(text, optionName);
    if (value <= 0) {
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

std::string readStringOption(int& index, int argc, char** argv, const std::string& optionName) {
    if (index + 1 >= argc) {
        throw std::invalid_argument(optionName + " requires a value");
    }
    return argv[++index];
}

lowmc::ShapeUnweightingProbeRule parseProbeRule(const std::string& text) {
    if (text == "interior") {
        return lowmc::ShapeUnweightingProbeRule::Interior;
    }
    if (text == "edge") {
        return lowmc::ShapeUnweightingProbeRule::Edge;
    }
    throw std::invalid_argument("--probe-rule must be interior or edge");
}

const char* probeRuleName(lowmc::ShapeUnweightingProbeRule rule) {
    switch (rule) {
        case lowmc::ShapeUnweightingProbeRule::Interior:
            return "interior";
        case lowmc::ShapeUnweightingProbeRule::Edge:
            return "edge";
    }
    return "unknown";
}

lowmc::ShapeUnweightingProbePolicy parseProbePolicy(const std::string& text) {
    if (text == "fixed") {
        return lowmc::ShapeUnweightingProbePolicy::Fixed;
    }
    if (text == "boosted-cells") {
        return lowmc::ShapeUnweightingProbePolicy::BoostedCells;
    }
    throw std::invalid_argument("--probe-policy must be fixed or boosted-cells");
}

const char* probePolicyName(lowmc::ShapeUnweightingProbePolicy policy) {
    switch (policy) {
        case lowmc::ShapeUnweightingProbePolicy::Fixed:
            return "fixed";
        case lowmc::ShapeUnweightingProbePolicy::BoostedCells:
            return "boosted-cells";
    }
    return "unknown";
}

lowmc::EnvelopeFailurePolicy parseFailurePolicy(const std::string& text) {
    if (text == "diagnostic") {
        return lowmc::EnvelopeFailurePolicy::Diagnostic;
    }
    if (text == "strict") {
        return lowmc::EnvelopeFailurePolicy::Strict;
    }
    throw std::invalid_argument("--failure-policy must be diagnostic or strict");
}

const char* failurePolicyName(lowmc::EnvelopeFailurePolicy policy) {
    switch (policy) {
        case lowmc::EnvelopeFailurePolicy::Diagnostic:
            return "diagnostic";
        case lowmc::EnvelopeFailurePolicy::Strict:
            return "strict";
    }
    return "unknown";
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

Options parseOptions(int argc, char** argv) {
    Options options;
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument == "--help" || argument == "-h") {
            printUsage(std::cout, argv[0]);
            std::exit(EXIT_SUCCESS);
        }
        if (argument == "--trials") {
            options.trials = parseNonNegativeInt(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--calibration-trials") {
            options.calibrationTrials =
                parseNonNegativeInt(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--seed") {
            options.seed = parseUInt64(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--topology") {
            options.topology = readStringOption(index, argc, argv, argument);
            if (options.topology != "two-pion" && options.topology != "event4" &&
                options.topology != "forward4" && options.topology != "multi8") {
                throw std::invalid_argument(
                    "--topology must be two-pion, event4, forward4, or multi8");
            }
            continue;
        }
        if (argument == "--particles") {
            options.particlesPath = readStringOption(index, argc, argv, argument);
            continue;
        }
        if (argument == "--momentum") {
            options.momentum = parsePositiveDouble(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--particle-eta") {
            options.particleEta = parseDouble(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--particle-phi") {
            options.particlePhi = parseDouble(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--second-momentum") {
            options.secondMomentum =
                parsePositiveDouble(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--second-particle-eta") {
            options.secondParticleEta =
                parseDouble(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--second-particle-phi") {
            options.secondParticlePhi =
                parseDouble(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--safety-factor") {
            options.safetyFactor =
                parsePositiveDouble(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--probe-bins") {
            options.probeBins = parsePositiveInt(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--probe-rule") {
            options.probeRule = parseProbeRule(readStringOption(index, argc, argv, argument));
            continue;
        }
        if (argument == "--probe-policy") {
            options.probePolicy = parseProbePolicy(readStringOption(index, argc, argv, argument));
            continue;
        }
        if (argument == "--boosted-probe-bins") {
            options.boostedProbeBins =
                parsePositiveInt(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--boosted-safety-factor") {
            options.boostedCellSafetyFactor =
                parsePositiveDouble(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--boosted-neighbor-scale") {
            options.boostedCellNeighborScale =
                parseDouble(readStringOption(index, argc, argv, argument), argument);
            if (options.boostedCellNeighborScale < 0.0) {
                throw std::invalid_argument(argument + " must be non-negative");
            }
            continue;
        }
        if (argument == "--targeted-refinement-levels") {
            options.targetedRefinementLevels =
                parseNonNegativeInt(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--targeted-refinement-neighbor-scale") {
            options.targetedRefinementNeighborScale =
                parseDouble(readStringOption(index, argc, argv, argument), argument);
            if (options.targetedRefinementNeighborScale < 0.0) {
                throw std::invalid_argument(argument + " must be non-negative");
            }
            continue;
        }
        if (argument == "--failure-policy") {
            options.failurePolicy = parseFailurePolicy(readStringOption(index, argc, argv, argument));
            continue;
        }
        if (argument == "--kt-min") {
            options.kTMin = parsePositiveDouble(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--kt-max") {
            options.kTMax = parsePositiveDouble(readStringOption(index, argc, argv, argument), argument);
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
        if (argument == "--grid-eta-bins") {
            options.gridConfig.etaBins =
                parsePositiveInt(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--grid-phi-bins") {
            options.gridConfig.phiBins =
                parsePositiveInt(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--grid-refinement-levels") {
            options.gridConfig.refinementLevels =
                parseNonNegativeInt(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--grid-refine-fraction") {
            options.gridConfig.refineFraction =
                parsePositiveDouble(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--grid-uniform-fraction") {
            options.gridConfig.uniformFraction =
                parseDouble(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--samples-csv") {
            options.samplesCsv = readStringOption(index, argc, argv, argument);
            continue;
        }
        throw std::invalid_argument("unknown option: " + argument);
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

lowmc::ChargedParticle makePion(
    double momentum,
    double particleEta,
    double particlePhi,
    double charge) {
    const double mass = 0.13957039;
    const double pT = momentum / std::cosh(particleEta);
    const double px = pT * std::cos(particlePhi);
    const double py = pT * std::sin(particlePhi);
    const double pz = pT * std::sinh(particleEta);
    const double energy = std::sqrt(momentum * momentum + mass * mass);
    return {lowmc::FourVector{energy, px, py, pz}, charge, lowmc::InOutSign::Outgoing};
}

std::vector<lowmc::ChargedParticle> makeEvent4Pions() {
    return {
        makePion(80.0, 0.20, 0.70, +1.0),
        makePion(55.0, -0.80, 2.40, -1.0),
        makePion(25.0, 1.60, 4.30, +1.0),
        makePion(12.0, -2.20, 5.40, -1.0),
    };
}

std::vector<lowmc::ChargedParticle> makeForward4Pions() {
    return {
        makePion(350.0, 4.50, 0.40, +1.0),
        makePion(120.0, 2.80, 1.35, -1.0),
        makePion(45.0, -0.60, 3.30, +1.0),
        makePion(30.0, -2.00, 5.10, -1.0),
    };
}

std::vector<lowmc::ChargedParticle> makeMulti8Pions() {
    return {
        makePion(140.0, 1.20, 0.25, +1.0),
        makePion(95.0, -1.10, 0.95, -1.0),
        makePion(70.0, 2.40, 1.80, +1.0),
        makePion(50.0, -2.70, 2.55, -1.0),
        makePion(35.0, 0.30, 3.35, +1.0),
        makePion(25.0, -0.40, 4.10, -1.0),
        makePion(18.0, 3.10, 5.00, +1.0),
        makePion(12.0, -3.20, 5.75, -1.0),
    };
}

std::vector<lowmc::ChargedParticle> makeToyPions(const Options& options) {
    if (!options.particlesPath.empty()) {
        return lowmc_diagnostics::readChargedParticlesFromTextFile(options.particlesPath);
    }
    if (options.topology == "event4") {
        return makeEvent4Pions();
    }
    if (options.topology == "forward4") {
        return makeForward4Pions();
    }
    if (options.topology == "multi8") {
        return makeMulti8Pions();
    }

    if (!options.secondMomentum && !options.secondParticleEta && !options.secondParticlePhi) {
        return makeBackToBackPions(
            options.momentum, options.particleEta, options.particlePhi);
    }

    const double secondMomentum = options.secondMomentum.value_or(options.momentum);
    const double secondEta = options.secondParticleEta.value_or(-options.particleEta);
    const double secondPhi =
        options.secondParticlePhi.value_or(options.particlePhi + std::acos(-1.0));

    return {
        makePion(options.momentum, options.particleEta, options.particlePhi, +1.0),
        makePion(secondMomentum, secondEta, secondPhi, -1.0),
    };
}

void writeSamplesCsv(
    const std::string& path,
    const std::vector<lowmc::ShapeUnweightedPhotonCandidate>& candidates) {
    std::ofstream output(path);
    if (!output) {
        throw std::runtime_error("could not open samples CSV for writing: " + path);
    }
    output << "kT,eta,phi,weight,E,px,py,pz\n";
    for (const lowmc::ShapeUnweightedPhotonCandidate& candidate : candidates) {
        output << candidate.kT << ','
               << candidate.eta << ','
               << candidate.phi << ','
               << candidate.photon.weight << ','
               << candidate.photon.momentum.E << ','
               << candidate.photon.momentum.px << ','
               << candidate.photon.momentum.py << ','
               << candidate.photon.momentum.pz << '\n';
    }
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Options options = parseOptions(argc, argv);
        lowmc::ShapeUnweightingConfig config;
        config.phaseSpace = lowmc::PhotonPhaseSpace{
            options.kTMin,
            options.kTMax,
            options.etaMin,
            options.etaMax,
            options.phiMin,
            options.phiMax,
        };
        config.trials = options.trials;
        config.seed = options.seed;
        config.gridConfig = options.gridConfig;
        config.calibrationTrials = options.calibrationTrials;
        config.probeBins = options.probeBins;
        config.probeRule = options.probeRule;
        config.probePolicy = options.probePolicy;
        config.boostedProbeBins = options.boostedProbeBins;
        config.boostedCellSafetyFactor = options.boostedCellSafetyFactor;
        config.boostedCellNeighborScale = options.boostedCellNeighborScale;
        config.targetedRefinementLevels = options.targetedRefinementLevels;
        config.targetedRefinementNeighborScale =
            options.targetedRefinementNeighborScale;
        config.safetyFactor = options.safetyFactor;
        config.failurePolicy = options.failurePolicy;

        const std::vector<lowmc::ChargedParticle> selectedParticles = makeToyPions(options);
        const lowmc::ShapeUnweightingResult result =
            lowmc::generateShapeUnweightedPhotonCandidates(selectedParticles, config);

        if (!options.samplesCsv.empty()) {
            writeSamplesCsv(options.samplesCsv, result.candidates);
        }

        std::cout << "status=" << statusName(result.status)
                  << " failure_policy=" << failurePolicyName(options.failurePolicy)
                  << " trials=" << result.envelope.trialCount
                  << " calibration_trials=" << result.envelope.calibrationCount
                  << " accepted_candidates=" << result.candidates.size()
                  << " accepted_count=" << result.envelope.acceptedCount
                  << " failed_envelope_count=" << result.envelope.failedEnvelopeCount
                  << " max_trial_bound_ratio=" << result.envelope.maxTrialBoundRatio
                  << " observed_acceptance_efficiency="
                  << result.envelope.observedAcceptanceEfficiency
                  << " expected_acceptance_efficiency="
                  << result.envelope.expectedAcceptanceEfficiency
                  << " envelope_integral=" << result.envelope.envelopeIntegral
                  << " probe_evaluations=" << result.envelope.probeEvaluations
                  << " grid_cells=" << result.envelope.gridCellCount
                  << " suspect_cell_count=" << result.envelope.suspectCellCount
                  << " promoted_cell_count=" << result.envelope.promotedCellCount
                  << " fallback_count=" << result.envelope.fallbackCount
                  << " safety_factor=" << result.envelope.safetyFactor
                  << " boosted_safety_factor=" << result.envelope.boostedCellSafetyFactor
                  << " boosted_neighbor_scale=" << result.envelope.boostedCellNeighborScale
                  << " targeted_refinement_levels=" << result.envelope.targetedRefinementLevels
                  << " targeted_refinement_neighbor_scale="
                  << result.envelope.targetedRefinementNeighborScale
                  << " probe_bins=" << result.envelope.probeBins
                  << " boosted_probe_bins=" << result.envelope.boostedProbeBins
                  << " probe_rule=" << probeRuleName(result.envelope.probeRule)
                  << " probe_policy=" << probePolicyName(result.envelope.probePolicy)
                  << " topology=" << (options.particlesPath.empty() ? options.topology : "particles")
                  << " particle_file="
                  << (options.particlesPath.empty() ? "none" : options.particlesPath)
                  << " particle_count=" << selectedParticles.size()
                  << " momentum=" << options.momentum
                  << " particle_eta=" << options.particleEta
                  << " particle_phi=" << options.particlePhi
                  << " second_momentum=" << options.secondMomentum.value_or(options.momentum)
                  << " second_particle_eta="
                  << options.secondParticleEta.value_or(-options.particleEta)
                  << " second_particle_phi="
                  << options.secondParticlePhi.value_or(options.particlePhi + std::acos(-1.0))
                  << " grid_build_kt=" << result.gridBuildKT
                  << " signed_charge_sum=" << result.chargeConservation.netCharge
                  << '\n';
    } catch (const std::exception& error) {
        std::cerr << "FAILED: " << error.what() << "\n\n";
        printUsage(std::cerr, argv[0]);
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
