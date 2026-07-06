#include "ParticleInput.h"

#include "lowmc/LowMC.h"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr double DiagnosticKT = 0.02;

lowmc::FourVector photonFromKtEtaPhi(double kT, double eta, double phi) {
    const double px = kT * std::cos(phi);
    const double py = kT * std::sin(phi);
    const double pz = kT * std::sinh(eta);
    const double energy = kT * std::cosh(eta);
    return lowmc::FourVector{energy, px, py, pz};
}

struct Options {
    std::string outputPath = "two_particle_eta_phi.csv";
    std::string topology = "two-pion";
    std::string particlesPath;
    double momentum = 1.0;
    double particleEta = 0.0;
    double particlePhi = 0.5 * std::acos(-1.0);
    std::optional<double> secondMomentum;
    std::optional<double> secondParticleEta;
    std::optional<double> secondParticlePhi;
    double etaMin = -5.0;
    double etaMax = 5.0;
    double phiMin = 0.0;
    double phiMax = 2.0 * std::acos(-1.0);
    int etaBins = 101;
    int phiBins = 121;
};

const char* validationErrorMessage(lowmc::PhotonPhaseSpaceValidationError error) {
    switch (error) {
        case lowmc::PhotonPhaseSpaceValidationError::None:
            return "phase space is valid";
        case lowmc::PhotonPhaseSpaceValidationError::NonFiniteBound:
            return "phase-space bounds must be finite";
        case lowmc::PhotonPhaseSpaceValidationError::NonPositiveKTMin:
            return "kTMin must be positive";
        case lowmc::PhotonPhaseSpaceValidationError::EmptyKTRange:
            return "kTMin must be smaller than kTMax";
        case lowmc::PhotonPhaseSpaceValidationError::EmptyEtaRange:
            return "etaMin must be smaller than etaMax";
        case lowmc::PhotonPhaseSpaceValidationError::EmptyPhiRange:
            return "phiMin must be smaller than phiMax";
    }
    return "unknown phase-space validation error";
}

lowmc::PhotonPhaseSpace makeDiagnosticPhaseSpace(const Options& options) {
    return lowmc::PhotonPhaseSpace{
        DiagnosticKT,
        std::nextafter(DiagnosticKT, std::numeric_limits<double>::infinity()),
        options.etaMin,
        options.etaMax,
        options.phiMin,
        options.phiMax,
    };
}

void printUsage(const char* programName) {
    std::cout << "Usage: " << programName
              << " [output.csv] [--momentum GeV] [--particle-eta eta]\n"
              << "       [--topology two-pion|event4|forward4|multi8]\n"
              << "       [--particles file]\n"
              << "       [--particle-phi phi] [--second-momentum GeV]\n"
              << "       [--second-particle-eta eta] [--second-particle-phi phi]\n"
              << "       [--eta-min eta] [--eta-max eta]\n"
              << "       [--phi-min phi] [--phi-max phi]\n"
              << "       [--eta-bins n] [--phi-bins n]\n\n"
              << "Create an eta,phi,S CSV for a charged-pion toy event.\n"
              << "Supplying any --second-* option switches to an explicit two-pion topology.\n"
              << "--topology event4 uses a fixed charge-conserving four-pion topology.\n"
              << "--topology forward4 uses a four-pion topology with one very forward particle.\n"
              << "--topology multi8 uses a fixed charge-conserving eight-pion topology.\n"
              << "--particles reads columns: charge inOutSign mass p eta phi.\n"
              << "The default momentum is 1 GeV and the default direction is eta=0, phi=pi/2.\n"
              << "The default grid is 101 eta points by 121 phi points.\n";
}

double parseDouble(const std::string& text, const std::string& optionName) {
    std::size_t parsed = 0;
    const double value = std::stod(text, &parsed);
    if (parsed != text.size()) {
        throw std::invalid_argument(optionName + " must be a number");
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

int parseGridBins(const std::string& text, const std::string& optionName) {
    std::size_t parsed = 0;
    const int value = std::stoi(text, &parsed);
    if (parsed != text.size() || value < 2) {
        throw std::invalid_argument(optionName + " must be an integer >= 2");
    }
    return value;
}

double readOptionValue(int& index, int argc, char** argv, const std::string& optionName) {
    if (index + 1 >= argc) {
        throw std::invalid_argument(optionName + " requires a value");
    }
    return parseDouble(argv[++index], optionName);
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
            printUsage(argv[0]);
            std::exit(EXIT_SUCCESS);
        }

        if (argument == "--momentum") {
            if (index + 1 >= argc) {
                throw std::invalid_argument("--momentum requires a value in GeV");
            }
            options.momentum = parsePositiveDouble(argv[++index], "--momentum");
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
        if (argument == "--particle-eta") {
            options.particleEta = readOptionValue(index, argc, argv, "--particle-eta");
            continue;
        }
        if (argument == "--particle-phi") {
            options.particlePhi = readOptionValue(index, argc, argv, "--particle-phi");
            continue;
        }
        if (argument == "--second-momentum") {
            if (index + 1 >= argc) {
                throw std::invalid_argument("--second-momentum requires a value in GeV");
            }
            options.secondMomentum = parsePositiveDouble(argv[++index], "--second-momentum");
            continue;
        }
        if (argument == "--second-particle-eta") {
            options.secondParticleEta = readOptionValue(index, argc, argv, "--second-particle-eta");
            continue;
        }
        if (argument == "--second-particle-phi") {
            options.secondParticlePhi = readOptionValue(index, argc, argv, "--second-particle-phi");
            continue;
        }
        if (argument == "--eta-min") {
            options.etaMin = readOptionValue(index, argc, argv, "--eta-min");
            continue;
        }
        if (argument == "--eta-max") {
            options.etaMax = readOptionValue(index, argc, argv, "--eta-max");
            continue;
        }
        if (argument == "--phi-min") {
            options.phiMin = readOptionValue(index, argc, argv, "--phi-min");
            continue;
        }
        if (argument == "--phi-max") {
            options.phiMax = readOptionValue(index, argc, argv, "--phi-max");
            continue;
        }
        if (argument == "--eta-bins") {
            options.etaBins =
                parseGridBins(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (argument == "--phi-bins") {
            options.phiBins =
                parseGridBins(readStringOption(index, argc, argv, argument), argument);
            continue;
        }

        if (!argument.empty() && argument[0] == '-') {
            throw std::invalid_argument("unknown option: " + argument);
        }

        if (options.outputPath != "two_particle_eta_phi.csv") {
            throw std::invalid_argument("only one output CSV path can be specified");
        }
        options.outputPath = argument;
    }

    return options;
}

lowmc::LowCurrent makeBackToBackToyCurrent(
    double momentum, double particleEta, double particlePhi) {
    const double mass = 0.13957039;
    const double pt = momentum / std::cosh(particleEta);
    const double px = pt * std::cos(particlePhi);
    const double py = pt * std::sin(particlePhi);
    const double pz = pt * std::sinh(particleEta);
    const double energy = std::sqrt(momentum * momentum + mass * mass);

    return lowmc::LowCurrent({
        {lowmc::FourVector{energy, px, py, pz}, +1.0, lowmc::InOutSign::Outgoing},
        {lowmc::FourVector{energy, -px, -py, -pz}, -1.0, lowmc::InOutSign::Outgoing},
    });
}

lowmc::ChargedParticle makePion(
    double momentum,
    double particleEta,
    double particlePhi,
    double charge) {
    const double mass = 0.13957039;
    const double pt = momentum / std::cosh(particleEta);
    const double px = pt * std::cos(particlePhi);
    const double py = pt * std::sin(particlePhi);
    const double pz = pt * std::sinh(particleEta);
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

lowmc::LowCurrent makeToyCurrent(const Options& options) {
    if (!options.particlesPath.empty()) {
        return lowmc::LowCurrent(
            lowmc_diagnostics::readChargedParticlesFromTextFile(options.particlesPath));
    }
    if (options.topology == "event4") {
        return lowmc::LowCurrent(makeEvent4Pions());
    }
    if (options.topology == "forward4") {
        return lowmc::LowCurrent(makeForward4Pions());
    }
    if (options.topology == "multi8") {
        return lowmc::LowCurrent(makeMulti8Pions());
    }

    if (!options.secondMomentum && !options.secondParticleEta && !options.secondParticlePhi) {
        return makeBackToBackToyCurrent(
            options.momentum, options.particleEta, options.particlePhi);
    }

    const double secondMomentum = options.secondMomentum.value_or(options.momentum);
    const double secondEta = options.secondParticleEta.value_or(-options.particleEta);
    const double secondPhi =
        options.secondParticlePhi.value_or(options.particlePhi + std::acos(-1.0));

    return lowmc::LowCurrent({
        makePion(options.momentum, options.particleEta, options.particlePhi, +1.0),
        makePion(secondMomentum, secondEta, secondPhi, -1.0),
    });
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    try {
        options = parseOptions(argc, argv);
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << "\n\n";
        printUsage(argv[0]);
        return EXIT_FAILURE;
    }

    const lowmc::PhotonPhaseSpace phaseSpace = makeDiagnosticPhaseSpace(options);
    const lowmc::PhotonPhaseSpaceValidation validation =
        lowmc::validatePhotonPhaseSpace(phaseSpace);
    if (!validation.valid()) {
        std::cerr << "Invalid photon phase space: "
                  << validationErrorMessage(validation.error) << '\n';
        return EXIT_FAILURE;
    }

    std::ofstream output(options.outputPath);
    if (!output) {
        std::cerr << "Could not open output file: " << options.outputPath << '\n';
        return EXIT_FAILURE;
    }

    const lowmc::LowCurrent current = makeToyCurrent(options);

    output << "eta,phi,S\n";
    output << std::setprecision(17);

    for (int etaBin = 0; etaBin < options.etaBins; ++etaBin) {
        const double eta =
            phaseSpace.etaMin + (phaseSpace.etaMax - phaseSpace.etaMin) *
                                    static_cast<double>(etaBin) /
                                    static_cast<double>(options.etaBins - 1);

        for (int phiBin = 0; phiBin < options.phiBins; ++phiBin) {
            const double phi =
                phaseSpace.phiMin + (phaseSpace.phiMax - phaseSpace.phiMin) *
                                        static_cast<double>(phiBin) /
                                        static_cast<double>(options.phiBins - 1);
            const lowmc::FourVector photon = photonFromKtEtaPhi(DiagnosticKT, eta, phi);
            const double intensity = current.intensity(photon);
            output << eta << ',' << phi << ',' << intensity << '\n';
        }
    }

    if (!output) {
        std::cerr << "Failed while writing output file: " << options.outputPath << '\n';
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
