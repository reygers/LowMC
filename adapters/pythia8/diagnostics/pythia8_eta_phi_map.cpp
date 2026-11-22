#include "lowmc/LowMC.h"
#include "lowmc/pythia8/Pythia8Adapter.h"

#include <Pythia8/Pythia.h>

#include <cmath>
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

constexpr double Pi = 3.14159265358979323846;

struct Options {
    std::string gridPath;
    std::string particlesPath;
    int seed = 12345;
    double eCM = 13000.0;
    double kT = 0.02;
    double etaMin = -3.0;
    double etaMax = 3.0;
    double phiMin = -Pi;
    double phiMax = Pi;
    int etaBins = 241;
    int phiBins = 321;
};

void printUsage(std::ostream& out, const char* program) {
    out << "Usage: " << program << " GRID.csv PARTICLES.csv [options]\n"
        << "\n"
        << "Generate one non-diffractive Pythia8 pp event and evaluate the\n"
        << "event-specific leading-Low angular yield on an eta-phi grid.\n"
        << "\n"
        << "Options:\n"
        << "  --seed N          Pythia random seed (default: 12345)\n"
        << "  --ecm GEV         Center-of-mass energy (default: 13000)\n"
        << "  --kt GEV          Fixed kT used for evaluation (default: 0.02)\n"
        << "  --eta-min X       Minimum eta (default: -3)\n"
        << "  --eta-max X       Maximum eta (default: 3)\n"
        << "  --eta-bins N      Number of eta grid points (default: 241)\n"
        << "  --phi-bins N      Number of phi grid points (default: 321)\n";
}

std::string nextValue(int& index, int argc, char** argv, const std::string& option) {
    if (++index >= argc) {
        throw std::invalid_argument(option + " requires a value");
    }
    return argv[index];
}

Options parseOptions(int argc, char** argv) {
    if (argc < 3) {
        printUsage(std::cerr, argv[0]);
        throw std::invalid_argument("GRID.csv and PARTICLES.csv are required");
    }

    Options options;
    options.gridPath = argv[1];
    options.particlesPath = argv[2];
    for (int index = 3; index < argc; ++index) {
        const std::string option = argv[index];
        if (option == "--help" || option == "-h") {
            printUsage(std::cout, argv[0]);
            std::exit(EXIT_SUCCESS);
        }
        if (option == "--seed") {
            options.seed = std::stoi(nextValue(index, argc, argv, option));
        } else if (option == "--ecm") {
            options.eCM = std::stod(nextValue(index, argc, argv, option));
        } else if (option == "--kt") {
            options.kT = std::stod(nextValue(index, argc, argv, option));
        } else if (option == "--eta-min") {
            options.etaMin = std::stod(nextValue(index, argc, argv, option));
        } else if (option == "--eta-max") {
            options.etaMax = std::stod(nextValue(index, argc, argv, option));
        } else if (option == "--eta-bins") {
            options.etaBins = std::stoi(nextValue(index, argc, argv, option));
        } else if (option == "--phi-bins") {
            options.phiBins = std::stoi(nextValue(index, argc, argv, option));
        } else {
            throw std::invalid_argument("unknown option: " + option);
        }
    }

    if (options.seed <= 0 || options.eCM <= 0.0 || options.kT <= 0.0 ||
        options.etaMin >= options.etaMax || options.etaBins < 2 || options.phiBins < 2) {
        throw std::invalid_argument("invalid numeric option");
    }
    return options;
}

lowmc::FourVector photonMomentum(double kT, double eta, double phi) {
    return lowmc::FourVector{
        kT * std::cosh(eta),
        kT * std::cos(phi),
        kT * std::sin(phi),
        kT * std::sinh(eta),
    };
}

double pseudorapidity(const lowmc::FourVector& momentum) {
    const double pt = momentum.pt();
    if (pt == 0.0) {
        return std::copysign(std::numeric_limits<double>::infinity(), momentum.pz);
    }
    return std::asinh(momentum.pz / pt);
}

double wrappedPhi(const lowmc::FourVector& momentum) {
    return std::atan2(momentum.py, momentum.px);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Options options = parseOptions(argc, argv);

        Pythia8::Pythia pythia(LOWMC_PYTHIA8_XMLDOC, false);
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

        int attempts = 0;
        while (!pythia.next()) {
            if (++attempts > 100) {
                throw std::runtime_error("Pythia did not produce an accepted event");
            }
        }

        const lowmc::AdapterParticleSelectionResult selection =
            lowmc::pythia8::selectChargedParticles(pythia.event);
        const lowmc::LowCurrent current(selection.selectedParticles);
        const lowmc::ChargeConservationDiagnostic charge = current.chargeConservation();
        if (!charge.conserved) {
            throw std::runtime_error("selected event does not conserve charge");
        }

        std::ofstream grid(options.gridPath);
        if (!grid) {
            throw std::runtime_error("could not open grid output");
        }
        grid << "eta,phi,dN_dlnkT_deta_dphi\n" << std::setprecision(12);
        for (int etaIndex = 0; etaIndex < options.etaBins; ++etaIndex) {
            const double eta = options.etaMin +
                (options.etaMax - options.etaMin) * etaIndex / (options.etaBins - 1.0);
            for (int phiIndex = 0; phiIndex < options.phiBins; ++phiIndex) {
                const double phi = options.phiMin +
                    (options.phiMax - options.phiMin) * phiIndex / (options.phiBins - 1.0);
                const lowmc::FourVector photon = photonMomentum(options.kT, eta, phi);
                const double angularYield =
                    options.kT * options.kT * lowmc::lowPhotonDensity(current, photon);
                grid << eta << ',' << phi << ',' << angularYield << '\n';
            }
        }

        std::ofstream particles(options.particlesPath);
        if (!particles) {
            throw std::runtime_error("could not open particle output");
        }
        particles << "eta,phi,pt,charge,in_out\n" << std::setprecision(12);
        for (const lowmc::ChargedParticle& particle : selection.selectedParticles) {
            if (particle.inOutSign == lowmc::InOutSign::Incoming || particle.momentum.pt() == 0.0) {
                continue;
            }
            particles << pseudorapidity(particle.momentum) << ','
                      << wrappedPhi(particle.momentum) << ','
                      << particle.momentum.pt() << ','
                      << particle.charge << ",outgoing\n";
        }

        std::cout << "seed=" << options.seed
                  << " selected_charged=" << selection.selectedParticles.size()
                  << " outgoing_charged=" << selection.diagnostic.selectedOutgoingCharged
                  << " signed_charge_sum=" << selection.diagnostic.signedChargeSum
                  << " grid=" << options.etaBins << 'x' << options.phiBins << '\n';
    } catch (const std::exception& error) {
        std::cerr << "FAILED: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
