#pragma once

#include "lowmc/LowMC.h"

#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace lowmc_diagnostics {

inline lowmc::InOutSign parseInOutSign(const std::string& text) {
    if (text == "outgoing" || text == "out" || text == "+") {
        return lowmc::InOutSign::Outgoing;
    }
    if (text == "incoming" || text == "in" || text == "-") {
        return lowmc::InOutSign::Incoming;
    }
    throw std::invalid_argument("inOutSign must be incoming or outgoing");
}

inline lowmc::ChargedParticle makeChargedParticleFromMassPEtaPhi(
    double charge,
    lowmc::InOutSign sign,
    double mass,
    double momentum,
    double eta,
    double phi) {
    if (!std::isfinite(charge) || charge == 0.0) {
        throw std::invalid_argument("charge must be finite and non-zero");
    }
    if (!std::isfinite(mass) || mass < 0.0) {
        throw std::invalid_argument("mass must be finite and non-negative");
    }
    if (!std::isfinite(momentum) || momentum < 0.0) {
        throw std::invalid_argument("momentum must be finite and non-negative");
    }
    if (!std::isfinite(eta) || !std::isfinite(phi)) {
        throw std::invalid_argument("eta and phi must be finite");
    }

    const double pT = momentum / std::cosh(eta);
    const double px = pT * std::cos(phi);
    const double py = pT * std::sin(phi);
    const double pz = pT * std::sinh(eta);
    const double energy = std::sqrt(momentum * momentum + mass * mass);
    return {lowmc::FourVector{energy, px, py, pz}, charge, sign};
}

inline std::vector<lowmc::ChargedParticle> readChargedParticlesFromTextFile(
    const std::string& path) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("could not open particle input file: " + path);
    }

    std::vector<lowmc::ChargedParticle> particles;
    std::string line;
    int lineNumber = 0;
    while (std::getline(input, line)) {
        ++lineNumber;
        const std::size_t comment = line.find('#');
        if (comment != std::string::npos) {
            line.erase(comment);
        }

        std::istringstream stream(line);
        double charge = 0.0;
        std::string signText;
        double mass = 0.0;
        double momentum = 0.0;
        double eta = 0.0;
        double phi = 0.0;
        if (!(stream >> charge >> signText >> mass >> momentum >> eta >> phi)) {
            if (line.find_first_not_of(" \t\r\n") == std::string::npos) {
                continue;
            }
            throw std::runtime_error(
                path + ":" + std::to_string(lineNumber) +
                ": expected columns: charge inOutSign mass p eta phi");
        }

        std::string extra;
        if (stream >> extra) {
            throw std::runtime_error(
                path + ":" + std::to_string(lineNumber) +
                ": unexpected extra column after phi");
        }

        try {
            particles.push_back(makeChargedParticleFromMassPEtaPhi(
                charge, parseInOutSign(signText), mass, momentum, eta, phi));
        } catch (const std::exception& error) {
            throw std::runtime_error(
                path + ":" + std::to_string(lineNumber) + ": " + error.what());
        }
    }

    if (particles.empty()) {
        throw std::runtime_error("particle input file contains no particles: " + path);
    }
    return particles;
}

}  // namespace lowmc_diagnostics
