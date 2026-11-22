#include "lowmc/Balek.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace lowmc {
namespace {

double pi() {
    return std::acos(-1.0);
}

void validateParameters(const BalekParameters& parameters) {
    if (!std::isfinite(parameters.sigmaInelastic) ||
        parameters.sigmaInelastic <= 0.0 ||
        !std::isfinite(parameters.meanChargeTransferSquared) ||
        parameters.meanChargeTransferSquared < 0.0 ||
        !std::isfinite(parameters.electronMass) ||
        parameters.electronMass <= 0.0 ||
        !std::isfinite(parameters.alpha) || parameters.alpha <= 0.0) {
        throw std::invalid_argument("invalid Balek parameters");
    }
}

void validateMass(double mass, double electronMass) {
    if (!std::isfinite(mass) || mass < 2.0 * electronMass) {
        throw std::invalid_argument(
            "Balek pair mass must be finite and at or above threshold");
    }
}

double normalization(const BalekParameters& parameters) {
    return parameters.sigmaInelastic * 4.0 * parameters.alpha * parameters.alpha *
        parameters.meanChargeTransferSquared / (3.0 * pi() * pi());
}

double leptonMassFactor(double mass, double electronMass) {
    const double ratio = electronMass * electronMass / (mass * mass);
    return (1.0 + 2.0 * ratio) *
        std::sqrt(std::max(0.0, 1.0 - 4.0 * ratio));
}

}  // namespace

double balekDifferentialCrossSection(
    double mass,
    double qT,
    const BalekParameters& parameters) {
    validateParameters(parameters);
    validateMass(mass, parameters.electronMass);
    if (!std::isfinite(qT) || qT < 0.0) {
        throw std::invalid_argument("Balek qT must be finite and non-negative");
    }
    return normalization(parameters) * leptonMassFactor(mass, parameters.electronMass) *
        qT / (mass * (qT * qT + mass * mass));
}

double balekMassProjection(
    double mass,
    double qTMin,
    double qTMax,
    const BalekParameters& parameters) {
    validateParameters(parameters);
    validateMass(mass, parameters.electronMass);
    if (!std::isfinite(qTMin) || !std::isfinite(qTMax) || qTMin < 0.0 ||
        qTMax <= qTMin) {
        throw std::invalid_argument("invalid Balek qT integration range");
    }
    const double mass2 = mass * mass;
    const double logarithm =
        std::log((qTMax * qTMax + mass2) / (qTMin * qTMin + mass2));
    return normalization(parameters) * leptonMassFactor(mass, parameters.electronMass) *
        logarithm / (2.0 * mass);
}

double balekQTProjectionMassless(
    double qT,
    double massMin,
    double massMax,
    const BalekParameters& parameters) {
    validateParameters(parameters);
    validateMass(massMin, parameters.electronMass);
    validateMass(massMax, parameters.electronMass);
    if (!std::isfinite(qT) || qT < 0.0 || massMax <= massMin) {
        throw std::invalid_argument("invalid Balek mass integration range or qT");
    }
    if (qT == 0.0) {
        return 0.0;
    }
    const double qT2 = qT * qT;
    const double min2 = massMin * massMin;
    const double max2 = massMax * massMax;
    const double logarithm =
        std::log(max2 * (qT2 + min2) / (min2 * (qT2 + max2)));
    return normalization(parameters) * logarithm / (2.0 * qT);
}

}  // namespace lowmc
