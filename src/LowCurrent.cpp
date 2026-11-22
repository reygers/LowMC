#include "lowmc/LowCurrent.h"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace lowmc {

LowCurrent::LowCurrent(std::vector<ChargedParticle> particles)
    : particles_(std::move(particles)) {}

const std::vector<ChargedParticle>& LowCurrent::particles() const {
    return particles_;
}

FourVector LowCurrent::current(const FourVector& photonMomentum) const {
    return current(photonMomentum, VirtualCurrentModel::LeadingLow);
}

FourVector LowCurrent::current(
    const FourVector& photonMomentum,
    VirtualCurrentModel model) const {
    FourVector result;

    for (const ChargedParticle& particle : particles_) {
        const double sign = static_cast<double>(signValue(particle.inOutSign));
        if (model == VirtualCurrentModel::LichardEq2_8) {
            const double denominator =
                2.0 * particle.momentum.dot(photonMomentum) +
                sign * photonMomentum.mass2();
            if (!std::isfinite(denominator) || denominator == 0.0) {
                throw std::invalid_argument(
                    "Lichard virtual-current denominator must be finite and nonzero");
            }
            const FourVector numerator =
                2.0 * particle.momentum + sign * photonMomentum;
            result += sign * particle.charge / denominator * numerator;
        } else {
            const double denominator = particle.momentum.dot(photonMomentum);
            const double coefficient = sign * particle.charge / denominator;
            result += coefficient * particle.momentum;
        }
    }

    return result;
}

double LowCurrent::intensity(const FourVector& photonMomentum) const {
    return intensity(photonMomentum, VirtualCurrentModel::LeadingLow);
}

double LowCurrent::intensity(
    const FourVector& photonMomentum,
    VirtualCurrentModel model) const {
    const FourVector J = current(photonMomentum, model);
    return -J.dot(J);
}

double LowCurrent::intensityPairSum(const FourVector& photonMomentum) const {
    double sum = 0.0;

    for (const ChargedParticle& first : particles_) {
        const double firstDenominator = first.momentum.dot(photonMomentum);
        const double firstCharge =
            static_cast<double>(signValue(first.inOutSign)) * first.charge;

        for (const ChargedParticle& second : particles_) {
            const double secondDenominator = second.momentum.dot(photonMomentum);
            const double secondCharge =
                static_cast<double>(signValue(second.inOutSign)) * second.charge;

            sum += firstCharge * secondCharge *
                   (-first.momentum.dot(second.momentum)) /
                   (firstDenominator * secondDenominator);
        }
    }

    return sum;
}

ChargeConservationDiagnostic LowCurrent::chargeConservation(double tolerance) const {
    ChargeConservationDiagnostic diagnostic;

    for (const ChargedParticle& particle : particles_) {
        diagnostic.netCharge +=
            static_cast<double>(signValue(particle.inOutSign)) * particle.charge;
    }

    diagnostic.conserved = std::abs(diagnostic.netCharge) <= tolerance;
    return diagnostic;
}

double lowPhotonDensity(
    const LowCurrent& current,
    const FourVector& photonMomentum,
    double alpha) {
    const double twoPi = 2.0 * std::acos(-1.0);
    return alpha / (twoPi * twoPi) * current.intensity(photonMomentum);
}

}  // namespace lowmc
