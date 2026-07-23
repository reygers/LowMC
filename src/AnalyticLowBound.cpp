#include "lowmc/AnalyticLowBound.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

#include "lowmc/LowCurrent.h"

namespace lowmc {
namespace {

[[nodiscard]] double twoPi() {
    return 2.0 * std::acos(-1.0);
}

[[nodiscard]] bool finiteMomentum(const FourVector& momentum) {
    return std::isfinite(momentum.E) &&
           std::isfinite(momentum.px) &&
           std::isfinite(momentum.py) &&
           std::isfinite(momentum.pz);
}

[[nodiscard]] AnalyticLowLegKinematics makeLegKinematics(
    const ChargedParticle& particle) {
    AnalyticLowLegKinematics result;
    const FourVector& momentum = particle.momentum;
    if (!finiteMomentum(momentum) || !(momentum.E > 0.0) ||
        !std::isfinite(particle.charge) || particle.charge == 0.0) {
        return result;
    }

    const long double energy = momentum.E;
    const long double px = momentum.px;
    const long double py = momentum.py;
    const long double pz = momentum.pz;
    const long double massSquared =
        energy * energy - px * px - py * py - pz * pz;
    if (!(massSquared > 0.0L) ||
        massSquared >
            static_cast<long double>(std::numeric_limits<double>::max())) {
        return result;
    }

    result.massSquared = static_cast<double>(massSquared);
    result.mass = std::sqrt(result.massSquared);
    result.transverseMomentum = std::hypot(momentum.px, momentum.py);
    result.transverseMass =
        std::hypot(result.mass, result.transverseMomentum);
    if (!(result.transverseMass > 0.0) ||
        !std::isfinite(result.transverseMass)) {
        return result;
    }
    result.rapidity = std::asinh(momentum.pz / result.transverseMass);
    result.phi = std::atan2(momentum.py, momentum.px);
    result.globalMinimumDenominatorScale =
        result.massSquared /
        (result.transverseMass + result.transverseMomentum);
    result.valid =
        std::isfinite(result.rapidity) &&
        std::isfinite(result.phi) &&
        result.globalMinimumDenominatorScale > 0.0 &&
        std::isfinite(result.globalMinimumDenominatorScale);
    return result;
}

[[nodiscard]] double wrappedAngleDifference(double left, double right) {
    return std::remainder(left - right, twoPi());
}

[[nodiscard]] double intervalDistance(
    double value,
    double lower,
    double upper) {
    if (value >= lower && value <= upper) {
        return 0.0;
    }
    return std::min(std::abs(value - lower), std::abs(value - upper));
}

[[nodiscard]] double periodicIntervalDistance(
    double value,
    double lower,
    double upper) {
    if (upper - lower >= twoPi()) {
        return 0.0;
    }
    const double midpoint = 0.5 * (lower + upper);
    const double nearestShift =
        std::round((midpoint - value) / twoPi());
    double minimum = std::numeric_limits<double>::infinity();
    for (int offset = -1; offset <= 1; ++offset) {
        const double shifted =
            value + (nearestShift + static_cast<double>(offset)) * twoPi();
        minimum = std::min(
            minimum,
            intervalDistance(shifted, lower, upper));
    }
    return std::min(minimum, std::acos(-1.0));
}

[[nodiscard]] double positiveDenominatorScale(
    const AnalyticLowLegKinematics& leg,
    double deltaEta,
    double deltaPhi) {
    const double sinhHalf = std::sinh(0.5 * deltaEta);
    const double sinHalf = std::sin(0.5 * deltaPhi);
    return leg.globalMinimumDenominatorScale +
           2.0 * leg.transverseMass * sinhHalf * sinhHalf +
           2.0 * leg.transverseMomentum * sinHalf * sinHalf;
}

}  // namespace

AnalyticLowBound::AnalyticLowBound(
    std::vector<ChargedParticle> particles)
    : particles_(std::move(particles)) {
    legs_.reserve(particles_.size());
    for (const ChargedParticle& particle : particles_) {
        legs_.push_back(makeLegKinematics(particle));
        invalidParticleCount_ += legs_.back().valid ? 0 : 1;
    }

    const std::size_t pairCount = particles_.size() * particles_.size();
    signedPairCoefficients_.reserve(pairCount);
    absolutePairCoefficients_.reserve(pairCount);
    for (const ChargedParticle& first : particles_) {
        const double firstCharge =
            static_cast<double>(signValue(first.inOutSign)) *
            first.charge;
        for (const ChargedParticle& second : particles_) {
            const double secondCharge =
                static_cast<double>(signValue(second.inOutSign)) *
                second.charge;
            const double coefficient =
                firstCharge * secondCharge *
                (-first.momentum.dot(second.momentum));
            signedPairCoefficients_.push_back(coefficient);
            absolutePairCoefficients_.push_back(std::abs(coefficient));
        }
    }
}

const std::vector<ChargedParticle>& AnalyticLowBound::particles() const {
    return particles_;
}

const std::vector<AnalyticLowLegKinematics>&
AnalyticLowBound::legKinematics() const {
    return legs_;
}

int AnalyticLowBound::invalidParticleCount() const {
    return invalidParticleCount_;
}

bool AnalyticLowBound::valid() const {
    return !particles_.empty() && invalidParticleCount_ == 0;
}

double AnalyticLowBound::denominatorScale(
    std::size_t particleIndex,
    double eta,
    double phi) const {
    if (particleIndex >= legs_.size()) {
        throw std::out_of_range(
            "analytic Low-bound particle index is out of range");
    }
    if (!std::isfinite(eta) || !std::isfinite(phi)) {
        throw std::invalid_argument(
            "analytic Low-bound coordinates must be finite");
    }
    const AnalyticLowLegKinematics& leg = legs_[particleIndex];
    if (!leg.valid) {
        throw std::invalid_argument(
            "analytic Low-bound particle is not finite massive timelike input");
    }
    return positiveDenominatorScale(
        leg,
        eta - leg.rapidity,
        wrappedAngleDifference(phi, leg.phi));
}

double AnalyticLowBound::angularIntensity(
    double eta,
    double phi) const {
    if (!valid()) {
        throw std::invalid_argument(
            "analytic angular intensity requires finite massive timelike particles");
    }
    std::vector<double> denominators;
    denominators.reserve(legs_.size());
    for (std::size_t index = 0; index < legs_.size(); ++index) {
        denominators.push_back(denominatorScale(index, eta, phi));
    }

    double result = 0.0;
    std::size_t pairIndex = 0;
    for (std::size_t first = 0; first < legs_.size(); ++first) {
        for (std::size_t second = 0; second < legs_.size(); ++second) {
            result +=
                signedPairCoefficients_[pairIndex++] /
                (denominators[first] * denominators[second]);
        }
    }
    return result;
}

AnalyticLowAngularBound
AnalyticLowBound::globalAngularIntensityBound() const {
    if (!valid()) {
        return AnalyticLowAngularBound{
            std::numeric_limits<double>::infinity(),
            invalidParticleCount_,
            false,
        };
    }
    double result = 0.0;
    std::size_t pairIndex = 0;
    for (std::size_t first = 0; first < legs_.size(); ++first) {
        for (std::size_t second = 0; second < legs_.size(); ++second) {
            result +=
                absolutePairCoefficients_[pairIndex++] /
                (legs_[first].globalMinimumDenominatorScale *
                 legs_[second].globalMinimumDenominatorScale);
        }
    }
    return AnalyticLowAngularBound{
        result,
        0,
        std::isfinite(result) && result >= 0.0,
    };
}

AnalyticLowAngularBound
AnalyticLowBound::cellAngularIntensityBound(
    double etaMin,
    double etaMax,
    double phiMin,
    double phiMax) const {
    if (!std::isfinite(etaMin) || !std::isfinite(etaMax) ||
        !std::isfinite(phiMin) || !std::isfinite(phiMax) ||
        !(etaMin < etaMax) || !(phiMin < phiMax)) {
        throw std::invalid_argument(
            "analytic Low-bound cell must have finite ordered bounds");
    }
    if (!valid()) {
        return AnalyticLowAngularBound{
            std::numeric_limits<double>::infinity(),
            invalidParticleCount_,
            false,
        };
    }

    std::vector<double> minimumDenominators;
    minimumDenominators.reserve(legs_.size());
    for (const AnalyticLowLegKinematics& leg : legs_) {
        minimumDenominators.push_back(
            positiveDenominatorScale(
                leg,
                intervalDistance(leg.rapidity, etaMin, etaMax),
                periodicIntervalDistance(leg.phi, phiMin, phiMax)));
    }

    double result = 0.0;
    std::size_t pairIndex = 0;
    for (std::size_t first = 0; first < legs_.size(); ++first) {
        for (std::size_t second = 0; second < legs_.size(); ++second) {
            result +=
                absolutePairCoefficients_[pairIndex++] /
                (minimumDenominators[first] *
                 minimumDenominators[second]);
        }
    }
    return AnalyticLowAngularBound{
        result,
        0,
        std::isfinite(result) && result >= 0.0,
    };
}

double analyticLowWeightBound(
    const AnalyticLowAngularBound& angularBound,
    const PhotonPhaseSpace& phaseSpace,
    double etaPhiProposalDensity,
    double alpha) {
    if (!validatePhotonPhaseSpace(phaseSpace).valid()) {
        throw std::invalid_argument(
            "analytic Low weight bound requires valid photon phase space");
    }
    if (!std::isfinite(etaPhiProposalDensity) ||
        !(etaPhiProposalDensity > 0.0)) {
        throw std::invalid_argument(
            "analytic Low weight bound requires positive eta-phi proposal density");
    }
    if (!std::isfinite(alpha) || !(alpha > 0.0)) {
        throw std::invalid_argument(
            "analytic Low weight bound requires positive finite alpha");
    }
    if (!angularBound.valid || !std::isfinite(angularBound.value) ||
        angularBound.value < 0.0) {
        return std::numeric_limits<double>::infinity();
    }
    const double normalization =
        alpha /
        std::pow(2.0 * std::acos(-1.0), 2) *
        std::log(phaseSpace.kTMax / phaseSpace.kTMin);
    return normalization *
           angularBound.value /
           etaPhiProposalDensity;
}

}  // namespace lowmc
