#include "lowmc/Dielectron.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>

namespace lowmc {
namespace {

double pi() {
    return std::acos(-1.0);
}

void validateUnit(double value, const char* name) {
    if (!std::isfinite(value) || value < 0.0 || value >= 1.0) {
        throw std::invalid_argument(std::string(name) + " must be in [0, 1)");
    }
}

void validateConstants(double electronMass, double alpha) {
    if (!std::isfinite(electronMass) || electronMass <= 0.0) {
        throw std::invalid_argument("electron mass must be positive and finite");
    }
    if (!std::isfinite(alpha) || alpha <= 0.0) {
        throw std::invalid_argument("fine-structure constant must be positive and finite");
    }
}

FourVector boostFromRestFrame(const FourVector& rest, const FourVector& parent) {
    const double inverseEnergy = 1.0 / parent.E;
    const double bx = parent.px * inverseEnergy;
    const double by = parent.py * inverseEnergy;
    const double bz = parent.pz * inverseEnergy;
    const double beta2 = bx * bx + by * by + bz * bz;
    if (!(beta2 > 0.0)) {
        return rest;
    }

    const double gamma = parent.E / std::sqrt(parent.mass2());
    const double betaDotP = bx * rest.px + by * rest.py + bz * rest.pz;
    const double spatialScale =
        ((gamma - 1.0) * betaDotP / beta2) + gamma * rest.E;
    return FourVector{
        gamma * (rest.E + betaDotP),
        rest.px + spatialScale * bx,
        rest.py + spatialScale * by,
        rest.pz + spatialScale * bz,
    };
}

FourVector boostToRestFrame(const FourVector& vector, const FourVector& parent) {
    const double inverseEnergy = 1.0 / parent.E;
    const double bx = parent.px * inverseEnergy;
    const double by = parent.py * inverseEnergy;
    const double bz = parent.pz * inverseEnergy;
    const double beta2 = bx * bx + by * by + bz * bz;
    if (!(beta2 > 0.0)) {
        return vector;
    }

    const double gamma = parent.E / std::sqrt(parent.mass2());
    const double betaDotP = bx * vector.px + by * vector.py + bz * vector.pz;
    const double spatialScale =
        ((gamma - 1.0) * betaDotP / beta2) - gamma * vector.E;
    return FourVector{
        gamma * (vector.E - betaDotP),
        vector.px + spatialScale * bx,
        vector.py + spatialScale * by,
        vector.pz + spatialScale * bz,
    };
}

}  // namespace

DielectronPhaseSpaceValidation validateDielectronPhaseSpace(
    const DielectronPhaseSpace& phaseSpace,
    double electronMass) {
    if (!std::isfinite(electronMass) || electronMass <= 0.0 ||
        !std::isfinite(phaseSpace.massMin) ||
        !std::isfinite(phaseSpace.massMax) ||
        !std::isfinite(phaseSpace.qTMin) ||
        !std::isfinite(phaseSpace.qTMax) ||
        !std::isfinite(phaseSpace.rapidityMin) ||
        !std::isfinite(phaseSpace.rapidityMax) ||
        !std::isfinite(phaseSpace.phiMin) ||
        !std::isfinite(phaseSpace.phiMax)) {
        return {DielectronPhaseSpaceValidationError::NonFiniteBound};
    }
    if (phaseSpace.massMin < 2.0 * electronMass) {
        return {DielectronPhaseSpaceValidationError::BelowPairThreshold};
    }
    if (phaseSpace.massMax <= phaseSpace.massMin) {
        return {DielectronPhaseSpaceValidationError::EmptyMassRange};
    }
    if (phaseSpace.qTMin <= 0.0) {
        return {DielectronPhaseSpaceValidationError::NonPositiveQTMin};
    }
    if (phaseSpace.qTMax <= phaseSpace.qTMin) {
        return {DielectronPhaseSpaceValidationError::EmptyQTRange};
    }
    if (phaseSpace.rapidityMax <= phaseSpace.rapidityMin) {
        return {DielectronPhaseSpaceValidationError::EmptyRapidityRange};
    }
    if (phaseSpace.phiMax <= phaseSpace.phiMin) {
        return {DielectronPhaseSpaceValidationError::EmptyPhiRange};
    }
    return {};
}

FourVector virtualPhotonMomentumFromMassQTRapidityPhi(
    double mass,
    double qT,
    double rapidity,
    double phi) {
    if (!std::isfinite(mass) || mass <= 0.0 || !std::isfinite(qT) || qT < 0.0 ||
        !std::isfinite(rapidity) || !std::isfinite(phi)) {
        throw std::invalid_argument("invalid virtual-photon kinematics");
    }
    const double transverseMass = std::hypot(mass, qT);
    const FourVector momentum{
        transverseMass * std::cosh(rapidity),
        qT * std::cos(phi),
        qT * std::sin(phi),
        transverseMass * std::sinh(rapidity),
    };
    if (!std::isfinite(momentum.E) || !std::isfinite(momentum.px) ||
        !std::isfinite(momentum.py) || !std::isfinite(momentum.pz) ||
        !std::isfinite(momentum.mass2()) || momentum.mass2() <= 0.0) {
        throw std::invalid_argument("virtual-photon kinematics overflowed");
    }
    return momentum;
}

double krollWadaConversionDensityDMass2(
    double mass,
    double electronMass,
    double alpha) {
    validateConstants(electronMass, alpha);
    if (!std::isfinite(mass) || mass < 2.0 * electronMass) {
        throw std::invalid_argument("pair mass must be finite and at or above threshold");
    }
    const double mass2 = mass * mass;
    const double massRatio = electronMass * electronMass / mass2;
    const double beta = std::sqrt(std::max(0.0, 1.0 - 4.0 * massRatio));
    return alpha / (3.0 * pi() * mass2) * beta * (1.0 + 2.0 * massRatio);
}

double lowVirtualPhotonDensity(
    const LowCurrent& current,
    const FourVector& virtualPhotonMomentum,
    double alpha) {
    if (!std::isfinite(virtualPhotonMomentum.mass2()) ||
        virtualPhotonMomentum.mass2() <= 0.0) {
        throw std::invalid_argument("virtual-photon momentum must be timelike");
    }
    if (!std::isfinite(alpha) || alpha <= 0.0) {
        throw std::invalid_argument("fine-structure constant must be positive and finite");
    }
    const double twoPi = 2.0 * pi();
    return alpha / (twoPi * twoPi) * current.intensity(virtualPhotonMomentum);
}

WeightedDielectronSample makeWeightedDielectronSample(
    const LowCurrent& current,
    const DielectronPhaseSpace& phaseSpace,
    double unitMass,
    double unitQT,
    double unitRapidity,
    double unitPhi,
    double unitDecayCosTheta,
    double unitDecayPhi,
    double electronMass,
    double alpha) {
    validateConstants(electronMass, alpha);
    if (!validateDielectronPhaseSpace(phaseSpace, electronMass).valid()) {
        throw std::invalid_argument("invalid dielectron phase space");
    }
    validateUnit(unitMass, "mass unit value");
    validateUnit(unitQT, "qT unit value");
    validateUnit(unitRapidity, "rapidity unit value");
    validateUnit(unitPhi, "phi unit value");
    validateUnit(unitDecayCosTheta, "decay cos(theta) unit value");
    validateUnit(unitDecayPhi, "decay phi unit value");

    const double logMassRange = std::log(phaseSpace.massMax / phaseSpace.massMin);
    const double mass = phaseSpace.massMin * std::exp(unitMass * logMassRange);
    const double logQTRange = std::log(phaseSpace.qTMax / phaseSpace.qTMin);
    const double qT = phaseSpace.qTMin * std::exp(unitQT * logQTRange);
    const double rapidity = phaseSpace.rapidityMin +
        unitRapidity * (phaseSpace.rapidityMax - phaseSpace.rapidityMin);
    const double phi = phaseSpace.phiMin +
        unitPhi * (phaseSpace.phiMax - phaseSpace.phiMin);
    const FourVector virtualPhoton =
        virtualPhotonMomentumFromMassQTRapidityPhi(mass, qT, rapidity, phi);

    const double massProposalDensity = 1.0 / (2.0 * mass * mass * logMassRange);
    const double qTProposalDensity = 1.0 / (qT * qT * logQTRange);
    const double rapidityPhiProposalDensity =
        1.0 / ((phaseSpace.rapidityMax - phaseSpace.rapidityMin) *
               (phaseSpace.phiMax - phaseSpace.phiMin));
    const double proposalDensity =
        massProposalDensity * qTProposalDensity * rapidityPhiProposalDensity;

    const double leptonEnergy = 0.5 * mass;
    const double leptonMomentum =
        0.5 * std::sqrt(std::max(0.0, mass * mass - 4.0 * electronMass * electronMass));
    const double cosTheta = 2.0 * unitDecayCosTheta - 1.0;
    const double sinTheta = std::sqrt(std::max(0.0, 1.0 - cosTheta * cosTheta));
    const double decayPhi = 2.0 * pi() * unitDecayPhi;
    const double nx = sinTheta * std::cos(decayPhi);
    const double ny = sinTheta * std::sin(decayPhi);
    const double nz = cosTheta;
    const FourVector electronRest{
        leptonEnergy,
        leptonMomentum * nx,
        leptonMomentum * ny,
        leptonMomentum * nz,
    };
    const FourVector positronRest{
        leptonEnergy,
        -leptonMomentum * nx,
        -leptonMomentum * ny,
        -leptonMomentum * nz,
    };

    const FourVector currentRest = boostToRestFrame(current.current(virtualPhoton), virtualPhoton);
    const double currentSpatial2 = currentRest.px * currentRest.px +
        currentRest.py * currentRest.py + currentRest.pz * currentRest.pz;
    const double beta2 = std::max(0.0, 1.0 -
        4.0 * electronMass * electronMass / (mass * mass));
    double decayAngularWeight = 1.0;
    if (currentSpatial2 > 0.0) {
        const double projection =
            (currentRest.px * nx + currentRest.py * ny + currentRest.pz * nz);
        const double cos2ToCurrent = projection * projection / currentSpatial2;
        decayAngularWeight =
            (1.0 - beta2 * cos2ToCurrent) / (1.0 - beta2 / 3.0);
    }

    const double productionDensity =
        lowVirtualPhotonDensity(current, virtualPhoton, alpha);
    const double conversionDensity =
        krollWadaConversionDensityDMass2(mass, electronMass, alpha);
    const double weight =
        productionDensity * conversionDensity * decayAngularWeight / proposalDensity;

    return WeightedDielectronSample{
        Dielectron{
            boostFromRestFrame(electronRest, virtualPhoton),
            boostFromRestFrame(positronRest, virtualPhoton),
            weight,
        },
        virtualPhoton,
        mass,
        qT,
        rapidity,
        phi,
        proposalDensity,
        productionDensity,
        conversionDensity,
        decayAngularWeight,
    };
}

WeightedDielectronGenerationResult generateWeightedDielectrons(
    const std::vector<ChargedParticle>& selectedParticles,
    const WeightedDielectronGeneratorConfig& config,
    double electronMass,
    double alpha) {
    validateConstants(electronMass, alpha);
    if (!validateDielectronPhaseSpace(config.phaseSpace, electronMass).valid()) {
        throw std::invalid_argument("invalid dielectron phase space");
    }
    if (config.samples < 0) {
        throw std::invalid_argument("dielectron sample count must be non-negative");
    }
    if (config.oversamplingFactor <= 0) {
        throw std::invalid_argument("dielectron oversampling factor must be positive");
    }
    if (config.samples > 0 &&
        config.samples > std::numeric_limits<int>::max() / config.oversamplingFactor) {
        throw std::invalid_argument("dielectron oversampled sample count is too large");
    }

    const int totalSamples = config.samples * config.oversamplingFactor;
    const double inverseOversampling = 1.0 / static_cast<double>(config.oversamplingFactor);
    const LowCurrent current(selectedParticles);
    WeightedDielectronGenerationResult result;
    result.chargeConservation = current.chargeConservation();
    result.samples.reserve(static_cast<std::size_t>(totalSamples));

    std::mt19937_64 rng(config.seed);
    std::vector<double> weights;
    weights.reserve(static_cast<std::size_t>(totalSamples));
    for (int index = 0; index < totalSamples; ++index) {
        auto unit = [&rng]() { return std::generate_canonical<double, 53>(rng); };
        const double unitMass = unit();
        const double unitQT = unit();
        const double unitRapidity = unit();
        const double unitPhi = unit();
        const double unitDecayCosTheta = unit();
        const double unitDecayPhi = unit();
        WeightedDielectronSample sample = makeWeightedDielectronSample(
            current,
            config.phaseSpace,
            unitMass,
            unitQT,
            unitRapidity,
            unitPhi,
            unitDecayCosTheta,
            unitDecayPhi,
            electronMass,
            alpha);
        sample.dielectron.weight *= inverseOversampling;
        weights.push_back(sample.dielectron.weight);
        result.samples.push_back(sample);
    }
    result.weightSummary = summarizeWeights(weights);
    return result;
}

}  // namespace lowmc
