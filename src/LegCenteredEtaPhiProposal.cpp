#include "lowmc/LegCenteredEtaPhiProposal.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace lowmc {
namespace {

constexpr double BetaSeriesThreshold = 1.0e-8;

double twoPi() {
    return 2.0 * std::acos(-1.0);
}

double clampOpenCosine(double value) {
    const double lower = std::nextafter(-1.0, 0.0);
    const double upper = std::nextafter(1.0, 0.0);
    return std::max(lower, std::min(value, upper));
}

double componentSphereNormalization(double beta) {
    if (beta < BetaSeriesThreshold) {
        return 1.0 / (4.0 * std::acos(-1.0));
    }
    const double logarithm = std::log1p(beta) - std::log1p(-beta);
    return beta / (twoPi() * logarithm);
}

double sampleComponentCosine(double beta, double unitPolar) {
    if (beta < BetaSeriesThreshold) {
        return clampOpenCosine(2.0 * unitPolar - 1.0);
    }
    const double logarithm = std::log1p(beta) - std::log1p(-beta);
    const double cosine =
        (1.0 - (1.0 + beta) * std::exp(-unitPolar * logarithm)) / beta;
    return clampOpenCosine(cosine);
}

double wrapPhiFromPhaseMinimum(double phi, const PhotonPhaseSpace& phaseSpace) {
    const double period = twoPi();
    double wrapped = std::fmod(phi - phaseSpace.phiMin, period);
    if (wrapped < 0.0) {
        wrapped += period;
    }
    return phaseSpace.phiMin + wrapped;
}

void requireUnitInterval(double value, const char* label) {
    if (!std::isfinite(value) || value < 0.0 || value >= 1.0) {
        throw std::invalid_argument(std::string(label) +
                                    " unit random value must satisfy 0 <= u < 1");
    }
}

}  // namespace

void validateLegCenteredEtaPhiProposalConfig(
    const LegCenteredEtaPhiProposalConfig& config) {
    if (!std::isfinite(config.uniformFraction) ||
        config.uniformFraction <= 0.0 || config.uniformFraction >= 1.0) {
        throw std::invalid_argument(
            "leg-centred uniformFraction must satisfy 0 < f < 1");
    }
    if (!std::isfinite(config.minimumOneMinusBeta) ||
        config.minimumOneMinusBeta <= 0.0 || config.minimumOneMinusBeta >= 1.0) {
        throw std::invalid_argument(
            "leg-centred minimumOneMinusBeta must satisfy 0 < value < 1");
    }
}

LegCenteredEtaPhiProposal makeLegCenteredEtaPhiProposal(
    const PhotonPhaseSpace& phaseSpace,
    const std::vector<ChargedParticle>& particles,
    const LegCenteredEtaPhiProposalConfig& config) {
    if (!validatePhotonPhaseSpace(phaseSpace).valid()) {
        throw std::invalid_argument("invalid photon phase space");
    }
    validateLegCenteredEtaPhiProposalConfig(config);
    if (phaseSpace.phiMax - phaseSpace.phiMin >
        twoPi() + 32.0 * std::numeric_limits<double>::epsilon()) {
        throw std::invalid_argument(
            "leg-centred proposal requires a phi interval no wider than 2 pi");
    }

    std::vector<LegCenteredEtaPhiProposalComponent> components;
    double totalWeight = 0.0;
    for (const ChargedParticle& particle : particles) {
        const double energy = std::abs(particle.momentum.E);
        const double momentum = std::sqrt(
            particle.momentum.px * particle.momentum.px +
            particle.momentum.py * particle.momentum.py +
            particle.momentum.pz * particle.momentum.pz);
        const double componentWeight = std::abs(particle.charge);
        if (!std::isfinite(energy) || !std::isfinite(momentum) ||
            energy <= 0.0 || momentum <= 0.0 ||
            !std::isfinite(componentWeight) || componentWeight <= 0.0) {
            continue;
        }

        const double maximumBeta = 1.0 - config.minimumOneMinusBeta;
        const double beta = std::min(momentum / energy, maximumBeta);
        if (!std::isfinite(beta) || beta < 0.0) {
            continue;
        }
        components.push_back(LegCenteredEtaPhiProposalComponent{
            particle.momentum.px / momentum,
            particle.momentum.py / momentum,
            particle.momentum.pz / momentum,
            beta,
            componentSphereNormalization(beta),
            componentWeight,
        });
        totalWeight += componentWeight;
    }

    if (components.empty() || totalWeight <= 0.0) {
        throw std::invalid_argument("leg-centred proposal has no active components");
    }

    const double componentFraction = 1.0 - config.uniformFraction;
    for (LegCenteredEtaPhiProposalComponent& component : components) {
        component.mixtureWeight *= componentFraction / totalWeight;
    }

    return LegCenteredEtaPhiProposal{phaseSpace, config, components};
}

double legCenteredEtaPhiDensity(
    const LegCenteredEtaPhiProposal& proposal,
    double eta,
    double phi) {
    validateLegCenteredEtaPhiProposalConfig(proposal.config);
    if (!std::isfinite(eta) || !std::isfinite(phi)) {
        throw std::invalid_argument("leg-centred density requires finite eta and phi");
    }

    const double coshEta = std::cosh(eta);
    if (!std::isfinite(coshEta)) {
        return 0.0;
    }
    const double inverseCosh = 1.0 / coshEta;
    const double nx = std::cos(phi) * inverseCosh;
    const double ny = std::sin(phi) * inverseCosh;
    const double nz = std::tanh(eta);

    double sphereDensity = proposal.config.uniformFraction /
                           (4.0 * std::acos(-1.0));
    for (const LegCenteredEtaPhiProposalComponent& component : proposal.components) {
        const double cosine = component.nx * nx + component.ny * ny + component.nz * nz;
        const double denominator = 1.0 - component.beta * cosine;
        sphereDensity += component.mixtureWeight *
                         component.sphereNormalization / denominator;
    }

    return sphereDensity * inverseCosh * inverseCosh;
}

LegCenteredEtaPhiProposalPoint sampleLegCenteredEtaPhi(
    const LegCenteredEtaPhiProposal& proposal,
    double unitComponent,
    double unitPolar,
    double unitAzimuth) {
    validateLegCenteredEtaPhiProposalConfig(proposal.config);
    requireUnitInterval(unitComponent, "component");
    requireUnitInterval(unitPolar, "polar");
    requireUnitInterval(unitAzimuth, "azimuth");
    if (proposal.components.empty()) {
        throw std::invalid_argument("leg-centred proposal has no active components");
    }

    double nx = 0.0;
    double ny = 0.0;
    double nz = 0.0;
    if (unitComponent < proposal.config.uniformFraction) {
        nz = clampOpenCosine(2.0 * unitPolar - 1.0);
        const double transverse = std::sqrt(std::max(0.0, 1.0 - nz * nz));
        const double azimuth = twoPi() * unitAzimuth;
        nx = transverse * std::cos(azimuth);
        ny = transverse * std::sin(azimuth);
    } else {
        const double target = unitComponent - proposal.config.uniformFraction;
        double cumulative = 0.0;
        const LegCenteredEtaPhiProposalComponent* selected = &proposal.components.back();
        for (const LegCenteredEtaPhiProposalComponent& component : proposal.components) {
            cumulative += component.mixtureWeight;
            if (target < cumulative) {
                selected = &component;
                break;
            }
        }

        const double cosine = sampleComponentCosine(selected->beta, unitPolar);
        const double sine = std::sqrt(std::max(0.0, 1.0 - cosine * cosine));
        const double azimuth = twoPi() * unitAzimuth;

        double ux = 0.0;
        double uy = 0.0;
        double uz = 0.0;
        if (std::abs(selected->nz) < 0.9) {
            const double norm = std::hypot(selected->nx, selected->ny);
            ux = -selected->ny / norm;
            uy = selected->nx / norm;
        } else {
            const double norm = std::hypot(selected->ny, selected->nz);
            uy = -selected->nz / norm;
            uz = selected->ny / norm;
        }
        const double vx = selected->ny * uz - selected->nz * uy;
        const double vy = selected->nz * ux - selected->nx * uz;
        const double vz = selected->nx * uy - selected->ny * ux;

        nx = cosine * selected->nx +
             sine * (std::cos(azimuth) * ux + std::sin(azimuth) * vx);
        ny = cosine * selected->ny +
             sine * (std::cos(azimuth) * uy + std::sin(azimuth) * vy);
        nz = cosine * selected->nz +
             sine * (std::cos(azimuth) * uz + std::sin(azimuth) * vz);
        nz = clampOpenCosine(nz);
    }

    const double transverse = std::hypot(nx, ny);
    const double eta = std::asinh(nz / transverse);
    const double phi = wrapPhiFromPhaseMinimum(std::atan2(ny, nx), proposal.phaseSpace);
    const bool insideAcceptance =
        eta >= proposal.phaseSpace.etaMin && eta < proposal.phaseSpace.etaMax &&
        phi >= proposal.phaseSpace.phiMin && phi < proposal.phaseSpace.phiMax;

    return LegCenteredEtaPhiProposalPoint{
        eta,
        phi,
        legCenteredEtaPhiDensity(proposal, eta, phi),
        insideAcceptance,
    };
}

}  // namespace lowmc
