#pragma once

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "lowmc/ChargedParticle.h"
#include "lowmc/PhotonPhaseSpace.h"

namespace lowmc {

struct UniformEtaPhiProposalPoint {
    double eta = 0.0;
    double phi = 0.0;
    double density = 0.0;
};

struct AdaptiveEtaPhiProposalConfig {
    double uniformFraction = 0.2;
    double widthScale = 5.0;
    double minHalfWidth = 0.02;
    double maxHalfWidth = 0.8;
};

struct AdaptiveEtaPhiProposalComponent {
    double etaCenter = 0.0;
    double phiCenter = 0.0;
    double etaMin = 0.0;
    double etaMax = 0.0;
    double phiMin = 0.0;
    double phiMax = 0.0;
    double mixtureWeight = 0.0;
};

struct AdaptiveEtaPhiProposal {
    PhotonPhaseSpace phaseSpace;
    AdaptiveEtaPhiProposalConfig config;
    std::vector<AdaptiveEtaPhiProposalComponent> components;
};

struct AdaptiveEtaPhiProposalPoint {
    double eta = 0.0;
    double phi = 0.0;
    double density = 0.0;
};

struct SmoothEtaPhiProposalConfig {
    double uniformFraction = 0.2;
    double widthScale = 5.0;
    double minHalfWidth = 0.02;
    double maxHalfWidth = 0.8;
};

struct SmoothEtaPhiProposalComponent {
    double etaCenter = 0.0;
    double phiCenter = 0.0;
    double etaHalfWidth = 0.0;
    double phiHalfWidth = 0.0;
    double etaLowerOffset = 0.0;
    double etaUpperOffset = 0.0;
    double phiLowerOffset = 0.0;
    double phiUpperOffset = 0.0;
    double etaNormalization = 0.0;
    double phiNormalization = 0.0;
    double mixtureWeight = 0.0;
};

struct SmoothEtaPhiProposal {
    PhotonPhaseSpace phaseSpace;
    SmoothEtaPhiProposalConfig config;
    std::vector<SmoothEtaPhiProposalComponent> components;
};

struct SmoothEtaPhiProposalPoint {
    double eta = 0.0;
    double phi = 0.0;
    double density = 0.0;
};

inline void requireValidEtaPhiProposalRange(const PhotonPhaseSpace& phaseSpace) {
    if (!validatePhotonPhaseSpace(phaseSpace).valid()) {
        throw std::invalid_argument("invalid photon phase space");
    }
}

[[nodiscard]] inline double uniformEtaPhiDensity(const PhotonPhaseSpace& phaseSpace) {
    requireValidEtaPhiProposalRange(phaseSpace);
    return 1.0 /
           ((phaseSpace.etaMax - phaseSpace.etaMin) *
            (phaseSpace.phiMax - phaseSpace.phiMin));
}

[[nodiscard]] inline UniformEtaPhiProposalPoint sampleUniformEtaPhi(
    const PhotonPhaseSpace& phaseSpace,
    double unitEta,
    double unitPhi) {
    requireValidEtaPhiProposalRange(phaseSpace);

    if (!std::isfinite(unitEta) || unitEta < 0.0 || unitEta >= 1.0) {
        throw std::invalid_argument("eta unit random value must satisfy 0 <= u < 1");
    }
    if (!std::isfinite(unitPhi) || unitPhi < 0.0 || unitPhi >= 1.0) {
        throw std::invalid_argument("phi unit random value must satisfy 0 <= u < 1");
    }

    double eta = phaseSpace.etaMin + unitEta * (phaseSpace.etaMax - phaseSpace.etaMin);
    double phi = phaseSpace.phiMin + unitPhi * (phaseSpace.phiMax - phaseSpace.phiMin);
    if (eta >= phaseSpace.etaMax) {
        eta = std::nextafter(phaseSpace.etaMax, phaseSpace.etaMin);
    }
    if (phi >= phaseSpace.phiMax) {
        phi = std::nextafter(phaseSpace.phiMax, phaseSpace.phiMin);
    }

    return UniformEtaPhiProposalPoint{
        eta,
        phi,
        uniformEtaPhiDensity(phaseSpace),
    };
}

[[nodiscard]] inline double clampDouble(double value, double lower, double upper) {
    return std::max(lower, std::min(value, upper));
}

[[nodiscard]] inline double particlePseudorapidity(const ChargedParticle& particle) {
    const double pt = particle.momentum.pt();
    if (pt <= 0.0 || !std::isfinite(pt)) {
        throw std::invalid_argument("adaptive eta-phi proposal requires finite nonzero pT");
    }
    return std::asinh(particle.momentum.pz / pt);
}

[[nodiscard]] inline double particlePhiInRange(
    const ChargedParticle& particle,
    const PhotonPhaseSpace& phaseSpace) {
    const double twoPi = 2.0 * std::acos(-1.0);
    double phi = std::atan2(particle.momentum.py, particle.momentum.px);
    while (phi < phaseSpace.phiMin) {
        phi += twoPi;
    }
    while (phi >= phaseSpace.phiMax) {
        phi -= twoPi;
    }
    return phi;
}

[[nodiscard]] inline double particleAngularHalfWidth(
    const ChargedParticle& particle,
    const AdaptiveEtaPhiProposalConfig& config) {
    const double mass2 = std::max(0.0, particle.momentum.mass2());
    const double energy = std::abs(particle.momentum.E);
    if (energy <= 0.0 || !std::isfinite(energy)) {
        throw std::invalid_argument("adaptive eta-phi proposal requires finite positive energy");
    }

    const double massOverEnergy = std::sqrt(mass2) / energy;
    return clampDouble(
        config.widthScale * massOverEnergy,
        config.minHalfWidth,
        config.maxHalfWidth);
}

[[nodiscard]] inline double particleSmoothAngularHalfWidth(
    const ChargedParticle& particle,
    const SmoothEtaPhiProposalConfig& config) {
    const double mass2 = std::max(0.0, particle.momentum.mass2());
    const double energy = std::abs(particle.momentum.E);
    if (energy <= 0.0 || !std::isfinite(energy)) {
        throw std::invalid_argument("smooth eta-phi proposal requires finite positive energy");
    }

    const double massOverEnergy = std::sqrt(mass2) / energy;
    return clampDouble(
        config.widthScale * massOverEnergy,
        config.minHalfWidth,
        config.maxHalfWidth);
}

inline void validateAdaptiveEtaPhiProposalConfig(
    const AdaptiveEtaPhiProposalConfig& config) {
    if (!std::isfinite(config.uniformFraction) ||
        config.uniformFraction < 0.0 ||
        config.uniformFraction >= 1.0) {
        throw std::invalid_argument("adaptive uniformFraction must satisfy 0 <= f < 1");
    }
    if (!std::isfinite(config.widthScale) || config.widthScale <= 0.0) {
        throw std::invalid_argument("adaptive widthScale must be positive");
    }
    if (!std::isfinite(config.minHalfWidth) || config.minHalfWidth <= 0.0) {
        throw std::invalid_argument("adaptive minHalfWidth must be positive");
    }
    if (!std::isfinite(config.maxHalfWidth) || config.maxHalfWidth < config.minHalfWidth) {
        throw std::invalid_argument("adaptive maxHalfWidth must be at least minHalfWidth");
    }
}

inline void validateSmoothEtaPhiProposalConfig(
    const SmoothEtaPhiProposalConfig& config) {
    if (!std::isfinite(config.uniformFraction) ||
        config.uniformFraction < 0.0 ||
        config.uniformFraction >= 1.0) {
        throw std::invalid_argument("smooth uniformFraction must satisfy 0 <= f < 1");
    }
    if (!std::isfinite(config.widthScale) || config.widthScale <= 0.0) {
        throw std::invalid_argument("smooth widthScale must be positive");
    }
    if (!std::isfinite(config.minHalfWidth) || config.minHalfWidth <= 0.0) {
        throw std::invalid_argument("smooth minHalfWidth must be positive");
    }
    if (!std::isfinite(config.maxHalfWidth) || config.maxHalfWidth < config.minHalfWidth) {
        throw std::invalid_argument("smooth maxHalfWidth must be at least minHalfWidth");
    }
}

[[nodiscard]] inline double triangularKernelCdf(double offset, double halfWidth) {
    const double x = clampDouble(offset / halfWidth, -1.0, 1.0);
    if (x < 0.0) {
        return 0.5 + x + 0.5 * x * x;
    }
    return 0.5 + x - 0.5 * x * x;
}

[[nodiscard]] inline double triangularKernelPdf(double offset, double halfWidth) {
    const double scaled = std::abs(offset) / halfWidth;
    if (scaled > 1.0) {
        return 0.0;
    }
    return (1.0 - scaled) / halfWidth;
}

[[nodiscard]] inline double inverseTriangularKernelCdf(double cumulative, double halfWidth) {
    const double y = clampDouble(cumulative, 0.0, 1.0);
    if (y < 0.5) {
        return halfWidth * (-1.0 + std::sqrt(2.0 * y));
    }
    return halfWidth * (1.0 - std::sqrt(2.0 * (1.0 - y)));
}

[[nodiscard]] inline double sampleTruncatedTriangularKernel(
    double unit,
    double lowerOffset,
    double upperOffset,
    double halfWidth) {
    if (!std::isfinite(unit) || unit < 0.0 || unit >= 1.0) {
        throw std::invalid_argument("triangular-kernel unit random value must satisfy 0 <= u < 1");
    }
    const double lowerCdf = triangularKernelCdf(lowerOffset, halfWidth);
    const double upperCdf = triangularKernelCdf(upperOffset, halfWidth);
    const double cumulative = lowerCdf + unit * (upperCdf - lowerCdf);
    return inverseTriangularKernelCdf(cumulative, halfWidth);
}

[[nodiscard]] inline AdaptiveEtaPhiProposal makeAdaptiveEtaPhiProposal(
    const PhotonPhaseSpace& phaseSpace,
    const std::vector<ChargedParticle>& particles,
    const AdaptiveEtaPhiProposalConfig& config = {}) {
    requireValidEtaPhiProposalRange(phaseSpace);
    validateAdaptiveEtaPhiProposalConfig(config);

    std::vector<AdaptiveEtaPhiProposalComponent> components;
    double totalComponentWeight = 0.0;

    for (const ChargedParticle& particle : particles) {
        const double etaCenter = particlePseudorapidity(particle);
        const double phiCenter = particlePhiInRange(particle, phaseSpace);
        const double halfWidth = particleAngularHalfWidth(particle, config);

        const double etaMin = std::max(phaseSpace.etaMin, etaCenter - halfWidth);
        const double etaMax = std::min(phaseSpace.etaMax, etaCenter + halfWidth);
        const double phiMin = std::max(phaseSpace.phiMin, phiCenter - halfWidth);
        const double phiMax = std::min(phaseSpace.phiMax, phiCenter + halfWidth);
        if (etaMin >= etaMax || phiMin >= phiMax) {
            continue;
        }

        const double componentWeight = std::abs(particle.charge);
        if (componentWeight <= 0.0 || !std::isfinite(componentWeight)) {
            continue;
        }
        components.push_back(AdaptiveEtaPhiProposalComponent{
            etaCenter,
            phiCenter,
            etaMin,
            etaMax,
            phiMin,
            phiMax,
            componentWeight,
        });
        totalComponentWeight += componentWeight;
    }

    if (components.empty() || totalComponentWeight <= 0.0) {
        throw std::invalid_argument("adaptive eta-phi proposal has no active components");
    }

    const double componentFraction = 1.0 - config.uniformFraction;
    for (AdaptiveEtaPhiProposalComponent& component : components) {
        component.mixtureWeight = componentFraction * component.mixtureWeight /
                                  totalComponentWeight;
    }

    return AdaptiveEtaPhiProposal{phaseSpace, config, components};
}

[[nodiscard]] inline double adaptiveEtaPhiDensity(
    const AdaptiveEtaPhiProposal& proposal,
    double eta,
    double phi) {
    requireValidEtaPhiProposalRange(proposal.phaseSpace);
    validateAdaptiveEtaPhiProposalConfig(proposal.config);

    if (!std::isfinite(eta) || !std::isfinite(phi) ||
        eta < proposal.phaseSpace.etaMin || eta > proposal.phaseSpace.etaMax ||
        phi < proposal.phaseSpace.phiMin || phi > proposal.phaseSpace.phiMax) {
        throw std::invalid_argument("eta-phi point is outside proposal phase space");
    }

    double density = proposal.config.uniformFraction * uniformEtaPhiDensity(proposal.phaseSpace);
    for (const AdaptiveEtaPhiProposalComponent& component : proposal.components) {
        if (eta >= component.etaMin && eta <= component.etaMax &&
            phi >= component.phiMin && phi <= component.phiMax) {
            const double area =
                (component.etaMax - component.etaMin) *
                (component.phiMax - component.phiMin);
            density += component.mixtureWeight / area;
        }
    }

    return density;
}

[[nodiscard]] inline AdaptiveEtaPhiProposalPoint sampleAdaptiveEtaPhi(
    const AdaptiveEtaPhiProposal& proposal,
    double unitComponent,
    double unitEta,
    double unitPhi) {
    requireValidEtaPhiProposalRange(proposal.phaseSpace);
    validateAdaptiveEtaPhiProposalConfig(proposal.config);

    if (!std::isfinite(unitComponent) || unitComponent < 0.0 || unitComponent >= 1.0) {
        throw std::invalid_argument("component unit random value must satisfy 0 <= u < 1");
    }
    if (!std::isfinite(unitEta) || unitEta < 0.0 || unitEta >= 1.0) {
        throw std::invalid_argument("eta unit random value must satisfy 0 <= u < 1");
    }
    if (!std::isfinite(unitPhi) || unitPhi < 0.0 || unitPhi >= 1.0) {
        throw std::invalid_argument("phi unit random value must satisfy 0 <= u < 1");
    }

    double eta = 0.0;
    double phi = 0.0;
    if (unitComponent < proposal.config.uniformFraction) {
        const UniformEtaPhiProposalPoint uniform =
            sampleUniformEtaPhi(proposal.phaseSpace, unitEta, unitPhi);
        eta = uniform.eta;
        phi = uniform.phi;
    } else {
        const double target = unitComponent - proposal.config.uniformFraction;
        double cumulative = 0.0;
        const AdaptiveEtaPhiProposalComponent* selected = &proposal.components.back();
        for (const AdaptiveEtaPhiProposalComponent& component : proposal.components) {
            cumulative += component.mixtureWeight;
            if (target < cumulative) {
                selected = &component;
                break;
            }
        }

        eta = selected->etaMin + unitEta * (selected->etaMax - selected->etaMin);
        phi = selected->phiMin + unitPhi * (selected->phiMax - selected->phiMin);
        if (eta >= selected->etaMax) {
            eta = std::nextafter(selected->etaMax, selected->etaMin);
        }
        if (phi >= selected->phiMax) {
            phi = std::nextafter(selected->phiMax, selected->phiMin);
        }
    }

    return AdaptiveEtaPhiProposalPoint{
        eta,
        phi,
        adaptiveEtaPhiDensity(proposal, eta, phi),
    };
}

[[nodiscard]] inline SmoothEtaPhiProposal makeSmoothEtaPhiProposal(
    const PhotonPhaseSpace& phaseSpace,
    const std::vector<ChargedParticle>& particles,
    const SmoothEtaPhiProposalConfig& config = {}) {
    requireValidEtaPhiProposalRange(phaseSpace);
    validateSmoothEtaPhiProposalConfig(config);

    std::vector<SmoothEtaPhiProposalComponent> components;
    double totalComponentWeight = 0.0;

    for (const ChargedParticle& particle : particles) {
        const double etaCenter = particlePseudorapidity(particle);
        const double phiCenter = particlePhiInRange(particle, phaseSpace);
        const double halfWidth = particleSmoothAngularHalfWidth(particle, config);

        const double etaLowerOffset =
            std::max(-halfWidth, phaseSpace.etaMin - etaCenter);
        const double etaUpperOffset =
            std::min(halfWidth, phaseSpace.etaMax - etaCenter);
        const double phiLowerOffset =
            std::max(-halfWidth, phaseSpace.phiMin - phiCenter);
        const double phiUpperOffset =
            std::min(halfWidth, phaseSpace.phiMax - phiCenter);
        if (etaLowerOffset >= etaUpperOffset || phiLowerOffset >= phiUpperOffset) {
            continue;
        }

        const double etaNormalization =
            triangularKernelCdf(etaUpperOffset, halfWidth) -
            triangularKernelCdf(etaLowerOffset, halfWidth);
        const double phiNormalization =
            triangularKernelCdf(phiUpperOffset, halfWidth) -
            triangularKernelCdf(phiLowerOffset, halfWidth);
        if (etaNormalization <= 0.0 || phiNormalization <= 0.0) {
            continue;
        }

        const double componentWeight = std::abs(particle.charge);
        if (componentWeight <= 0.0 || !std::isfinite(componentWeight)) {
            continue;
        }
        components.push_back(SmoothEtaPhiProposalComponent{
            etaCenter,
            phiCenter,
            halfWidth,
            halfWidth,
            etaLowerOffset,
            etaUpperOffset,
            phiLowerOffset,
            phiUpperOffset,
            etaNormalization,
            phiNormalization,
            componentWeight,
        });
        totalComponentWeight += componentWeight;
    }

    if (components.empty() || totalComponentWeight <= 0.0) {
        throw std::invalid_argument("smooth eta-phi proposal has no active components");
    }

    const double componentFraction = 1.0 - config.uniformFraction;
    for (SmoothEtaPhiProposalComponent& component : components) {
        component.mixtureWeight = componentFraction * component.mixtureWeight /
                                  totalComponentWeight;
    }

    return SmoothEtaPhiProposal{phaseSpace, config, components};
}

[[nodiscard]] inline double smoothEtaPhiDensity(
    const SmoothEtaPhiProposal& proposal,
    double eta,
    double phi) {
    requireValidEtaPhiProposalRange(proposal.phaseSpace);
    validateSmoothEtaPhiProposalConfig(proposal.config);

    if (!std::isfinite(eta) || !std::isfinite(phi) ||
        eta < proposal.phaseSpace.etaMin || eta > proposal.phaseSpace.etaMax ||
        phi < proposal.phaseSpace.phiMin || phi > proposal.phaseSpace.phiMax) {
        throw std::invalid_argument("eta-phi point is outside proposal phase space");
    }

    double density = proposal.config.uniformFraction * uniformEtaPhiDensity(proposal.phaseSpace);
    for (const SmoothEtaPhiProposalComponent& component : proposal.components) {
        const double etaOffset = eta - component.etaCenter;
        const double phiOffset = phi - component.phiCenter;
        if (etaOffset < component.etaLowerOffset || etaOffset > component.etaUpperOffset ||
            phiOffset < component.phiLowerOffset || phiOffset > component.phiUpperOffset) {
            continue;
        }
        const double etaDensity =
            triangularKernelPdf(etaOffset, component.etaHalfWidth) /
            component.etaNormalization;
        const double phiDensity =
            triangularKernelPdf(phiOffset, component.phiHalfWidth) /
            component.phiNormalization;
        density += component.mixtureWeight * etaDensity * phiDensity;
    }

    return density;
}

[[nodiscard]] inline SmoothEtaPhiProposalPoint sampleSmoothEtaPhi(
    const SmoothEtaPhiProposal& proposal,
    double unitComponent,
    double unitEta,
    double unitPhi) {
    requireValidEtaPhiProposalRange(proposal.phaseSpace);
    validateSmoothEtaPhiProposalConfig(proposal.config);

    if (!std::isfinite(unitComponent) || unitComponent < 0.0 || unitComponent >= 1.0) {
        throw std::invalid_argument("component unit random value must satisfy 0 <= u < 1");
    }
    if (!std::isfinite(unitEta) || unitEta < 0.0 || unitEta >= 1.0) {
        throw std::invalid_argument("eta unit random value must satisfy 0 <= u < 1");
    }
    if (!std::isfinite(unitPhi) || unitPhi < 0.0 || unitPhi >= 1.0) {
        throw std::invalid_argument("phi unit random value must satisfy 0 <= u < 1");
    }

    double eta = 0.0;
    double phi = 0.0;
    if (unitComponent < proposal.config.uniformFraction) {
        const UniformEtaPhiProposalPoint uniform =
            sampleUniformEtaPhi(proposal.phaseSpace, unitEta, unitPhi);
        eta = uniform.eta;
        phi = uniform.phi;
    } else {
        const double target = unitComponent - proposal.config.uniformFraction;
        double cumulative = 0.0;
        const SmoothEtaPhiProposalComponent* selected = &proposal.components.back();
        for (const SmoothEtaPhiProposalComponent& component : proposal.components) {
            cumulative += component.mixtureWeight;
            if (target < cumulative) {
                selected = &component;
                break;
            }
        }

        eta = selected->etaCenter + sampleTruncatedTriangularKernel(
                                       unitEta,
                                       selected->etaLowerOffset,
                                       selected->etaUpperOffset,
                                       selected->etaHalfWidth);
        phi = selected->phiCenter + sampleTruncatedTriangularKernel(
                                       unitPhi,
                                       selected->phiLowerOffset,
                                       selected->phiUpperOffset,
                                       selected->phiHalfWidth);
        if (eta >= proposal.phaseSpace.etaMax) {
            eta = std::nextafter(proposal.phaseSpace.etaMax, proposal.phaseSpace.etaMin);
        }
        if (phi >= proposal.phaseSpace.phiMax) {
            phi = std::nextafter(proposal.phaseSpace.phiMax, proposal.phaseSpace.phiMin);
        }
    }

    return SmoothEtaPhiProposalPoint{
        eta,
        phi,
        smoothEtaPhiDensity(proposal, eta, phi),
    };
}

}  // namespace lowmc
