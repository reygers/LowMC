#include "lowmc/EventYield.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>

#include "lowmc/FixedKTWeightedSampler.h"

namespace lowmc {
namespace {

void validateIntegrationConfig(const EventYieldIntegrationConfig& config) {
    if (config.logKTBins <= 0 || config.etaBins <= 0 || config.phiBins <= 0) {
        throw std::invalid_argument("event-yield integration bins must be positive");
    }
    if (config.logKTBins > std::numeric_limits<int>::max() / config.etaBins ||
        config.logKTBins * config.etaBins >
            std::numeric_limits<int>::max() / config.phiBins) {
        throw std::invalid_argument("event-yield integration grid is too large");
    }
}

void validateAdaptiveIntegrationConfig(const AdaptiveEventYieldIntegrationConfig& config) {
    if (config.logKTBins <= 0 || config.etaSubBins <= 0 || config.phiSubBins <= 0) {
        throw std::invalid_argument("adaptive event-yield integration bins must be positive");
    }
    validateGridAdaptiveEtaPhiProposalConfig(config.gridConfig);
}

double midpointIntegral(
    const LowCurrent& current,
    const PhotonPhaseSpace& phaseSpace,
    const EventYieldIntegrationConfig& config,
    double alpha) {
    const double logKTMin = std::log(phaseSpace.kTMin);
    const double logKTMax = std::log(phaseSpace.kTMax);
    const double dLogKT = (logKTMax - logKTMin) / static_cast<double>(config.logKTBins);
    const double dEta =
        (phaseSpace.etaMax - phaseSpace.etaMin) / static_cast<double>(config.etaBins);
    const double dPhi =
        (phaseSpace.phiMax - phaseSpace.phiMin) / static_cast<double>(config.phiBins);

    double sum = 0.0;
    for (int iKT = 0; iKT < config.logKTBins; ++iKT) {
        const double logKT = logKTMin + (static_cast<double>(iKT) + 0.5) * dLogKT;
        const double kT = std::exp(logKT);
        const double measureJacobian = kT * kT;
        for (int iEta = 0; iEta < config.etaBins; ++iEta) {
            const double eta = phaseSpace.etaMin + (static_cast<double>(iEta) + 0.5) * dEta;
            for (int iPhi = 0; iPhi < config.phiBins; ++iPhi) {
                const double phi =
                    phaseSpace.phiMin + (static_cast<double>(iPhi) + 0.5) * dPhi;
                const FourVector photon = photonMomentumFromKTEtaPhi(kT, eta, phi);
                sum += lowPhotonDensity(current, photon, alpha) * measureJacobian;
            }
        }
    }

    return sum * dLogKT * dEta * dPhi;
}

double adaptiveIntegral(
    const LowCurrent& current,
    const PhotonPhaseSpace& phaseSpace,
    const AdaptiveEventYieldIntegrationConfig& config,
    double alpha,
    int& evaluations,
    int& gridCells) {
    const double logKTMin = std::log(phaseSpace.kTMin);
    const double logKTMax = std::log(phaseSpace.kTMax);
    const double dLogKT = (logKTMax - logKTMin) / static_cast<double>(config.logKTBins);

    double sum = 0.0;
    for (int iKT = 0; iKT < config.logKTBins; ++iKT) {
        const double logKT = logKTMin + (static_cast<double>(iKT) + 0.5) * dLogKT;
        const double kT = std::exp(logKT);
        const double measureJacobian = kT * kT;
        const GridAdaptiveEtaPhiProposal proposal =
            makeGridAdaptiveEtaPhiProposal(current, phaseSpace, kT, config.gridConfig);
        gridCells += static_cast<int>(proposal.cells.size());

        for (const GridAdaptiveEtaPhiCell& cell : proposal.cells) {
            const double dEta =
                (cell.etaMax - cell.etaMin) / static_cast<double>(config.etaSubBins);
            const double dPhi =
                (cell.phiMax - cell.phiMin) / static_cast<double>(config.phiSubBins);
            for (int etaIndex = 0; etaIndex < config.etaSubBins; ++etaIndex) {
                const double eta =
                    cell.etaMin + (static_cast<double>(etaIndex) + 0.5) * dEta;
                for (int phiIndex = 0; phiIndex < config.phiSubBins; ++phiIndex) {
                    const double phi =
                        cell.phiMin + (static_cast<double>(phiIndex) + 0.5) * dPhi;
                    const FourVector photon = photonMomentumFromKTEtaPhi(kT, eta, phi);
                    sum += lowPhotonDensity(current, photon, alpha) *
                           measureJacobian * dEta * dPhi;
                    ++evaluations;
                }
            }
        }
    }

    return sum * dLogKT;
}

double factorizedAdaptiveIntegral(
    const LowCurrent& current,
    const PhotonPhaseSpace& phaseSpace,
    const AdaptiveEventYieldIntegrationConfig& config,
    double alpha,
    int& evaluations,
    int& gridCells) {
    const double logKTInterval = std::log(phaseSpace.kTMax / phaseSpace.kTMin);
    const double kT = std::sqrt(phaseSpace.kTMin * phaseSpace.kTMax);
    const double measureJacobian = kT * kT;
    const GridAdaptiveEtaPhiProposal proposal =
        makeGridAdaptiveEtaPhiProposal(current, phaseSpace, kT, config.gridConfig);
    gridCells += static_cast<int>(proposal.cells.size());

    double angularSum = 0.0;
    for (const GridAdaptiveEtaPhiCell& cell : proposal.cells) {
        const double dEta =
            (cell.etaMax - cell.etaMin) / static_cast<double>(config.etaSubBins);
        const double dPhi =
            (cell.phiMax - cell.phiMin) / static_cast<double>(config.phiSubBins);
        for (int etaIndex = 0; etaIndex < config.etaSubBins; ++etaIndex) {
            const double eta = cell.etaMin + (static_cast<double>(etaIndex) + 0.5) * dEta;
            for (int phiIndex = 0; phiIndex < config.phiSubBins; ++phiIndex) {
                const double phi = cell.phiMin + (static_cast<double>(phiIndex) + 0.5) * dPhi;
                const FourVector photon = photonMomentumFromKTEtaPhi(kT, eta, phi);
                angularSum +=
                    lowPhotonDensity(current, photon, alpha) * measureJacobian * dEta * dPhi;
                ++evaluations;
            }
        }
    }

    return angularSum * logKTInterval;
}

EventYieldIntegrationConfig coarseConfig(EventYieldIntegrationConfig config) {
    config.logKTBins = std::max(1, config.logKTBins / 2);
    config.etaBins = std::max(1, config.etaBins / 2);
    config.phiBins = std::max(1, config.phiBins / 2);
    return config;
}

AdaptiveEventYieldIntegrationConfig coarseConfig(AdaptiveEventYieldIntegrationConfig config) {
    config.logKTBins = std::max(1, config.logKTBins / 2);
    config.gridConfig.refinementLevels =
        std::max(0, config.gridConfig.refinementLevels - 1);
    return config;
}

}  // namespace

EventYieldEstimate estimateEventYield(
    const LowCurrent& current,
    const PhotonPhaseSpace& phaseSpace,
    const EventYieldIntegrationConfig& config,
    double alpha) {
    if (!validatePhotonPhaseSpace(phaseSpace).valid()) {
        throw std::invalid_argument("invalid photon phase space");
    }
    validateIntegrationConfig(config);
    if (!std::isfinite(alpha) || alpha <= 0.0) {
        throw std::invalid_argument("fine-structure constant must be positive and finite");
    }

    const EventYieldIntegrationConfig coarse = coarseConfig(config);
    validateIntegrationConfig(coarse);

    EventYieldEstimate estimate;
    estimate.lambdaEvent = midpointIntegral(current, phaseSpace, config, alpha);
    estimate.coarseLambdaEvent = midpointIntegral(current, phaseSpace, coarse, alpha);
    estimate.absoluteDifference = std::abs(estimate.lambdaEvent - estimate.coarseLambdaEvent);
    estimate.relativeDifference =
        estimate.lambdaEvent != 0.0 ? estimate.absoluteDifference / std::abs(estimate.lambdaEvent)
                                    : 0.0;
    estimate.logKTBins = config.logKTBins;
    estimate.etaBins = config.etaBins;
    estimate.phiBins = config.phiBins;
    estimate.evaluations = config.logKTBins * config.etaBins * config.phiBins;
    estimate.status =
        std::isfinite(estimate.lambdaEvent) && estimate.lambdaEvent >= 0.0 ? "ok" : "invalid";
    return estimate;
}

EventYieldEstimate estimateEventYield(
    const std::vector<ChargedParticle>& selectedParticles,
    const PhotonPhaseSpace& phaseSpace,
    const EventYieldIntegrationConfig& config,
    double alpha) {
    return estimateEventYield(LowCurrent(selectedParticles), phaseSpace, config, alpha);
}

EventYieldEstimate estimateAdaptiveEventYield(
    const LowCurrent& current,
    const PhotonPhaseSpace& phaseSpace,
    const AdaptiveEventYieldIntegrationConfig& config,
    double alpha) {
    if (!validatePhotonPhaseSpace(phaseSpace).valid()) {
        throw std::invalid_argument("invalid photon phase space");
    }
    validateAdaptiveIntegrationConfig(config);
    if (!std::isfinite(alpha) || alpha <= 0.0) {
        throw std::invalid_argument("fine-structure constant must be positive and finite");
    }

    const AdaptiveEventYieldIntegrationConfig coarse = coarseConfig(config);
    validateAdaptiveIntegrationConfig(coarse);

    EventYieldEstimate estimate;
    estimate.method = "adaptive-log-kt-grid-eta-phi-midpoint";
    estimate.logKTBins = config.logKTBins;
    estimate.etaBins = config.etaSubBins;
    estimate.phiBins = config.phiSubBins;
    estimate.lambdaEvent =
        adaptiveIntegral(
            current,
            phaseSpace,
            config,
            alpha,
            estimate.evaluations,
            estimate.gridCells);

    int coarseEvaluations = 0;
    int coarseGridCells = 0;
    estimate.coarseLambdaEvent =
        adaptiveIntegral(
            current,
            phaseSpace,
            coarse,
            alpha,
            coarseEvaluations,
            coarseGridCells);
    estimate.absoluteDifference = std::abs(estimate.lambdaEvent - estimate.coarseLambdaEvent);
    estimate.relativeDifference =
        estimate.lambdaEvent != 0.0 ? estimate.absoluteDifference / std::abs(estimate.lambdaEvent)
                                    : 0.0;
    estimate.status =
        std::isfinite(estimate.lambdaEvent) && estimate.lambdaEvent >= 0.0 ? "ok" : "invalid";
    return estimate;
}

EventYieldEstimate estimateAdaptiveEventYield(
    const std::vector<ChargedParticle>& selectedParticles,
    const PhotonPhaseSpace& phaseSpace,
    const AdaptiveEventYieldIntegrationConfig& config,
    double alpha) {
    return estimateAdaptiveEventYield(LowCurrent(selectedParticles), phaseSpace, config, alpha);
}

EventYieldEstimate estimateFactorizedAdaptiveEventYield(
    const LowCurrent& current,
    const PhotonPhaseSpace& phaseSpace,
    const AdaptiveEventYieldIntegrationConfig& config,
    double alpha) {
    if (!validatePhotonPhaseSpace(phaseSpace).valid()) {
        throw std::invalid_argument("invalid photon phase space");
    }
    validateAdaptiveIntegrationConfig(config);
    if (!std::isfinite(alpha) || alpha <= 0.0) {
        throw std::invalid_argument("fine-structure constant must be positive and finite");
    }

    const AdaptiveEventYieldIntegrationConfig coarse = coarseConfig(config);
    validateAdaptiveIntegrationConfig(coarse);

    EventYieldEstimate estimate;
    estimate.method = "factorized-adaptive-eta-phi-midpoint";
    estimate.logKTBins = 1;
    estimate.etaBins = config.etaSubBins;
    estimate.phiBins = config.phiSubBins;
    estimate.lambdaEvent =
        factorizedAdaptiveIntegral(
            current,
            phaseSpace,
            config,
            alpha,
            estimate.evaluations,
            estimate.gridCells);

    int coarseEvaluations = 0;
    int coarseGridCells = 0;
    estimate.coarseLambdaEvent =
        factorizedAdaptiveIntegral(
            current,
            phaseSpace,
            coarse,
            alpha,
            coarseEvaluations,
            coarseGridCells);
    estimate.absoluteDifference = std::abs(estimate.lambdaEvent - estimate.coarseLambdaEvent);
    estimate.relativeDifference =
        estimate.lambdaEvent != 0.0 ? estimate.absoluteDifference / std::abs(estimate.lambdaEvent)
                                    : 0.0;
    estimate.status =
        std::isfinite(estimate.lambdaEvent) && estimate.lambdaEvent >= 0.0 ? "ok" : "invalid";
    return estimate;
}

EventYieldEstimate estimateFactorizedAdaptiveEventYield(
    const std::vector<ChargedParticle>& selectedParticles,
    const PhotonPhaseSpace& phaseSpace,
    const AdaptiveEventYieldIntegrationConfig& config,
    double alpha) {
    return estimateFactorizedAdaptiveEventYield(
        LowCurrent(selectedParticles),
        phaseSpace,
        config,
        alpha);
}

PhotonMultiplicityDiagnostic samplePhotonMultiplicity(
    double lambdaEvent,
    const PhotonMultiplicityConfig& config) {
    if (!std::isfinite(lambdaEvent) || lambdaEvent < 0.0) {
        throw std::invalid_argument("lambda_event must be finite and non-negative");
    }
    if (!std::isfinite(config.oversamplingFactor) || config.oversamplingFactor <= 0.0) {
        throw std::invalid_argument("oversampling factor must be positive and finite");
    }

    const double poissonMean = lambdaEvent * config.oversamplingFactor;
    if (!std::isfinite(poissonMean)) {
        throw std::invalid_argument("oversampled Poisson mean must be finite");
    }

    std::mt19937_64 rng(config.seed);
    std::poisson_distribution<int> distribution(poissonMean);
    const int sampledPhotons = distribution(rng);
    const double sampleWeight = 1.0 / config.oversamplingFactor;

    return PhotonMultiplicityDiagnostic{
        lambdaEvent,
        poissonMean,
        sampledPhotons,
        sampleWeight,
        sampleWeight * static_cast<double>(sampledPhotons),
        config.oversamplingFactor,
        config.seed,
    };
}

}  // namespace lowmc
