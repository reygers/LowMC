#include "lowmc/GridEtaPhiProposal.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace lowmc {
namespace {

[[nodiscard]] FourVector photonFromKTEtaPhi(double kT, double eta, double phi) {
    const double px = kT * std::cos(phi);
    const double py = kT * std::sin(phi);
    const double pz = kT * std::sinh(eta);
    const double energy = kT * std::cosh(eta);
    return FourVector{energy, px, py, pz};
}

[[nodiscard]] double cellArea(const GridAdaptiveEtaPhiCell& cell) {
    return (cell.etaMax - cell.etaMin) * (cell.phiMax - cell.phiMin);
}

[[nodiscard]] double totalArea(const PhotonPhaseSpace& phaseSpace) {
    return (phaseSpace.etaMax - phaseSpace.etaMin) *
           (phaseSpace.phiMax - phaseSpace.phiMin);
}

[[nodiscard]] double cellTargetMass(
    const LowCurrent& current,
    double kT,
    const GridAdaptiveEtaPhiCell& cell) {
    const double eta = 0.5 * (cell.etaMin + cell.etaMax);
    const double phi = 0.5 * (cell.phiMin + cell.phiMax);
    const double intensity = current.intensity(photonFromKTEtaPhi(kT, eta, phi));
    if (!std::isfinite(intensity) || intensity < 0.0) {
        throw std::runtime_error("grid-adaptive eta-phi proposal found invalid intensity");
    }
    return intensity * cellArea(cell);
}

[[nodiscard]] double cellVarianceScore(
    const LowCurrent& current,
    double kT,
    const GridAdaptiveEtaPhiCell& cell,
    int probeBins) {
    double sum = 0.0;
    double sumSquares = 0.0;
    int count = 0;
    for (int etaIndex = 0; etaIndex < probeBins; ++etaIndex) {
        const double eta =
            cell.etaMin + (cell.etaMax - cell.etaMin) *
                              (static_cast<double>(etaIndex) + 0.5) /
                              static_cast<double>(probeBins);
        for (int phiIndex = 0; phiIndex < probeBins; ++phiIndex) {
            const double phi =
                cell.phiMin + (cell.phiMax - cell.phiMin) *
                                  (static_cast<double>(phiIndex) + 0.5) /
                                  static_cast<double>(probeBins);
            const double intensity = current.intensity(photonFromKTEtaPhi(kT, eta, phi));
            if (!std::isfinite(intensity) || intensity < 0.0) {
                throw std::runtime_error(
                    "grid-adaptive eta-phi proposal found invalid intensity");
            }
            sum += intensity;
            sumSquares += intensity * intensity;
            ++count;
        }
    }

    if (count <= 1) {
        return 0.0;
    }
    const double mean = sum / static_cast<double>(count);
    const double variance = std::max(
        0.0,
        sumSquares / static_cast<double>(count) - mean * mean);
    const double area = cellArea(cell);
    return area * area * variance;
}

[[nodiscard]] std::vector<GridAdaptiveEtaPhiCell> splitCell(
    const GridAdaptiveEtaPhiCell& cell) {
    const double etaMid = 0.5 * (cell.etaMin + cell.etaMax);
    const double phiMid = 0.5 * (cell.phiMin + cell.phiMax);
    return {
        {cell.etaMin, etaMid, cell.phiMin, phiMid},
        {etaMid, cell.etaMax, cell.phiMin, phiMid},
        {cell.etaMin, etaMid, phiMid, cell.phiMax},
        {etaMid, cell.etaMax, phiMid, cell.phiMax},
    };
}

void updateTargetMasses(
    const LowCurrent& current,
    double kT,
    std::vector<GridAdaptiveEtaPhiCell>& cells) {
    for (GridAdaptiveEtaPhiCell& cell : cells) {
        cell.targetMass = cellTargetMass(current, kT, cell);
    }
}

void updateRefinementScores(
    const LowCurrent& current,
    double kT,
    const GridAdaptiveEtaPhiProposalConfig& config,
    std::vector<GridAdaptiveEtaPhiCell>& cells) {
    for (GridAdaptiveEtaPhiCell& cell : cells) {
        if (config.refinementScore == GridAdaptiveRefinementScore::Variance) {
            cell.refinementScore =
                cellVarianceScore(current, kT, cell, config.varianceProbeBins);
        } else {
            cell.refinementScore = cell.targetMass;
        }
    }
}

}  // namespace

void validateGridAdaptiveEtaPhiProposalConfig(
    const GridAdaptiveEtaPhiProposalConfig& config) {
    if (config.etaBins < 1 || config.phiBins < 1) {
        throw std::invalid_argument("grid-adaptive eta/phi bin counts must be positive");
    }
    if (config.refinementLevels < 0) {
        throw std::invalid_argument("grid-adaptive refinement levels must be non-negative");
    }
    if (!std::isfinite(config.refineFraction) ||
        config.refineFraction <= 0.0 ||
        config.refineFraction > 1.0) {
        throw std::invalid_argument("grid-adaptive refineFraction must satisfy 0 < f <= 1");
    }
    if (!std::isfinite(config.uniformFraction) ||
        config.uniformFraction < 0.0 ||
        config.uniformFraction >= 1.0) {
        throw std::invalid_argument("grid-adaptive uniformFraction must satisfy 0 <= f < 1");
    }
    if (config.varianceProbeBins < 2) {
        throw std::invalid_argument(
            "grid-adaptive varianceProbeBins must be at least two");
    }
}

GridAdaptiveEtaPhiProposal makeGridAdaptiveEtaPhiProposal(
    const LowCurrent& current,
    const PhotonPhaseSpace& phaseSpace,
    double kT,
    const GridAdaptiveEtaPhiProposalConfig& config) {
    if (!validatePhotonPhaseSpace(phaseSpace).valid()) {
        throw std::invalid_argument("invalid photon phase space");
    }
    validateGridAdaptiveEtaPhiProposalConfig(config);
    if (!std::isfinite(kT) || kT < phaseSpace.kTMin || kT > phaseSpace.kTMax) {
        throw std::invalid_argument("fixed kT must be finite and inside the phase-space range");
    }

    std::vector<GridAdaptiveEtaPhiCell> cells;
    cells.reserve(static_cast<std::size_t>(config.etaBins * config.phiBins));
    for (int etaBin = 0; etaBin < config.etaBins; ++etaBin) {
        const double etaMin =
            phaseSpace.etaMin + (phaseSpace.etaMax - phaseSpace.etaMin) *
                                    static_cast<double>(etaBin) /
                                    static_cast<double>(config.etaBins);
        const double etaMax =
            phaseSpace.etaMin + (phaseSpace.etaMax - phaseSpace.etaMin) *
                                    static_cast<double>(etaBin + 1) /
                                    static_cast<double>(config.etaBins);
        for (int phiBin = 0; phiBin < config.phiBins; ++phiBin) {
            const double phiMin =
                phaseSpace.phiMin + (phaseSpace.phiMax - phaseSpace.phiMin) *
                                        static_cast<double>(phiBin) /
                                        static_cast<double>(config.phiBins);
            const double phiMax =
                phaseSpace.phiMin + (phaseSpace.phiMax - phaseSpace.phiMin) *
                                        static_cast<double>(phiBin + 1) /
                                        static_cast<double>(config.phiBins);
            cells.push_back(GridAdaptiveEtaPhiCell{etaMin, etaMax, phiMin, phiMax});
        }
    }

    updateTargetMasses(current, kT, cells);
    updateRefinementScores(current, kT, config, cells);
    for (int level = 0; level < config.refinementLevels; ++level) {
        const std::size_t refineCount = std::max<std::size_t>(
            1,
            static_cast<std::size_t>(
                std::ceil(config.refineFraction * static_cast<double>(cells.size()))));
        std::vector<std::size_t> indices(cells.size());
        for (std::size_t index = 0; index < indices.size(); ++index) {
            indices[index] = index;
        }
        std::partial_sort(
            indices.begin(),
            indices.begin() + static_cast<std::ptrdiff_t>(refineCount),
            indices.end(),
            [&cells](std::size_t left, std::size_t right) {
                return cells[left].refinementScore > cells[right].refinementScore;
            });

        std::vector<bool> refine(cells.size(), false);
        for (std::size_t index = 0; index < refineCount; ++index) {
            refine[indices[index]] = true;
        }

        std::vector<GridAdaptiveEtaPhiCell> nextCells;
        nextCells.reserve(cells.size() + 3 * refineCount);
        for (std::size_t index = 0; index < cells.size(); ++index) {
            if (!refine[index]) {
                nextCells.push_back(cells[index]);
                continue;
            }
            std::vector<GridAdaptiveEtaPhiCell> children = splitCell(cells[index]);
            nextCells.insert(nextCells.end(), children.begin(), children.end());
        }
        cells = std::move(nextCells);
        updateTargetMasses(current, kT, cells);
        updateRefinementScores(current, kT, config, cells);
    }

    double totalTargetMass = 0.0;
    for (const GridAdaptiveEtaPhiCell& cell : cells) {
        totalTargetMass += cell.targetMass;
    }
    if (!(totalTargetMass > 0.0) || !std::isfinite(totalTargetMass)) {
        throw std::runtime_error("grid-adaptive eta-phi proposal has zero target mass");
    }

    const double adaptiveFraction = 1.0 - config.uniformFraction;
    for (GridAdaptiveEtaPhiCell& cell : cells) {
        cell.mixtureWeight = adaptiveFraction * cell.targetMass / totalTargetMass;
    }

    return GridAdaptiveEtaPhiProposal{phaseSpace, config, cells};
}

double gridAdaptiveEtaPhiDensity(
    const GridAdaptiveEtaPhiProposal& proposal,
    double eta,
    double phi) {
    if (!std::isfinite(eta) || !std::isfinite(phi) ||
        eta < proposal.phaseSpace.etaMin || eta > proposal.phaseSpace.etaMax ||
        phi < proposal.phaseSpace.phiMin || phi > proposal.phaseSpace.phiMax) {
        throw std::invalid_argument("eta-phi point is outside proposal phase space");
    }

    double density = proposal.config.uniformFraction / totalArea(proposal.phaseSpace);
    for (const GridAdaptiveEtaPhiCell& cell : proposal.cells) {
        if (eta >= cell.etaMin && eta <= cell.etaMax &&
            phi >= cell.phiMin && phi <= cell.phiMax) {
            density += cell.mixtureWeight / cellArea(cell);
            break;
        }
    }
    return density;
}

GridAdaptiveEtaPhiProposalPoint sampleGridAdaptiveEtaPhi(
    const GridAdaptiveEtaPhiProposal& proposal,
    double unitComponent,
    double unitCell,
    double unitEta,
    double unitPhi) {
    if (!std::isfinite(unitComponent) || unitComponent < 0.0 || unitComponent >= 1.0) {
        throw std::invalid_argument("component unit random value must satisfy 0 <= u < 1");
    }
    if (!std::isfinite(unitCell) || unitCell < 0.0 || unitCell >= 1.0) {
        throw std::invalid_argument("cell unit random value must satisfy 0 <= u < 1");
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
        eta = proposal.phaseSpace.etaMin +
              unitEta * (proposal.phaseSpace.etaMax - proposal.phaseSpace.etaMin);
        phi = proposal.phaseSpace.phiMin +
              unitPhi * (proposal.phaseSpace.phiMax - proposal.phaseSpace.phiMin);
        if (eta >= proposal.phaseSpace.etaMax) {
            eta = std::nextafter(proposal.phaseSpace.etaMax, proposal.phaseSpace.etaMin);
        }
        if (phi >= proposal.phaseSpace.phiMax) {
            phi = std::nextafter(proposal.phaseSpace.phiMax, proposal.phaseSpace.phiMin);
        }
    } else {
        const double target = unitCell * (1.0 - proposal.config.uniformFraction);
        double cumulative = 0.0;
        const GridAdaptiveEtaPhiCell* selected = &proposal.cells.back();
        for (const GridAdaptiveEtaPhiCell& cell : proposal.cells) {
            cumulative += cell.mixtureWeight;
            if (target < cumulative) {
                selected = &cell;
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

    return GridAdaptiveEtaPhiProposalPoint{
        eta,
        phi,
        gridAdaptiveEtaPhiDensity(proposal, eta, phi),
    };
}

}  // namespace lowmc
