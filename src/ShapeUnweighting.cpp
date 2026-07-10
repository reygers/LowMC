#include "lowmc/ShapeUnweighting.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iterator>
#include <limits>
#include <random>
#include <stdexcept>

#include "lowmc/KTProposal.h"
#include "lowmc/WeightedPhotonGenerator.h"

namespace lowmc {
namespace {

struct TrialRecord {
    GridAdaptiveWeightedPhotonSample sample;
    std::size_t cellIndex = 0;
};

struct EnvelopeCell {
    double bound = 0.0;
    double proposalMass = 0.0;
    double envelopeMass = 0.0;
};

struct PointCellIndex {
    double etaMin = 0.0;
    double etaMax = 0.0;
    double etaBinWidth = 0.0;
    int etaBins = 0;
    std::vector<std::vector<std::size_t>> etaBuckets;
};

struct EnvelopeCellSampler {
    std::vector<double> cumulativeMasses;
};

struct ProbePlan {
    std::vector<int> bins;
    std::vector<double> safetyFactors;
    int evidenceCellCount = 0;
    int promotedCellCount = 0;
    std::size_t probeEvaluations = 0;
};

[[nodiscard]] double secondsSince(std::chrono::steady_clock::time_point start) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}

[[nodiscard]] double cellArea(const GridAdaptiveEtaPhiCell& cell) {
    return (cell.etaMax - cell.etaMin) * (cell.phiMax - cell.phiMin);
}

[[nodiscard]] double totalEtaPhiArea(const PhotonPhaseSpace& phaseSpace) {
    return (phaseSpace.etaMax - phaseSpace.etaMin) *
           (phaseSpace.phiMax - phaseSpace.phiMin);
}

[[nodiscard]] int effectiveBaseProbeBins(const ShapeUnweightingConfig& config) {
    return config.baseProbeBins > 0 ? config.baseProbeBins : config.probeBins;
}

[[nodiscard]] double cellProposalMass(
    const GridAdaptiveEtaPhiProposal& proposal,
    const GridAdaptiveEtaPhiCell& cell) {
    return proposal.config.uniformFraction *
               cellArea(cell) / totalEtaPhiArea(proposal.phaseSpace) +
           cell.mixtureWeight;
}

[[nodiscard]] double etaPhiDensityForCell(
    const GridAdaptiveEtaPhiProposal& proposal,
    const GridAdaptiveEtaPhiCell& cell) {
    return proposal.config.uniformFraction / totalEtaPhiArea(proposal.phaseSpace) +
           cell.mixtureWeight / cellArea(cell);
}

[[nodiscard]] double probeCoordinate(
    double lower,
    double upper,
    int index,
    int bins,
    ShapeUnweightingProbeRule rule) {
    if (rule == ShapeUnweightingProbeRule::Interior) {
        return lower + (upper - lower) *
                           (static_cast<double>(index) + 0.5) /
                           static_cast<double>(bins);
    }

    if (bins == 1) {
        return 0.5 * (lower + upper);
    }
    const double fraction = static_cast<double>(index) / static_cast<double>(bins - 1);
    if (index == 0) {
        return std::nextafter(lower, upper);
    }
    if (index == bins - 1) {
        return std::nextafter(upper, lower);
    }
    return lower + (upper - lower) * fraction;
}

[[nodiscard]] double wrapPhiIntoRange(double phi, const PhotonPhaseSpace& phaseSpace) {
    const double period = phaseSpace.phiMax - phaseSpace.phiMin;
    if (!(period > 0.0) || !std::isfinite(period)) {
        return phi;
    }
    while (phi < phaseSpace.phiMin) {
        phi += period;
    }
    while (phi >= phaseSpace.phiMax) {
        phi -= period;
    }
    return phi;
}

[[nodiscard]] bool chargedParticleEtaPhi(
    const ChargedParticle& particle,
    const PhotonPhaseSpace& phaseSpace,
    double& eta,
    double& phi) {
    const double pT =
        std::hypot(particle.momentum.px, particle.momentum.py);
    if (!(pT > 0.0) || !std::isfinite(pT)) {
        return false;
    }
    eta = std::asinh(particle.momentum.pz / pT);
    phi = wrapPhiIntoRange(std::atan2(particle.momentum.py, particle.momentum.px), phaseSpace);
    return std::isfinite(eta) && std::isfinite(phi);
}

[[nodiscard]] double intervalDistance(double value, double lower, double upper) {
    if (value >= lower && value <= upper) {
        return 0.0;
    }
    return std::min(std::abs(value - lower), std::abs(value - upper));
}

[[nodiscard]] bool cellNearCoordinate(
    const GridAdaptiveEtaPhiCell& cell,
    double eta,
    double phi,
    const PhotonPhaseSpace& phaseSpace,
    double neighborScale) {
    const double etaWidth = cell.etaMax - cell.etaMin;
    const double phiWidth = cell.phiMax - cell.phiMin;
    const double etaDistance = intervalDistance(eta, cell.etaMin, cell.etaMax);
    const double phiDistance = intervalDistance(phi, cell.phiMin, cell.phiMax);
    const double phiPeriod = phaseSpace.phiMax - phaseSpace.phiMin;
    const double wrappedPhiDistance = std::min(
        phiDistance,
        std::min(
            intervalDistance(phi - phiPeriod, cell.phiMin, cell.phiMax),
            intervalDistance(phi + phiPeriod, cell.phiMin, cell.phiMax)));
    return etaDistance <= neighborScale * etaWidth &&
           wrappedPhiDistance <= neighborScale * phiWidth;
}

[[nodiscard]] std::vector<int> makeBaseProbeBins(
    const std::vector<ChargedParticle>& selectedParticles,
    const PhotonPhaseSpace& phaseSpace,
    const GridAdaptiveEtaPhiProposal& gridProposal,
    const ShapeUnweightingConfig& config) {
    std::vector<int> bins(gridProposal.cells.size(), effectiveBaseProbeBins(config));
    if (!config.targetedBaseProbes || effectiveBaseProbeBins(config) >= config.probeBins) {
        return bins;
    }

    for (const ChargedParticle& particle : selectedParticles) {
        double eta = 0.0;
        double phi = 0.0;
        if (!chargedParticleEtaPhi(particle, phaseSpace, eta, phi)) {
            continue;
        }
        if (eta < phaseSpace.etaMin || eta > phaseSpace.etaMax ||
            phi < phaseSpace.phiMin || phi >= phaseSpace.phiMax) {
            continue;
        }
        for (std::size_t index = 0; index < gridProposal.cells.size(); ++index) {
            if (cellNearCoordinate(
                    gridProposal.cells[index],
                    eta,
                    phi,
                    phaseSpace,
                    config.targetedBaseProbeNeighborScale)) {
                bins[index] = config.probeBins;
            }
        }
    }
    return bins;
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

void updateShapeUnweightingGridMasses(
    const LowCurrent& current,
    double kT,
    GridAdaptiveEtaPhiProposal& proposal) {
    double totalTargetMass = 0.0;
    for (GridAdaptiveEtaPhiCell& cell : proposal.cells) {
        const double eta = 0.5 * (cell.etaMin + cell.etaMax);
        const double phi = 0.5 * (cell.phiMin + cell.phiMax);
        const double intensity = current.intensity(photonMomentumFromKTEtaPhi(kT, eta, phi));
        if (!std::isfinite(intensity) || intensity < 0.0) {
            throw std::runtime_error("shape-unweighting targeted refinement found invalid intensity");
        }
        cell.targetMass = intensity * cellArea(cell);
        cell.refinementScore = cell.targetMass;
        totalTargetMass += cell.targetMass;
    }
    if (!(totalTargetMass > 0.0) || !std::isfinite(totalTargetMass)) {
        throw std::runtime_error("shape-unweighting targeted refinement has zero target mass");
    }

    const double adaptiveFraction = 1.0 - proposal.config.uniformFraction;
    for (GridAdaptiveEtaPhiCell& cell : proposal.cells) {
        cell.mixtureWeight = adaptiveFraction * cell.targetMass / totalTargetMass;
    }
}

void refineGridNearChargedParticles(
    const LowCurrent& current,
    const std::vector<ChargedParticle>& selectedParticles,
    double kT,
    int refinementLevels,
    double neighborScale,
    GridAdaptiveEtaPhiProposal& proposal) {
    if (refinementLevels <= 0) {
        return;
    }

    for (int level = 0; level < refinementLevels; ++level) {
        std::vector<bool> refine(proposal.cells.size(), false);
        for (const ChargedParticle& particle : selectedParticles) {
            double eta = 0.0;
            double phi = 0.0;
            if (!chargedParticleEtaPhi(particle, proposal.phaseSpace, eta, phi)) {
                continue;
            }
            if (eta < proposal.phaseSpace.etaMin || eta > proposal.phaseSpace.etaMax ||
                phi < proposal.phaseSpace.phiMin || phi >= proposal.phaseSpace.phiMax) {
                continue;
            }
            for (std::size_t index = 0; index < proposal.cells.size(); ++index) {
                if (cellNearCoordinate(
                        proposal.cells[index],
                        eta,
                        phi,
                        proposal.phaseSpace,
                        neighborScale)) {
                    refine[index] = true;
                }
            }
        }

        std::vector<GridAdaptiveEtaPhiCell> nextCells;
        nextCells.reserve(proposal.cells.size());
        for (std::size_t index = 0; index < proposal.cells.size(); ++index) {
            if (!refine[index]) {
                nextCells.push_back(proposal.cells[index]);
                continue;
            }
            std::vector<GridAdaptiveEtaPhiCell> children = splitCell(proposal.cells[index]);
            nextCells.insert(nextCells.end(), children.begin(), children.end());
        }
        proposal.cells = std::move(nextCells);
        updateShapeUnweightingGridMasses(current, kT, proposal);
    }
}

[[nodiscard]] std::size_t findCellIndex(
    const GridAdaptiveEtaPhiProposal& proposal,
    double eta,
    double phi) {
    for (std::size_t index = 0; index < proposal.cells.size(); ++index) {
        const GridAdaptiveEtaPhiCell& cell = proposal.cells[index];
        if (eta >= cell.etaMin && eta <= cell.etaMax &&
            phi >= cell.phiMin && phi <= cell.phiMax) {
            return index;
        }
    }
    throw std::runtime_error("sampled eta-phi point is outside all grid cells");
}

[[nodiscard]] PointCellIndex makePointCellIndex(const GridAdaptiveEtaPhiProposal& proposal) {
    PointCellIndex index;
    index.etaMin = proposal.phaseSpace.etaMin;
    index.etaMax = proposal.phaseSpace.etaMax;
    index.etaBins = std::max(
        1,
        static_cast<int>(std::sqrt(static_cast<double>(proposal.cells.size()))));
    index.etaBinWidth = (index.etaMax - index.etaMin) / static_cast<double>(index.etaBins);
    index.etaBuckets.resize(static_cast<std::size_t>(index.etaBins));

    for (std::size_t cellIndex = 0; cellIndex < proposal.cells.size(); ++cellIndex) {
        const GridAdaptiveEtaPhiCell& cell = proposal.cells[cellIndex];
        const int firstBin = std::max(
            0,
            static_cast<int>(std::floor((cell.etaMin - index.etaMin) / index.etaBinWidth)));
        const int lastBin = std::min(
            index.etaBins - 1,
            static_cast<int>(std::floor(
                (std::nextafter(cell.etaMax, cell.etaMin) - index.etaMin) /
                index.etaBinWidth)));
        for (int bin = firstBin; bin <= lastBin; ++bin) {
            index.etaBuckets[static_cast<std::size_t>(bin)].push_back(cellIndex);
        }
    }
    return index;
}

[[nodiscard]] std::size_t findCellIndex(
    const GridAdaptiveEtaPhiProposal& proposal,
    const PointCellIndex& index,
    double eta,
    double phi) {
    if (index.etaBins <= 0 || !(index.etaBinWidth > 0.0) ||
        eta < index.etaMin || eta > index.etaMax) {
        return findCellIndex(proposal, eta, phi);
    }

    int etaBin = static_cast<int>(std::floor((eta - index.etaMin) / index.etaBinWidth));
    etaBin = std::clamp(etaBin, 0, index.etaBins - 1);
    const std::vector<std::size_t>& candidates = index.etaBuckets[static_cast<std::size_t>(etaBin)];
    std::size_t selected = proposal.cells.size();
    for (const std::size_t cellIndex : candidates) {
        const GridAdaptiveEtaPhiCell& cell = proposal.cells[cellIndex];
        if (eta >= cell.etaMin && eta <= cell.etaMax &&
            phi >= cell.phiMin && phi <= cell.phiMax) {
            selected = std::min(selected, cellIndex);
        }
    }
    if (selected != proposal.cells.size()) {
        return selected;
    }
    return findCellIndex(proposal, eta, phi);
}

[[nodiscard]] TrialRecord makeTrialRecord(
    const LowCurrent& current,
    const PhotonPhaseSpace& phaseSpace,
    const GridAdaptiveEtaPhiProposal& gridProposal,
    const PointCellIndex& pointCellIndex,
    std::mt19937_64& rng,
    double alpha) {
    const double unitKT = std::generate_canonical<double, 53>(rng);
    const double unitEta = std::generate_canonical<double, 53>(rng);
    const double unitPhi = std::generate_canonical<double, 53>(rng);
    const double unitComponent = std::generate_canonical<double, 53>(rng);
    const double unitCell = std::generate_canonical<double, 53>(rng);

    GridAdaptiveWeightedPhotonSample sample =
        makeGridAdaptiveWeightedPhotonSample(
            current,
            phaseSpace,
            gridProposal,
            unitKT,
            unitComponent,
            unitCell,
            unitEta,
            unitPhi,
            alpha);
    const std::size_t cellIndex =
        findCellIndex(
            gridProposal,
            pointCellIndex,
            sample.etaPhiProposal.eta,
            sample.etaPhiProposal.phi);
    return TrialRecord{sample, cellIndex};
}

[[nodiscard]] EnvelopeCellSampler makeEnvelopeCellSampler(
    const std::vector<EnvelopeCell>& envelopeCells) {
    EnvelopeCellSampler sampler;
    sampler.cumulativeMasses.reserve(envelopeCells.size());
    double cumulative = 0.0;
    for (const EnvelopeCell& cell : envelopeCells) {
        cumulative += cell.envelopeMass;
        sampler.cumulativeMasses.push_back(cumulative);
    }
    return sampler;
}

[[nodiscard]] std::size_t sampleEnvelopeCell(
    const std::vector<EnvelopeCell>& envelopeCells,
    double totalEnvelopeMass,
    double unitCell) {
    if (!std::isfinite(unitCell) || unitCell < 0.0 || unitCell >= 1.0) {
        throw std::invalid_argument("envelope cell unit random value must satisfy 0 <= u < 1");
    }
    const double target = unitCell * totalEnvelopeMass;
    double cumulative = 0.0;
    std::size_t selected = envelopeCells.size() - 1;
    for (std::size_t index = 0; index < envelopeCells.size(); ++index) {
        cumulative += envelopeCells[index].envelopeMass;
        if (target < cumulative) {
            selected = index;
            break;
        }
    }
    return selected;
}

[[nodiscard]] std::size_t sampleEnvelopeCell(
    const EnvelopeCellSampler& sampler,
    double totalEnvelopeMass,
    double unitCell) {
    if (!std::isfinite(unitCell) || unitCell < 0.0 || unitCell >= 1.0) {
        throw std::invalid_argument("envelope cell unit random value must satisfy 0 <= u < 1");
    }
    const double target = unitCell * totalEnvelopeMass;
    const auto selected = std::upper_bound(
        sampler.cumulativeMasses.begin(),
        sampler.cumulativeMasses.end(),
        target);
    if (selected == sampler.cumulativeMasses.end()) {
        return sampler.cumulativeMasses.size() - 1;
    }
    return static_cast<std::size_t>(
        std::distance(sampler.cumulativeMasses.begin(), selected));
}

[[nodiscard]] TrialRecord makeEnvelopeTrialRecord(
    const LowCurrent& current,
    const PhotonPhaseSpace& phaseSpace,
    const GridAdaptiveEtaPhiProposal& gridProposal,
    const EnvelopeCellSampler& envelopeCellSampler,
    const PointCellIndex& pointCellIndex,
    double totalEnvelopeMass,
    std::mt19937_64& rng,
    double alpha) {
    const double unitKT = std::generate_canonical<double, 53>(rng);
    const double unitCell = std::generate_canonical<double, 53>(rng);
    const double unitEta = std::generate_canonical<double, 53>(rng);
    const double unitPhi = std::generate_canonical<double, 53>(rng);

    const std::size_t cellIndex =
        sampleEnvelopeCell(envelopeCellSampler, totalEnvelopeMass, unitCell);
    const GridAdaptiveEtaPhiCell& cell = gridProposal.cells[cellIndex];

    const LogUniformKTProposalPoint kTProposal = sampleLogUniformKT(phaseSpace, unitKT);
    double eta = cell.etaMin + unitEta * (cell.etaMax - cell.etaMin);
    double phi = cell.phiMin + unitPhi * (cell.phiMax - cell.phiMin);
    if (eta >= cell.etaMax) {
        eta = std::nextafter(cell.etaMax, cell.etaMin);
    }
    if (phi >= cell.phiMax) {
        phi = std::nextafter(cell.phiMax, cell.phiMin);
    }

    const std::size_t densityCellIndex = findCellIndex(gridProposal, pointCellIndex, eta, phi);
    const double etaPhiDensity =
        etaPhiDensityForCell(gridProposal, gridProposal.cells[densityCellIndex]);
    const double proposalDensity =
        kTProposal.densityPhysicalMeasure * etaPhiDensity;
    const FourVector momentum = photonMomentumFromKTEtaPhi(kTProposal.kT, eta, phi);
    const WeightedPhotonSample weightedPhoton = makeWeightedPhotonSample(
        current,
        phaseSpace,
        WeightedPhotonProposal{momentum, proposalDensity},
        alpha);

    return TrialRecord{
        GridAdaptiveWeightedPhotonSample{
            kTProposal,
            GridAdaptiveEtaPhiProposalPoint{eta, phi, etaPhiDensity},
            proposalDensity,
            weightedPhoton,
        },
        cellIndex,
    };
}

[[nodiscard]] double exactWeightAtEtaPhi(
    const LowCurrent& current,
    const PhotonPhaseSpace& phaseSpace,
    const GridAdaptiveEtaPhiProposal& gridProposal,
    const PointCellIndex& pointCellIndex,
    double kT,
    double eta,
    double phi,
    double alpha) {
    const std::size_t densityCellIndex = findCellIndex(gridProposal, pointCellIndex, eta, phi);
    const double proposalDensity =
        logUniformKTDensityPhysicalMeasure(phaseSpace, kT) *
        etaPhiDensityForCell(gridProposal, gridProposal.cells[densityCellIndex]);
    const WeightedPhotonSample sample =
        makeWeightedPhotonSample(
            current,
            phaseSpace,
            WeightedPhotonProposal{
                photonMomentumFromKTEtaPhi(kT, eta, phi),
                proposalDensity,
            },
            alpha);
    return sample.photon.weight;
}

[[nodiscard]] double exactWeightAtEtaPhiInKnownCell(
    const LowCurrent& current,
    const PhotonPhaseSpace& phaseSpace,
    const GridAdaptiveEtaPhiProposal& gridProposal,
    std::size_t cellIndex,
    double kT,
    double eta,
    double phi,
    double alpha) {
    const double proposalDensity =
        logUniformKTDensityPhysicalMeasure(phaseSpace, kT) *
        etaPhiDensityForCell(gridProposal, gridProposal.cells[cellIndex]);
    const WeightedPhotonSample sample =
        makeWeightedPhotonSample(
            current,
            phaseSpace,
            WeightedPhotonProposal{
                photonMomentumFromKTEtaPhi(kT, eta, phi),
                proposalDensity,
            },
            alpha);
    return sample.photon.weight;
}

[[nodiscard]] std::vector<EnvelopeCell> makeEnvelopeCells(
    const GridAdaptiveEtaPhiProposal& gridProposal,
    const std::vector<double>& probedCellMaxWeights,
    const std::vector<double>& calibrationCellMaxWeights,
    const std::vector<double>& cellSafetyFactors) {
    std::vector<EnvelopeCell> envelopeCells;
    envelopeCells.reserve(gridProposal.cells.size());
    for (std::size_t index = 0; index < gridProposal.cells.size(); ++index) {
        const double cellMax = std::max(
            probedCellMaxWeights[index],
            calibrationCellMaxWeights[index]);
        const double bound = cellSafetyFactors[index] * cellMax;
        const double proposalMass =
            cellProposalMass(gridProposal, gridProposal.cells[index]);
        const double envelopeMass =
            (std::isfinite(bound) && bound > 0.0)
                ? bound * proposalMass
                : 0.0;
        envelopeCells.push_back(EnvelopeCell{bound, proposalMass, envelopeMass});
    }
    return envelopeCells;
}

[[nodiscard]] double totalEnvelopeMass(const std::vector<EnvelopeCell>& envelopeCells) {
    double total = 0.0;
    for (const EnvelopeCell& cell : envelopeCells) {
        total += cell.envelopeMass;
    }
    return total;
}

[[nodiscard]] double probeSingleCellMaxWeight(
    const LowCurrent& current,
    const PhotonPhaseSpace& phaseSpace,
    const GridAdaptiveEtaPhiProposal& gridProposal,
    const PointCellIndex& pointCellIndex,
    std::size_t cellIndex,
    int probeBins,
    ShapeUnweightingProbeRule probeRule,
    double alpha) {
    const double kT = std::sqrt(phaseSpace.kTMin * phaseSpace.kTMax);
    const GridAdaptiveEtaPhiCell& cell = gridProposal.cells[cellIndex];
    double maxWeight = 0.0;
    for (int etaIndex = 0; etaIndex < probeBins; ++etaIndex) {
        const double eta =
            probeCoordinate(cell.etaMin, cell.etaMax, etaIndex, probeBins, probeRule);
        for (int phiIndex = 0; phiIndex < probeBins; ++phiIndex) {
            const double phi =
                probeCoordinate(cell.phiMin, cell.phiMax, phiIndex, probeBins, probeRule);
            const double weight =
                probeRule == ShapeUnweightingProbeRule::Interior
                    ? exactWeightAtEtaPhiInKnownCell(
                          current,
                          phaseSpace,
                          gridProposal,
                          cellIndex,
                          kT,
                          eta,
                          phi,
                          alpha)
                    : exactWeightAtEtaPhi(
                          current,
                          phaseSpace,
                          gridProposal,
                          pointCellIndex,
                          kT,
                          eta,
                          phi,
                          alpha);
            maxWeight = std::max(maxWeight, weight);
        }
    }
    return maxWeight;
}

[[nodiscard]] std::vector<double> probeCellMaxWeights(
    const LowCurrent& current,
    const PhotonPhaseSpace& phaseSpace,
    const GridAdaptiveEtaPhiProposal& gridProposal,
    const PointCellIndex& pointCellIndex,
    const std::vector<int>& cellProbeBins,
    ShapeUnweightingProbeRule probeRule,
    double alpha) {
    std::vector<double> cellMaxWeights(gridProposal.cells.size(), 0.0);
    for (std::size_t cellIndex = 0; cellIndex < gridProposal.cells.size(); ++cellIndex) {
        cellMaxWeights[cellIndex] =
            probeSingleCellMaxWeight(
                current,
                phaseSpace,
                gridProposal,
                pointCellIndex,
                cellIndex,
                cellProbeBins[cellIndex],
                probeRule,
                alpha);
    }
    return cellMaxWeights;
}

[[nodiscard]] ProbePlan makeProbePlan(
    const std::vector<ChargedParticle>& selectedParticles,
    const PhotonPhaseSpace& phaseSpace,
    const GridAdaptiveEtaPhiProposal& gridProposal,
    const ShapeUnweightingConfig& config,
    const std::vector<int>& baseProbeBins,
    const std::vector<double>& probedCellMaxWeights,
    const std::vector<double>& calibrationCellMaxWeights) {
    ProbePlan plan;
    plan.bins = baseProbeBins;
    plan.safetyFactors.assign(gridProposal.cells.size(), config.safetyFactor);

    if (config.evidenceTargetedBaseProbes &&
        effectiveBaseProbeBins(config) < config.probeBins) {
        for (std::size_t index = 0; index < gridProposal.cells.size(); ++index) {
            if (calibrationCellMaxWeights[index] <= probedCellMaxWeights[index]) {
                continue;
            }
            if (plan.bins[index] < config.probeBins) {
                plan.bins[index] = config.probeBins;
                ++plan.evidenceCellCount;
            }
        }
    }

    if (config.probePolicy == ShapeUnweightingProbePolicy::BoostedCells) {
        std::vector<bool> promote(gridProposal.cells.size(), false);
        for (std::size_t index = 0; index < gridProposal.cells.size(); ++index) {
            if (calibrationCellMaxWeights[index] > probedCellMaxWeights[index]) {
                promote[index] = true;
            }
        }

        for (const ChargedParticle& particle : selectedParticles) {
            double eta = 0.0;
            double phi = 0.0;
            if (!chargedParticleEtaPhi(particle, phaseSpace, eta, phi)) {
                continue;
            }
            if (eta < phaseSpace.etaMin || eta > phaseSpace.etaMax ||
                phi < phaseSpace.phiMin || phi >= phaseSpace.phiMax) {
                continue;
            }
            for (std::size_t index = 0; index < gridProposal.cells.size(); ++index) {
                if (cellNearCoordinate(
                        gridProposal.cells[index],
                        eta,
                        phi,
                        phaseSpace,
                        config.boostedCellNeighborScale)) {
                    promote[index] = true;
                }
            }
        }

        for (std::size_t index = 0; index < gridProposal.cells.size(); ++index) {
            if (!promote[index]) {
                continue;
            }
            plan.bins[index] = std::max(baseProbeBins[index], config.boostedProbeBins);
            plan.safetyFactors[index] =
                std::max(config.safetyFactor, config.boostedCellSafetyFactor);
            ++plan.promotedCellCount;
        }
    }

    for (std::size_t index = 0; index < plan.bins.size(); ++index) {
        const int bins = plan.bins[index];
        plan.probeEvaluations +=
            static_cast<std::size_t>(baseProbeBins[index]) *
            static_cast<std::size_t>(baseProbeBins[index]);
        if (bins != baseProbeBins[index]) {
            plan.probeEvaluations +=
                static_cast<std::size_t>(bins) * static_cast<std::size_t>(bins);
        }
    }
    return plan;
}

void validateShapeUnweightingConfig(const ShapeUnweightingConfig& config) {
    if (!validatePhotonPhaseSpace(config.phaseSpace).valid()) {
        throw std::invalid_argument("invalid photon phase space");
    }
    if (config.trials < 0) {
        throw std::invalid_argument("shape-unweighting trial count must be non-negative");
    }
    if (config.calibrationTrials < 0) {
        throw std::invalid_argument("shape-unweighting calibration trial count must be non-negative");
    }
    if (config.probeBins <= 0) {
        throw std::invalid_argument("shape-unweighting probeBins must be positive");
    }
    if (config.baseProbeBins < 0) {
        throw std::invalid_argument("shape-unweighting baseProbeBins must be non-negative");
    }
    if (!std::isfinite(config.targetedBaseProbeNeighborScale) ||
        config.targetedBaseProbeNeighborScale < 0.0) {
        throw std::invalid_argument(
            "shape-unweighting targetedBaseProbeNeighborScale must be non-negative");
    }
    if (config.boostedProbeBins <= 0) {
        throw std::invalid_argument("shape-unweighting boostedProbeBins must be positive");
    }
    if (!std::isfinite(config.safetyFactor) || config.safetyFactor <= 0.0) {
        throw std::invalid_argument("shape-unweighting safetyFactor must be positive");
    }
    if (!std::isfinite(config.boostedCellSafetyFactor) ||
        config.boostedCellSafetyFactor <= 0.0) {
        throw std::invalid_argument(
            "shape-unweighting boostedCellSafetyFactor must be positive");
    }
    if (!std::isfinite(config.boostedCellNeighborScale) ||
        config.boostedCellNeighborScale < 0.0) {
        throw std::invalid_argument(
            "shape-unweighting boostedCellNeighborScale must be non-negative");
    }
    if (config.targetedRefinementLevels < 0) {
        throw std::invalid_argument(
            "shape-unweighting targetedRefinementLevels must be non-negative");
    }
    if (!std::isfinite(config.targetedRefinementNeighborScale) ||
        config.targetedRefinementNeighborScale < 0.0) {
        throw std::invalid_argument(
            "shape-unweighting targetedRefinementNeighborScale must be non-negative");
    }
    validateGridAdaptiveEtaPhiProposalConfig(config.gridConfig);
}

}  // namespace

ShapeUnweightingResult generateShapeUnweightedPhotonCandidates(
    const std::vector<ChargedParticle>& selectedParticles,
    const ShapeUnweightingConfig& config,
    double alpha) {
    validateShapeUnweightingConfig(config);

    const LowCurrent current(selectedParticles);
    ShapeUnweightingResult result;
    result.chargeConservation = current.chargeConservation();
    result.gridBuildKT = std::sqrt(config.phaseSpace.kTMin * config.phaseSpace.kTMax);
    result.envelope.trialCount = config.trials;
    result.envelope.calibrationCount = config.calibrationTrials;
    result.envelope.safetyFactor = config.safetyFactor;
    result.envelope.boostedCellSafetyFactor = config.boostedCellSafetyFactor;
    result.envelope.boostedCellNeighborScale = config.boostedCellNeighborScale;
    result.envelope.targetedRefinementLevels = config.targetedRefinementLevels;
    result.envelope.targetedRefinementNeighborScale = config.targetedRefinementNeighborScale;
    result.envelope.probeBins = config.probeBins;
    result.envelope.baseProbeBins = effectiveBaseProbeBins(config);
    result.envelope.targetedBaseProbeNeighborScale = config.targetedBaseProbeNeighborScale;
    result.envelope.evidenceBaseProbeCellCount = 0;
    result.envelope.boostedProbeBins = config.boostedProbeBins;
    result.envelope.probeRule = config.probeRule;
    result.envelope.probePolicy = config.probePolicy;

    if (config.trials == 0) {
        result.status = ShapeUnweightingStatus::NoTrials;
        return result;
    }

    const auto gridStart = std::chrono::steady_clock::now();
    GridAdaptiveEtaPhiProposal gridProposal =
        makeGridAdaptiveEtaPhiProposal(
            current,
            config.phaseSpace,
            result.gridBuildKT,
            config.gridConfig);
    refineGridNearChargedParticles(
        current,
        selectedParticles,
        result.gridBuildKT,
        config.targetedRefinementLevels,
        config.targetedRefinementNeighborScale,
        gridProposal);
    result.envelope.gridCellCount = gridProposal.cells.size();
    const PointCellIndex pointCellIndex = makePointCellIndex(gridProposal);
    const std::vector<int> baseProbeBins =
        makeBaseProbeBins(selectedParticles, config.phaseSpace, gridProposal, config);
    result.envelope.fullBaseProbeCellCount =
        static_cast<int>(
            std::count(baseProbeBins.begin(), baseProbeBins.end(), config.probeBins));
    result.envelope.gridSeconds = secondsSince(gridStart);

    const auto baseProbeStart = std::chrono::steady_clock::now();
    const std::vector<double> probedCellMaxWeights =
        probeCellMaxWeights(
            current,
            config.phaseSpace,
            gridProposal,
            pointCellIndex,
            baseProbeBins,
            config.probeRule,
            alpha);
    result.envelope.baseProbeSeconds = secondsSince(baseProbeStart);

    std::vector<double> calibrationCellMaxWeights(gridProposal.cells.size(), 0.0);
    std::mt19937_64 calibrationRng(config.seed);
    const auto calibrationStart = std::chrono::steady_clock::now();
    for (int index = 0; index < config.calibrationTrials; ++index) {
        const TrialRecord record =
            makeTrialRecord(
                current,
                config.phaseSpace,
                gridProposal,
                pointCellIndex,
                calibrationRng,
                alpha);
        calibrationCellMaxWeights[record.cellIndex] = std::max(
            calibrationCellMaxWeights[record.cellIndex],
            record.sample.weightedPhoton.photon.weight);
    }
    result.envelope.calibrationSeconds = secondsSince(calibrationStart);

    const auto probePlanStart = std::chrono::steady_clock::now();
    for (std::size_t index = 0; index < gridProposal.cells.size(); ++index) {
        if (calibrationCellMaxWeights[index] > probedCellMaxWeights[index]) {
            ++result.envelope.suspectCellCount;
        }
    }

    ProbePlan probePlan =
        makeProbePlan(
            selectedParticles,
            config.phaseSpace,
            gridProposal,
            config,
            baseProbeBins,
            probedCellMaxWeights,
            calibrationCellMaxWeights);
    result.envelope.probePlanSeconds = secondsSince(probePlanStart);
    std::vector<double> finalProbedCellMaxWeights = probedCellMaxWeights;
    const auto boostedProbeStart = std::chrono::steady_clock::now();
    if (config.probePolicy == ShapeUnweightingProbePolicy::BoostedCells) {
        for (std::size_t index = 0; index < gridProposal.cells.size(); ++index) {
            if (probePlan.bins[index] == baseProbeBins[index]) {
                continue;
            }
            finalProbedCellMaxWeights[index] =
                probeSingleCellMaxWeight(
                    current,
                    config.phaseSpace,
                    gridProposal,
                    pointCellIndex,
                    index,
                    probePlan.bins[index],
                    config.probeRule,
                    alpha);
        }
    }
    result.envelope.boostedProbeSeconds = secondsSince(boostedProbeStart);
    result.envelope.evidenceBaseProbeCellCount = probePlan.evidenceCellCount;
    result.envelope.promotedCellCount = probePlan.promotedCellCount;
    result.envelope.probeEvaluations = probePlan.probeEvaluations;

    const auto envelopeBuildStart = std::chrono::steady_clock::now();
    const std::vector<EnvelopeCell> envelopeCells =
        makeEnvelopeCells(
            gridProposal,
            finalProbedCellMaxWeights,
            calibrationCellMaxWeights,
            probePlan.safetyFactors);
    const double envelopeMass = totalEnvelopeMass(envelopeCells);
    if (!(envelopeMass > 0.0) || !std::isfinite(envelopeMass)) {
        throw std::runtime_error("shape-unweighting envelope has zero total mass");
    }
    result.envelope.envelopeIntegral = envelopeMass;
    const EnvelopeCellSampler envelopeCellSampler = makeEnvelopeCellSampler(envelopeCells);
    result.envelope.envelopeBuildSeconds = secondsSince(envelopeBuildStart);

    result.candidates.reserve(static_cast<std::size_t>(config.trials));
    std::mt19937_64 trialRng(config.seed + 1);
    std::mt19937_64 acceptRng(config.seed + 2);
    double cappedAcceptanceSum = 0.0;

    const auto trialStart = std::chrono::steady_clock::now();
    for (int index = 0; index < config.trials; ++index) {
        const TrialRecord record =
            makeEnvelopeTrialRecord(
                current,
                config.phaseSpace,
                gridProposal,
                envelopeCellSampler,
                pointCellIndex,
                envelopeMass,
                trialRng,
                alpha);
        const double weight = record.sample.weightedPhoton.photon.weight;
        double bound = envelopeCells[record.cellIndex].bound;
        if (!(bound > 0.0) || !std::isfinite(bound)) {
            ++result.envelope.fallbackCount;
            bound = std::numeric_limits<double>::infinity();
        }

        const double ratio = weight / bound;
        if (std::isfinite(ratio)) {
            result.envelope.maxTrialBoundRatio =
                std::max(result.envelope.maxTrialBoundRatio, ratio);
        }
        if (weight > bound) {
            ++result.envelope.failedEnvelopeCount;
        }

        const double acceptanceProbability =
            std::isfinite(bound) ? std::min(1.0, ratio) : 0.0;
        cappedAcceptanceSum += acceptanceProbability;

        const double unit = std::generate_canonical<double, 53>(acceptRng);
        if (unit < acceptanceProbability) {
            Photon photon = record.sample.weightedPhoton.photon;
            photon.weight = 1.0;
            result.candidates.push_back(ShapeUnweightedPhotonCandidate{
                photon,
                record.sample.kTProposal.kT,
                record.sample.etaPhiProposal.eta,
                record.sample.etaPhiProposal.phi,
            });
        }
    }
    result.envelope.trialSeconds = secondsSince(trialStart);

    const auto finalizationStart = std::chrono::steady_clock::now();
    result.envelope.acceptedCount = static_cast<int>(result.candidates.size());
    result.envelope.observedAcceptanceEfficiency =
        static_cast<double>(result.envelope.acceptedCount) / static_cast<double>(config.trials);
    result.envelope.expectedAcceptanceEfficiency =
        cappedAcceptanceSum / static_cast<double>(config.trials);

    if (result.envelope.failedEnvelopeCount > 0) {
        result.status = ShapeUnweightingStatus::EnvelopeFailed;
        if (config.failurePolicy == EnvelopeFailurePolicy::Strict) {
            result.candidates.clear();
            result.envelope.acceptedCount = 0;
            result.envelope.observedAcceptanceEfficiency = 0.0;
        }
    } else {
        result.status = ShapeUnweightingStatus::Ok;
    }
    result.envelope.finalizationSeconds = secondsSince(finalizationStart);

    return result;
}

ShapeUnweightingIndexConsistencyDiagnostic validateShapeUnweightingIndexConsistency(
    const std::vector<ChargedParticle>& selectedParticles,
    const ShapeUnweightingConfig& config,
    double alpha) {
    validateShapeUnweightingConfig(config);
    if (!std::isfinite(alpha) || alpha <= 0.0) {
        throw std::invalid_argument("fine-structure constant must be positive and finite");
    }

    const LowCurrent current(selectedParticles);
    const double gridBuildKT = std::sqrt(config.phaseSpace.kTMin * config.phaseSpace.kTMax);
    GridAdaptiveEtaPhiProposal gridProposal =
        makeGridAdaptiveEtaPhiProposal(
            current,
            config.phaseSpace,
            gridBuildKT,
            config.gridConfig);
    refineGridNearChargedParticles(
        current,
        selectedParticles,
        gridBuildKT,
        config.targetedRefinementLevels,
        config.targetedRefinementNeighborScale,
        gridProposal);
    const PointCellIndex pointCellIndex = makePointCellIndex(gridProposal);

    ShapeUnweightingIndexConsistencyDiagnostic diagnostic;
    diagnostic.consistent = true;
    diagnostic.gridCellCount = gridProposal.cells.size();

    auto checkPoint = [&](double eta, double phi) {
        const std::size_t reference = findCellIndex(gridProposal, eta, phi);
        const std::size_t indexed = findCellIndex(gridProposal, pointCellIndex, eta, phi);
        ++diagnostic.pointLookupChecks;
        if (reference != indexed) {
            diagnostic.consistent = false;
        }
    };

    for (const GridAdaptiveEtaPhiCell& cell : gridProposal.cells) {
        checkPoint(0.5 * (cell.etaMin + cell.etaMax), 0.5 * (cell.phiMin + cell.phiMax));
        checkPoint(std::nextafter(cell.etaMin, cell.etaMax), 0.5 * (cell.phiMin + cell.phiMax));
        checkPoint(std::nextafter(cell.etaMax, cell.etaMin), 0.5 * (cell.phiMin + cell.phiMax));
        checkPoint(0.5 * (cell.etaMin + cell.etaMax), std::nextafter(cell.phiMin, cell.phiMax));
        checkPoint(0.5 * (cell.etaMin + cell.etaMax), std::nextafter(cell.phiMax, cell.phiMin));
    }

    const double probeKT = std::sqrt(config.phaseSpace.kTMin * config.phaseSpace.kTMax);
    const int weightProbeBins = std::min(config.probeBins, 3);
    for (std::size_t cellIndex = 0; cellIndex < gridProposal.cells.size(); ++cellIndex) {
        const GridAdaptiveEtaPhiCell& cell = gridProposal.cells[cellIndex];
        for (int etaIndex = 0; etaIndex < weightProbeBins; ++etaIndex) {
            const double eta = probeCoordinate(
                cell.etaMin,
                cell.etaMax,
                etaIndex,
                weightProbeBins,
                ShapeUnweightingProbeRule::Interior);
            for (int phiIndex = 0; phiIndex < weightProbeBins; ++phiIndex) {
                const double phi = probeCoordinate(
                    cell.phiMin,
                    cell.phiMax,
                    phiIndex,
                    weightProbeBins,
                    ShapeUnweightingProbeRule::Interior);
                const double reference = exactWeightAtEtaPhi(
                    current,
                    config.phaseSpace,
                    gridProposal,
                    pointCellIndex,
                    probeKT,
                    eta,
                    phi,
                    alpha);
                const double knownCell = exactWeightAtEtaPhiInKnownCell(
                    current,
                    config.phaseSpace,
                    gridProposal,
                    cellIndex,
                    probeKT,
                    eta,
                    phi,
                    alpha);
                ++diagnostic.knownCellProbeWeightChecks;
                const double scale = std::max({1.0, std::abs(reference), std::abs(knownCell)});
                if (std::abs(reference - knownCell) > 1.0e-12 * scale) {
                    diagnostic.consistent = false;
                }
            }
        }
    }

    std::vector<EnvelopeCell> envelopeCells;
    envelopeCells.reserve(gridProposal.cells.size());
    for (std::size_t index = 0; index < gridProposal.cells.size(); ++index) {
        const double envelopeMass =
            static_cast<double>((index % 17) + 1) *
            cellProposalMass(gridProposal, gridProposal.cells[index]);
        envelopeCells.push_back(EnvelopeCell{1.0, envelopeMass, envelopeMass});
    }
    const double envelopeMass = totalEnvelopeMass(envelopeCells);
    const EnvelopeCellSampler sampler = makeEnvelopeCellSampler(envelopeCells);

    auto checkEnvelope = [&](double unitCell) {
        const std::size_t reference = sampleEnvelopeCell(envelopeCells, envelopeMass, unitCell);
        const std::size_t indexed = sampleEnvelopeCell(sampler, envelopeMass, unitCell);
        ++diagnostic.envelopeSamplingChecks;
        if (reference != indexed) {
            diagnostic.consistent = false;
        }
    };

    checkEnvelope(0.0);
    checkEnvelope(std::nextafter(1.0, 0.0));
    double cumulative = 0.0;
    for (std::size_t index = 0; index < std::min<std::size_t>(gridProposal.cells.size(), 128);
         ++index) {
        cumulative += envelopeCells[index].envelopeMass;
        const double unit = cumulative / envelopeMass;
        if (unit > 0.0 && unit < 1.0) {
            checkEnvelope(std::nextafter(unit, 0.0));
            checkEnvelope(unit);
            checkEnvelope(std::nextafter(unit, 1.0));
        }
    }

    return diagnostic;
}

}  // namespace lowmc
