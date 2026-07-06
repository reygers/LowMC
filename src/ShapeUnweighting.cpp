#include "lowmc/ShapeUnweighting.h"

#include <algorithm>
#include <cmath>
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

struct ProbePlan {
    std::vector<int> bins;
    std::vector<double> safetyFactors;
    int promotedCellCount = 0;
    std::size_t probeEvaluations = 0;
};

[[nodiscard]] double cellArea(const GridAdaptiveEtaPhiCell& cell) {
    return (cell.etaMax - cell.etaMin) * (cell.phiMax - cell.phiMin);
}

[[nodiscard]] double totalEtaPhiArea(const PhotonPhaseSpace& phaseSpace) {
    return (phaseSpace.etaMax - phaseSpace.etaMin) *
           (phaseSpace.phiMax - phaseSpace.phiMin);
}

[[nodiscard]] double cellProposalMass(
    const GridAdaptiveEtaPhiProposal& proposal,
    const GridAdaptiveEtaPhiCell& cell) {
    return proposal.config.uniformFraction *
               cellArea(cell) / totalEtaPhiArea(proposal.phaseSpace) +
           cell.mixtureWeight;
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

[[nodiscard]] TrialRecord makeTrialRecord(
    const LowCurrent& current,
    const PhotonPhaseSpace& phaseSpace,
    const GridAdaptiveEtaPhiProposal& gridProposal,
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
        findCellIndex(gridProposal, sample.etaPhiProposal.eta, sample.etaPhiProposal.phi);
    return TrialRecord{sample, cellIndex};
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

[[nodiscard]] TrialRecord makeEnvelopeTrialRecord(
    const LowCurrent& current,
    const PhotonPhaseSpace& phaseSpace,
    const GridAdaptiveEtaPhiProposal& gridProposal,
    const std::vector<EnvelopeCell>& envelopeCells,
    double totalEnvelopeMass,
    std::mt19937_64& rng,
    double alpha) {
    const double unitKT = std::generate_canonical<double, 53>(rng);
    const double unitCell = std::generate_canonical<double, 53>(rng);
    const double unitEta = std::generate_canonical<double, 53>(rng);
    const double unitPhi = std::generate_canonical<double, 53>(rng);

    const std::size_t cellIndex =
        sampleEnvelopeCell(envelopeCells, totalEnvelopeMass, unitCell);
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

    const double etaPhiDensity = gridAdaptiveEtaPhiDensity(gridProposal, eta, phi);
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
    double kT,
    double eta,
    double phi,
    double alpha) {
    const double proposalDensity =
        logUniformKTDensityPhysicalMeasure(phaseSpace, kT) *
        gridAdaptiveEtaPhiDensity(gridProposal, eta, phi);
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
    const GridAdaptiveEtaPhiCell& cell,
    int probeBins,
    ShapeUnweightingProbeRule probeRule,
    double alpha) {
    const double kT = std::sqrt(phaseSpace.kTMin * phaseSpace.kTMax);
    double maxWeight = 0.0;
    for (int etaIndex = 0; etaIndex < probeBins; ++etaIndex) {
        const double eta =
            probeCoordinate(cell.etaMin, cell.etaMax, etaIndex, probeBins, probeRule);
        for (int phiIndex = 0; phiIndex < probeBins; ++phiIndex) {
            const double phi =
                probeCoordinate(cell.phiMin, cell.phiMax, phiIndex, probeBins, probeRule);
            maxWeight = std::max(
                maxWeight,
                exactWeightAtEtaPhi(current, phaseSpace, gridProposal, kT, eta, phi, alpha));
        }
    }
    return maxWeight;
}

[[nodiscard]] std::vector<double> probeCellMaxWeights(
    const LowCurrent& current,
    const PhotonPhaseSpace& phaseSpace,
    const GridAdaptiveEtaPhiProposal& gridProposal,
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
                gridProposal.cells[cellIndex],
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
    const std::vector<double>& probedCellMaxWeights,
    const std::vector<double>& calibrationCellMaxWeights) {
    ProbePlan plan;
    plan.bins.assign(gridProposal.cells.size(), config.probeBins);
    plan.safetyFactors.assign(gridProposal.cells.size(), config.safetyFactor);

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
            plan.bins[index] = std::max(config.probeBins, config.boostedProbeBins);
            plan.safetyFactors[index] =
                std::max(config.safetyFactor, config.boostedCellSafetyFactor);
            ++plan.promotedCellCount;
        }
    }

    for (int bins : plan.bins) {
        plan.probeEvaluations +=
            static_cast<std::size_t>(config.probeBins) *
            static_cast<std::size_t>(config.probeBins);
        if (bins != config.probeBins) {
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
    result.envelope.boostedProbeBins = config.boostedProbeBins;
    result.envelope.probeRule = config.probeRule;
    result.envelope.probePolicy = config.probePolicy;

    if (config.trials == 0) {
        result.status = ShapeUnweightingStatus::NoTrials;
        return result;
    }

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
    const std::vector<int> baseProbeBins(gridProposal.cells.size(), config.probeBins);

    const std::vector<double> probedCellMaxWeights =
        probeCellMaxWeights(
            current,
            config.phaseSpace,
            gridProposal,
            baseProbeBins,
            config.probeRule,
            alpha);

    std::vector<double> calibrationCellMaxWeights(gridProposal.cells.size(), 0.0);
    std::mt19937_64 calibrationRng(config.seed);
    for (int index = 0; index < config.calibrationTrials; ++index) {
        const TrialRecord record =
            makeTrialRecord(current, config.phaseSpace, gridProposal, calibrationRng, alpha);
        calibrationCellMaxWeights[record.cellIndex] = std::max(
            calibrationCellMaxWeights[record.cellIndex],
            record.sample.weightedPhoton.photon.weight);
    }

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
            probedCellMaxWeights,
            calibrationCellMaxWeights);
    std::vector<double> finalProbedCellMaxWeights = probedCellMaxWeights;
    if (config.probePolicy == ShapeUnweightingProbePolicy::BoostedCells) {
        for (std::size_t index = 0; index < gridProposal.cells.size(); ++index) {
            if (probePlan.bins[index] == config.probeBins) {
                continue;
            }
            finalProbedCellMaxWeights[index] =
                probeSingleCellMaxWeight(
                    current,
                    config.phaseSpace,
                    gridProposal,
                    gridProposal.cells[index],
                    probePlan.bins[index],
                    config.probeRule,
                    alpha);
        }
    }
    result.envelope.promotedCellCount = probePlan.promotedCellCount;
    result.envelope.probeEvaluations = probePlan.probeEvaluations;

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

    result.candidates.reserve(static_cast<std::size_t>(config.trials));
    std::mt19937_64 trialRng(config.seed + 1);
    std::mt19937_64 acceptRng(config.seed + 2);
    double cappedAcceptanceSum = 0.0;

    for (int index = 0; index < config.trials; ++index) {
        const TrialRecord record =
            makeEnvelopeTrialRecord(
                current,
                config.phaseSpace,
                gridProposal,
                envelopeCells,
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

    return result;
}

}  // namespace lowmc
