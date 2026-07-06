#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

#include "lowmc/Photon.h"
#include "lowmc/WeightedPhotonSampler.h"

namespace lowmc {

struct WeightSummary {
    std::size_t count = 0;
    double sumWeight = 0.0;
    double minWeight = 0.0;
    double maxWeight = 0.0;
    double meanWeight = 0.0;
    double meanSquaredWeight = 0.0;
    double varianceWeight = 0.0;
    double effectiveSampleSize = 0.0;
};

[[nodiscard]] inline WeightSummary summarizeWeights(const std::vector<double>& weights) {
    WeightSummary summary;
    summary.count = weights.size();

    if (weights.empty()) {
        return summary;
    }

    summary.minWeight = std::numeric_limits<double>::infinity();
    summary.maxWeight = 0.0;

    double sumSquaredWeight = 0.0;
    for (double weight : weights) {
        if (!std::isfinite(weight) || weight < 0.0) {
            throw std::invalid_argument("weights must be finite and non-negative");
        }

        summary.sumWeight += weight;
        sumSquaredWeight += weight * weight;
        summary.minWeight = std::min(summary.minWeight, weight);
        summary.maxWeight = std::max(summary.maxWeight, weight);
    }

    summary.meanWeight = summary.sumWeight / static_cast<double>(summary.count);
    summary.meanSquaredWeight = sumSquaredWeight / static_cast<double>(summary.count);
    summary.varianceWeight =
        std::max(0.0, summary.meanSquaredWeight - summary.meanWeight * summary.meanWeight);

    if (sumSquaredWeight > 0.0) {
        summary.effectiveSampleSize = summary.sumWeight * summary.sumWeight / sumSquaredWeight;
    }

    return summary;
}

[[nodiscard]] inline WeightSummary summarizePhotons(const std::vector<Photon>& photons) {
    std::vector<double> weights;
    weights.reserve(photons.size());

    for (const Photon& photon : photons) {
        weights.push_back(photon.weight);
    }

    return summarizeWeights(weights);
}

[[nodiscard]] inline WeightSummary summarizeWeightedPhotonSamples(
    const std::vector<WeightedPhotonSample>& samples) {
    std::vector<double> weights;
    weights.reserve(samples.size());

    for (const WeightedPhotonSample& sample : samples) {
        weights.push_back(sample.photon.weight);
    }

    return summarizeWeights(weights);
}

}  // namespace lowmc
