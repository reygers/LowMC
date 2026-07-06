#pragma once

#include <cmath>
#include <stdexcept>

#include "lowmc/PhotonPhaseSpace.h"

namespace lowmc {

struct LogUniformKTProposalPoint {
    double kT = 0.0;
    double densityDKT = 0.0;
    double densityPhysicalMeasure = 0.0;
};

inline void requireValidKTProposalRange(const PhotonPhaseSpace& phaseSpace) {
    if (!validatePhotonPhaseSpace(phaseSpace).valid()) {
        throw std::invalid_argument("invalid photon phase space");
    }
}

[[nodiscard]] inline double logUniformKTDensityDKT(
    const PhotonPhaseSpace& phaseSpace,
    double kT) {
    requireValidKTProposalRange(phaseSpace);

    if (!std::isfinite(kT) || kT < phaseSpace.kTMin || kT > phaseSpace.kTMax) {
        throw std::invalid_argument("kT must be finite and inside the phase-space range");
    }

    return 1.0 / (kT * std::log(phaseSpace.kTMax / phaseSpace.kTMin));
}

[[nodiscard]] inline double logUniformKTDensityPhysicalMeasure(
    const PhotonPhaseSpace& phaseSpace,
    double kT) {
    return logUniformKTDensityDKT(phaseSpace, kT) / kT;
}

[[nodiscard]] inline LogUniformKTProposalPoint sampleLogUniformKT(
    const PhotonPhaseSpace& phaseSpace,
    double unit) {
    requireValidKTProposalRange(phaseSpace);

    if (!std::isfinite(unit) || unit < 0.0 || unit >= 1.0) {
        throw std::invalid_argument("unit random value must satisfy 0 <= u < 1");
    }

    const double logMin = std::log(phaseSpace.kTMin);
    const double logMax = std::log(phaseSpace.kTMax);
    double kT = std::exp(logMin + unit * (logMax - logMin));
    if (kT >= phaseSpace.kTMax) {
        kT = std::nextafter(phaseSpace.kTMax, phaseSpace.kTMin);
    }

    return LogUniformKTProposalPoint{
        kT,
        logUniformKTDensityDKT(phaseSpace, kT),
        logUniformKTDensityPhysicalMeasure(phaseSpace, kT),
    };
}

}  // namespace lowmc
