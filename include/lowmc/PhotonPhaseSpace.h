#pragma once

#include <cmath>

namespace lowmc {

struct PhotonPhaseSpace {
    double kTMin = 0.0;
    double kTMax = 0.0;
    double etaMin = 0.0;
    double etaMax = 0.0;
    double phiMin = 0.0;
    double phiMax = 0.0;
};

enum class PhotonPhaseSpaceValidationError {
    None,
    NonFiniteBound,
    NonPositiveKTMin,
    EmptyKTRange,
    EmptyEtaRange,
    EmptyPhiRange,
};

struct PhotonPhaseSpaceValidation {
    PhotonPhaseSpaceValidationError error = PhotonPhaseSpaceValidationError::None;

    [[nodiscard]] constexpr bool valid() const {
        return error == PhotonPhaseSpaceValidationError::None;
    }
};

[[nodiscard]] inline PhotonPhaseSpaceValidation
validatePhotonPhaseSpace(const PhotonPhaseSpace& phaseSpace) {
    if (!std::isfinite(phaseSpace.kTMin) || !std::isfinite(phaseSpace.kTMax) ||
        !std::isfinite(phaseSpace.etaMin) || !std::isfinite(phaseSpace.etaMax) ||
        !std::isfinite(phaseSpace.phiMin) || !std::isfinite(phaseSpace.phiMax)) {
        return {PhotonPhaseSpaceValidationError::NonFiniteBound};
    }

    if (phaseSpace.kTMin <= 0.0) {
        return {PhotonPhaseSpaceValidationError::NonPositiveKTMin};
    }

    if (phaseSpace.kTMax <= phaseSpace.kTMin) {
        return {PhotonPhaseSpaceValidationError::EmptyKTRange};
    }

    if (phaseSpace.etaMax <= phaseSpace.etaMin) {
        return {PhotonPhaseSpaceValidationError::EmptyEtaRange};
    }

    if (phaseSpace.phiMax <= phaseSpace.phiMin) {
        return {PhotonPhaseSpaceValidationError::EmptyPhiRange};
    }

    return {};
}

}  // namespace lowmc
