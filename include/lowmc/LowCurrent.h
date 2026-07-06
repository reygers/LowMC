#pragma once

#include <vector>

#include "lowmc/ChargedParticle.h"
#include "lowmc/FourVector.h"

namespace lowmc {

inline constexpr double FineStructureConstant = 1.0 / 137.035999084;

struct ChargeConservationDiagnostic {
    double netCharge = 0.0;
    bool conserved = true;
};

class LowCurrent {
public:
    explicit LowCurrent(std::vector<ChargedParticle> particles);

    [[nodiscard]] const std::vector<ChargedParticle>& particles() const;

    [[nodiscard]] FourVector current(const FourVector& photonMomentum) const;
    [[nodiscard]] double intensity(const FourVector& photonMomentum) const;
    [[nodiscard]] double intensityPairSum(const FourVector& photonMomentum) const;

    [[nodiscard]] ChargeConservationDiagnostic chargeConservation(
        double tolerance = 1.0e-12) const;

private:
    std::vector<ChargedParticle> particles_;
};

[[nodiscard]] double lowPhotonDensity(
    const LowCurrent& current,
    const FourVector& photonMomentum,
    double alpha = FineStructureConstant);

}  // namespace lowmc
