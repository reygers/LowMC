#pragma once

#include <cstdint>
#include <vector>

#include "lowmc/ChargedParticle.h"
#include "lowmc/LowCurrent.h"
#include "lowmc/SamplerDiagnostics.h"

namespace lowmc {

inline constexpr double ElectronMassGeV = 0.00051099895000;

// The generated density is defined with respect to
// dM^2 qT dqT dy dphi.  Rapidity (not pseudorapidity) is used because the
// virtual photon is massive.
struct DielectronPhaseSpace {
    double massMin = 0.0;
    double massMax = 0.0;
    double qTMin = 0.0;
    double qTMax = 0.0;
    double rapidityMin = 0.0;
    double rapidityMax = 0.0;
    double phiMin = 0.0;
    double phiMax = 0.0;
};

enum class DielectronPhaseSpaceValidationError {
    None,
    NonFiniteBound,
    BelowPairThreshold,
    EmptyMassRange,
    NonPositiveQTMin,
    EmptyQTRange,
    EmptyRapidityRange,
    EmptyPhiRange,
};

struct DielectronPhaseSpaceValidation {
    DielectronPhaseSpaceValidationError error =
        DielectronPhaseSpaceValidationError::None;

    [[nodiscard]] constexpr bool valid() const {
        return error == DielectronPhaseSpaceValidationError::None;
    }
};

struct Dielectron {
    FourVector electron;
    FourVector positron;
    double weight = 1.0;

    [[nodiscard]] constexpr FourVector pairMomentum() const {
        return electron + positron;
    }
};

struct WeightedDielectronSample {
    Dielectron dielectron;
    FourVector virtualPhoton;
    double mass = 0.0;
    double qT = 0.0;
    double rapidity = 0.0;
    double phi = 0.0;
    // Marginal proposal density with respect to dM^2 qT dqT dy dphi.
    // The isotropic angular proposal cancels inside decayAngularWeight.
    double proposalDensity = 0.0;
    double productionDensity = 0.0;
    double conversionDensityDMass2 = 0.0;
    double decayAngularWeight = 1.0;
};

struct WeightedDielectronGeneratorConfig {
    DielectronPhaseSpace phaseSpace;
    int samples = 1;
    std::uint64_t seed = 12345;
    int oversamplingFactor = 1;
};

struct WeightedDielectronGenerationResult {
    std::vector<WeightedDielectronSample> samples;
    WeightSummary weightSummary;
    ChargeConservationDiagnostic chargeConservation;
};

[[nodiscard]] DielectronPhaseSpaceValidation validateDielectronPhaseSpace(
    const DielectronPhaseSpace& phaseSpace,
    double electronMass = ElectronMassGeV);

[[nodiscard]] FourVector virtualPhotonMomentumFromMassQTRapidityPhi(
    double mass,
    double qT,
    double rapidity,
    double phi);

// dP(gamma* -> e+e-) / dM^2 in the angle-integrated pointlike
// Kroll-Wada approximation.
[[nodiscard]] double krollWadaConversionDensityDMass2(
    double mass,
    double electronMass = ElectronMassGeV,
    double alpha = FineStructureConstant);

[[nodiscard]] double lowVirtualPhotonDensity(
    const LowCurrent& current,
    const FourVector& virtualPhotonMomentum,
    double alpha = FineStructureConstant);

// Unit inputs are independent uniform values in [0, 1).  Mass and qT are
// sampled logarithmically.  Decay angles are proposed isotropically and the
// returned weight contains the leading-current polarization correction.
[[nodiscard]] WeightedDielectronSample makeWeightedDielectronSample(
    const LowCurrent& current,
    const DielectronPhaseSpace& phaseSpace,
    double unitMass,
    double unitQT,
    double unitRapidity,
    double unitPhi,
    double unitDecayCosTheta,
    double unitDecayPhi,
    double electronMass = ElectronMassGeV,
    double alpha = FineStructureConstant);

[[nodiscard]] WeightedDielectronGenerationResult generateWeightedDielectrons(
    const std::vector<ChargedParticle>& selectedParticles,
    const WeightedDielectronGeneratorConfig& config,
    double electronMass = ElectronMassGeV,
    double alpha = FineStructureConstant);

}  // namespace lowmc
