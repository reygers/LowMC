#pragma once

#include "lowmc/Dielectron.h"

namespace lowmc {

struct BalekParameters {
    // The returned cross sections use the same unit as sigmaInelastic.
    double sigmaInelastic = 57.8;
    double meanChargeTransferSquared = 1.32;
    double electronMass = ElectronMassGeV;
    double alpha = FineStructureConstant;
};

// BPP Eq. (9a), including the exact lepton-mass factor from their Eq. (3):
// d^3 sigma / (dy dqT dM).
[[nodiscard]] double balekDifferentialCrossSection(
    double mass,
    double qT,
    const BalekParameters& parameters = {});

// Integral of balekDifferentialCrossSection over qT.  This is the mass
// projection used for the left panel of ALICE Fig. 4.
[[nodiscard]] double balekMassProjection(
    double mass,
    double qTMin,
    double qTMax,
    const BalekParameters& parameters = {});

// Analytic massless-electron integral over M.  Above the ALICE threshold
// (M >= 0.15 GeV) the omitted electron-mass correction is negligible.
[[nodiscard]] double balekQTProjectionMassless(
    double qT,
    double massMin,
    double massMax,
    const BalekParameters& parameters = {});

}  // namespace lowmc
