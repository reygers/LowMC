#pragma once

#include <cstddef>
#include <vector>

#include "lowmc/ChargedParticle.h"
#include "lowmc/PhotonPhaseSpace.h"

namespace lowmc {

struct AnalyticLowLegKinematics {
    double massSquared = 0.0;
    double mass = 0.0;
    double transverseMomentum = 0.0;
    double transverseMass = 0.0;
    double rapidity = 0.0;
    double phi = 0.0;
    double globalMinimumDenominatorScale = 0.0;
    bool valid = false;
};

struct AnalyticLowAngularBound {
    double value = 0.0;
    int invalidParticleCount = 0;
    bool valid = false;
};

class AnalyticLowBound {
public:
    explicit AnalyticLowBound(std::vector<ChargedParticle> particles);

    [[nodiscard]] const std::vector<ChargedParticle>& particles() const;
    [[nodiscard]] const std::vector<AnalyticLowLegKinematics>& legKinematics() const;
    [[nodiscard]] int invalidParticleCount() const;
    [[nodiscard]] bool valid() const;

    [[nodiscard]] double denominatorScale(
        std::size_t particleIndex,
        double eta,
        double phi) const;

    [[nodiscard]] double angularIntensity(double eta, double phi) const;

    [[nodiscard]] AnalyticLowAngularBound globalAngularIntensityBound() const;

    [[nodiscard]] AnalyticLowAngularBound cellAngularIntensityBound(
        double etaMin,
        double etaMax,
        double phiMin,
        double phiMax) const;

private:
    std::vector<ChargedParticle> particles_;
    std::vector<AnalyticLowLegKinematics> legs_;
    std::vector<double> signedPairCoefficients_;
    std::vector<double> absolutePairCoefficients_;
    int invalidParticleCount_ = 0;
};

[[nodiscard]] double analyticLowWeightBound(
    const AnalyticLowAngularBound& angularBound,
    const PhotonPhaseSpace& phaseSpace,
    double etaPhiProposalDensity,
    double alpha);

}  // namespace lowmc
