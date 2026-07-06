#pragma once

#include <string>
#include <vector>

#include <HepMC3/Units.h>

#include "lowmc/AdapterSupport.h"
#include "lowmc/ShapeUnweighting.h"

namespace HepMC3 {
class GenEvent;
class GenParticle;
}  // namespace HepMC3

namespace lowmc::hepmc3 {

struct HepMC3ParticleClassification {
    bool isIncomingBeam = false;
    bool isFinalStable = false;
    bool isBeamRemnant = false;
    std::string sourceLabel;
};

struct HepMC3LowPhotonOutputOptions {
    double enhancementFactor = 1.0;
    std::string version;
    double lambdaEvent = 0.0;
    bool hasLambdaEvent = false;
    std::string validationClass;
    std::string shapeStatus;
    double maxTrialBoundRatio = 0.0;
    bool hasMaxTrialBoundRatio = false;
};

struct HepMC3LowPhotonOutputResult {
    int photonsAdded = 0;
};

[[nodiscard]] const char* hepmc3AdapterBuildMarker();

[[nodiscard]] AdapterParticleCandidate makeAdapterParticleCandidate(
    const HepMC3::GenParticle& particle,
    HepMC3::Units::MomentumUnit momentumUnit,
    const HepMC3ParticleClassification& classification = {});

[[nodiscard]] std::vector<AdapterParticleCandidate> collectAdapterParticleCandidates(
    const HepMC3::GenEvent& event);

[[nodiscard]] AdapterParticleSelectionResult selectChargedParticles(
    const HepMC3::GenEvent& event,
    const ParticleSelectionConfig& config = {});

[[nodiscard]] HepMC3LowPhotonOutputResult appendUnweightedLowPhotons(
    HepMC3::GenEvent& event,
    const std::vector<ShapeUnweightedPhotonCandidate>& photons,
    const HepMC3LowPhotonOutputOptions& options = {});

}  // namespace lowmc::hepmc3
