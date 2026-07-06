#pragma once

#include <string>
#include <vector>

#include "lowmc/AdapterSupport.h"

namespace Pythia8 {
class Event;
class Particle;
}  // namespace Pythia8

namespace lowmc::pythia8 {

struct Pythia8ParticleClassification {
    bool isIncomingBeam = false;
    bool isFinalStable = false;
    bool isBeamRemnant = false;
    std::string sourceLabel;
};

[[nodiscard]] const char* pythia8AdapterBuildMarker();

[[nodiscard]] AdapterParticleCandidate makeAdapterParticleCandidate(
    const Pythia8::Particle& particle,
    const Pythia8ParticleClassification& classification = {});

[[nodiscard]] std::vector<AdapterParticleCandidate> collectAdapterParticleCandidates(
    const Pythia8::Event& event);

[[nodiscard]] AdapterParticleSelectionResult selectChargedParticles(
    const Pythia8::Event& event,
    const ParticleSelectionConfig& config = {});

}  // namespace lowmc::pythia8

