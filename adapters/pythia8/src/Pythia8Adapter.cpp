#include "lowmc/pythia8/Pythia8Adapter.h"

#include <Pythia8/Event.h>

#include <string>
#include <vector>

namespace {

[[nodiscard]] bool isConventionalBeamIndex(int index) {
    return index == 1 || index == 2;
}

[[nodiscard]] std::string sourceLabel(int index, const Pythia8::Particle& particle) {
    return "pythia8:index=" + std::to_string(index) +
           ",id=" + std::to_string(particle.id()) +
           ",status=" + std::to_string(particle.status());
}

}  // namespace

namespace lowmc::pythia8 {

const char* pythia8AdapterBuildMarker() {
    return "lowmc_pythia8_adapter";
}

AdapterParticleCandidate makeAdapterParticleCandidate(
    const Pythia8::Particle& particle,
    const Pythia8ParticleClassification& classification) {
    return AdapterParticleCandidate{
        FourVector{particle.e(), particle.px(), particle.py(), particle.pz()},
        particle.id(),
        classification.isIncomingBeam,
        classification.isFinalStable,
        classification.isBeamRemnant,
        classification.sourceLabel,
    };
}

std::vector<AdapterParticleCandidate> collectAdapterParticleCandidates(
    const Pythia8::Event& event) {
    std::vector<AdapterParticleCandidate> candidates;

    for (int index = 0; index < event.size(); ++index) {
        const Pythia8::Particle& particle = event[index];
        const bool isIncomingBeam = isConventionalBeamIndex(index);
        const bool isFinalStable = !isIncomingBeam && particle.isFinal();

        if (!isIncomingBeam && !isFinalStable) {
            continue;
        }

        candidates.push_back(makeAdapterParticleCandidate(
            particle,
            Pythia8ParticleClassification{
                isIncomingBeam,
                isFinalStable,
                false,
                sourceLabel(index, particle),
            }));
    }

    return candidates;
}

AdapterParticleSelectionResult selectChargedParticles(
    const Pythia8::Event& event,
    const ParticleSelectionConfig& config) {
    return selectAdapterParticles(collectAdapterParticleCandidates(event), config);
}

}  // namespace lowmc::pythia8

