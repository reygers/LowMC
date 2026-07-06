#include "lowmc/pythia8/Pythia8Adapter.h"

#include <Pythia8/Event.h>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "FAILED: " << message << '\n';
        std::exit(EXIT_FAILURE);
    }
}

Pythia8::Event makeToyEvent() {
    Pythia8::Event event;
    event.reset();
    event.append(2212, -12, 0, 0, 0, 0, 0, 0, 0.0, 0.0, 6500.0, 6500.0, 0.938);
    event.append(2212, -12, 0, 0, 0, 0, 0, 0, 0.0, 0.0, -6500.0, 6500.0, 0.938);
    event.append(2212, 1, 0, 0, 0, 0, 0, 0, 1.0, 0.0, 0.0, 1.37, 0.938);
    event.append(2212, 1, 0, 0, 0, 0, 0, 0, -1.0, 0.0, 0.0, 1.37, 0.938);
    event.append(22, 1, 0, 0, 0, 0, 0, 0, 0.0, 0.5, 0.0, 0.5, 0.0);
    event.append(321, -91, 0, 0, 0, 0, 0, 0, 0.0, -0.5, 0.0, 0.6, 0.494);
    return event;
}

void testBuildMarker() {
    require(std::string(lowmc::pythia8::pythia8AdapterBuildMarker()) ==
                "lowmc_pythia8_adapter",
            "unexpected build marker");
}

void testCandidateConversion() {
    Pythia8::Event event = makeToyEvent();
    const lowmc::AdapterParticleCandidate candidate =
        lowmc::pythia8::makeAdapterParticleCandidate(
            event[3],
            lowmc::pythia8::Pythia8ParticleClassification{
                false,
                true,
                false,
                "toy-final-pion",
            });

    require(candidate.pdgId == 2212, "candidate pdg id");
    require(candidate.isFinalStable, "candidate final-stable flag");
    require(candidate.sourceLabel == "toy-final-pion", "candidate source label");
    require(std::abs(candidate.momentum.E - 1.37) < 1.0e-12, "candidate E");
    require(std::abs(candidate.momentum.px - 1.0) < 1.0e-12, "candidate px");
}

void testCollectAndSelectToyEvent() {
    Pythia8::Event event = makeToyEvent();
    const std::vector<lowmc::AdapterParticleCandidate> candidates =
        lowmc::pythia8::collectAdapterParticleCandidates(event);
    require(candidates.size() == 5, "candidate count excludes unstable intermediate");
    require(candidates[0].isIncomingBeam, "first beam flag");
    require(candidates[1].isIncomingBeam, "second beam flag");
    require(candidates[2].isFinalStable, "first final flag");

    const lowmc::AdapterParticleSelectionResult result =
        lowmc::pythia8::selectChargedParticles(event);
    require(result.selectedParticles.size() == 4, "selected charged count");
    require(result.diagnostic.selectedIncomingCharged == 2, "selected incoming count");
    require(result.diagnostic.selectedOutgoingCharged == 2, "selected outgoing count");
    require(result.diagnostic.ignoredNeutralCandidates == 1, "ignored neutral count");
    require(result.diagnostic.rejectedChargedCount() == 0, "rejected charged count");
    require(std::abs(result.diagnostic.signedChargeSum) < 1.0e-12, "signed charge sum");
    require(std::abs(result.diagnostic.ordinaryChargeSum - 4.0) < 1.0e-12,
            "ordinary charge sum");
}

void testCanDisableIncomingBeams() {
    Pythia8::Event event = makeToyEvent();
    lowmc::ParticleSelectionConfig config;
    config.includeIncomingBeams = false;

    const lowmc::AdapterParticleSelectionResult result =
        lowmc::pythia8::selectChargedParticles(event, config);
    require(result.selectedParticles.size() == 2, "selected charged count without beams");
    require(result.diagnostic.selectedIncomingCharged == 0, "incoming count disabled");
    require(result.diagnostic.selectedOutgoingCharged == 2, "outgoing count without beams");
    require(result.diagnostic.rejectedChargedCount() == 2, "beam rejection count");
}

}  // namespace

int main() {
    testBuildMarker();
    testCandidateConversion();
    testCollectAndSelectToyEvent();
    testCanDisableIncomingBeams();
    return EXIT_SUCCESS;
}
