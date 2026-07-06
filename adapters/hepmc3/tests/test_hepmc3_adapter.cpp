#include "lowmc/hepmc3/HepMC3Adapter.h"

#include <HepMC3/FourVector.h>
#include <HepMC3/Attribute.h>
#include <HepMC3/GenEvent.h>
#include <HepMC3/GenParticle.h>
#include <HepMC3/GenVertex.h>
#include <HepMC3/Units.h>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "FAILED: " << message << '\n';
        std::exit(EXIT_FAILURE);
    }
}

void requireNear(double actual, double expected, double tolerance, const std::string& message) {
    if (std::abs(actual - expected) > tolerance) {
        std::cerr << "FAILED: " << message << "\n  actual: " << actual
                  << "\n  expected: " << expected << "\n  tolerance: " << tolerance << '\n';
        std::exit(EXIT_FAILURE);
    }
}

template <class AttributeType>
std::shared_ptr<AttributeType> requireAttribute(
    const std::shared_ptr<AttributeType>& attribute,
    const std::string& message) {
    require(static_cast<bool>(attribute), message);
    return attribute;
}

void testCandidateCopiesPDGAndClassification() {
    const HepMC3::GenParticle particle(
        HepMC3::FourVector{1.0, 2.0, 3.0, 4.0},
        211,
        1);
    const lowmc::hepmc3::HepMC3ParticleClassification classification{
        false,
        true,
        false,
        "status=1,index=7",
    };

    const lowmc::AdapterParticleCandidate candidate =
        lowmc::hepmc3::makeAdapterParticleCandidate(
            particle,
            HepMC3::Units::GEV,
            classification);

    require(candidate.pdgId == 211, "HepMC3 adapter should copy PDG id");
    require(!candidate.isIncomingBeam, "HepMC3 adapter should copy incoming flag");
    require(candidate.isFinalStable, "HepMC3 adapter should copy final-stable flag");
    require(!candidate.isBeamRemnant, "HepMC3 adapter should copy beam-remnant flag");
    require(
        candidate.sourceLabel == "status=1,index=7",
        "HepMC3 adapter should copy source label");
    requireNear(candidate.momentum.E, 4.0, 1.0e-12, "HepMC3 adapter should map energy");
    requireNear(candidate.momentum.px, 1.0, 1.0e-12, "HepMC3 adapter should map px");
    requireNear(candidate.momentum.py, 2.0, 1.0e-12, "HepMC3 adapter should map py");
    requireNear(candidate.momentum.pz, 3.0, 1.0e-12, "HepMC3 adapter should map pz");
}

void testCandidateConvertsMeVToGeV() {
    const HepMC3::GenParticle particle(
        HepMC3::FourVector{1000.0, 2000.0, 3000.0, 4000.0},
        2212,
        4);
    const lowmc::hepmc3::HepMC3ParticleClassification classification{
        true,
        false,
        false,
        "beam",
    };

    const lowmc::AdapterParticleCandidate candidate =
        lowmc::hepmc3::makeAdapterParticleCandidate(
            particle,
            HepMC3::Units::MEV,
            classification);

    require(candidate.isIncomingBeam, "HepMC3 adapter should copy incoming beam flag");
    requireNear(candidate.momentum.E, 4.0, 1.0e-12, "HepMC3 adapter should convert E to GeV");
    requireNear(candidate.momentum.px, 1.0, 1.0e-12, "HepMC3 adapter should convert px to GeV");
    requireNear(candidate.momentum.py, 2.0, 1.0e-12, "HepMC3 adapter should convert py to GeV");
    requireNear(candidate.momentum.pz, 3.0, 1.0e-12, "HepMC3 adapter should convert pz to GeV");
}

void testCollectCandidatesFromHandBuiltEvent() {
    HepMC3::GenEvent event(HepMC3::Units::MEV, HepMC3::Units::MM);

    auto beamPlus = std::make_shared<HepMC3::GenParticle>(
        HepMC3::FourVector{0.0, 0.0, 7000.0, 7000.0},
        2212,
        4);
    auto beamMinus = std::make_shared<HepMC3::GenParticle>(
        HepMC3::FourVector{0.0, 0.0, -7000.0, 7000.0},
        2212,
        4);
    event.set_beam_particles(beamPlus, beamMinus);

    auto finalPion = std::make_shared<HepMC3::GenParticle>(
        HepMC3::FourVector{1000.0, 2000.0, 3000.0, 4000.0},
        211,
        1);
    auto finalPhoton = std::make_shared<HepMC3::GenParticle>(
        HepMC3::FourVector{500.0, 0.0, 0.0, 500.0},
        22,
        1);
    auto decayingPion = std::make_shared<HepMC3::GenParticle>(
        HepMC3::FourVector{0.0, 0.0, 1000.0, 1200.0},
        211,
        1);

    auto hardVertex = std::make_shared<HepMC3::GenVertex>();
    hardVertex->add_particle_in(beamPlus);
    hardVertex->add_particle_in(beamMinus);
    hardVertex->add_particle_out(finalPion);
    hardVertex->add_particle_out(finalPhoton);
    hardVertex->add_particle_out(decayingPion);
    event.add_vertex(hardVertex);

    auto decayProduct = std::make_shared<HepMC3::GenParticle>(
        HepMC3::FourVector{0.0, 0.0, 500.0, 600.0},
        111,
        1);
    auto decayVertex = std::make_shared<HepMC3::GenVertex>();
    decayVertex->add_particle_in(decayingPion);
    event.add_vertex(decayVertex);
    decayVertex->add_particle_out(decayProduct);

    require(
        static_cast<bool>(decayingPion->end_vertex()),
        "test setup should give decaying pion an end vertex");

    const std::vector<lowmc::AdapterParticleCandidate> candidates =
        lowmc::hepmc3::collectAdapterParticleCandidates(event);

    require(
        candidates.size() == 5,
        "HepMC3 event collector should return beams and finals, actual=" +
            std::to_string(candidates.size()));
    require(candidates[0].isIncomingBeam, "first beam should be incoming");
    require(candidates[1].isIncomingBeam, "second beam should be incoming");
    require(!candidates[0].isFinalStable, "beam should not be final stable");
    require(candidates[2].isFinalStable, "final pion should be final stable");
    require(candidates[3].isFinalStable, "final photon should be final stable");
    require(candidates[4].isFinalStable, "decay product should be final stable");
    require(candidates[0].pdgId == 2212, "first beam PDG id");
    require(candidates[1].pdgId == 2212, "second beam PDG id");
    require(candidates[2].pdgId == 211, "final pion PDG id");
    require(candidates[3].pdgId == 22, "final photon PDG id");
    require(candidates[4].pdgId == 111, "decay product PDG id");
    requireNear(candidates[0].momentum.E, 7.0, 1.0e-12, "beam energy converted to GeV");
    requireNear(candidates[1].momentum.pz, -7.0, 1.0e-12, "beam pz converted to GeV");
    requireNear(candidates[2].momentum.E, 4.0, 1.0e-12, "final pion energy converted to GeV");
    requireNear(candidates[2].momentum.px, 1.0, 1.0e-12, "final pion px converted to GeV");
    require(
        candidates[2].sourceLabel.find("pid=211") != std::string::npos,
        "source label should contain PDG id");
    require(
        candidates[2].sourceLabel.find("status=1") != std::string::npos,
        "source label should contain status");

    const lowmc::AdapterParticleSelectionResult selected =
        lowmc::hepmc3::selectChargedParticles(event);
    require(
        selected.selectedParticles.size() == 3,
        "selection should keep charged beams and final pion");
    require(
        selected.diagnostic.ignoredNeutralCandidates == 2,
        "selection should count neutral final candidates");
    require(
        selected.diagnostic.selectedIncomingCharged == 2,
        "selection should count selected incoming beams");
    require(
        selected.diagnostic.selectedOutgoingCharged == 1,
        "selection should count selected final pion");

    lowmc::ParticleSelectionConfig config;
    config.includeIncomingBeams = false;
    const lowmc::AdapterParticleSelectionResult withoutIncoming =
        lowmc::hepmc3::selectChargedParticles(event, config);
    require(
        withoutIncoming.selectedParticles.size() == 1,
        "selection config should allow disabling incoming beams");
    require(
        withoutIncoming.diagnostic.selectedIncomingCharged == 0,
        "disabled incoming beams should not be selected");
    require(
        withoutIncoming.diagnostic.selectedOutgoingCharged == 1,
        "disabled incoming beams should leave final pion selected");
    require(
        withoutIncoming.diagnostic.ignoredNeutralCandidates == 2,
        "disabled incoming beams should preserve neutral counts");
    require(
        withoutIncoming.diagnostic.rejectedChargedCount() == 2,
        "disabled incoming beams should be reported as rejected charged candidates");
    require(
        withoutIncoming.diagnostic.rejectedChargedCandidates[0].reason ==
            lowmc::RejectedParticleReason::OutsideSelection,
        "disabled incoming beams should use outside-selection rejection reason");
    require(
        withoutIncoming.diagnostic.rejectedChargedCandidates[0].sourceLabel.find("pid=2212") !=
            std::string::npos,
        "rejected incoming beam should retain HepMC3 source label");
}

void testAppendUnweightedLowPhotons() {
    HepMC3::GenEvent event(HepMC3::Units::MEV, HepMC3::Units::MM);
    event.weights().push_back(2.5);

    auto beamPlus = std::make_shared<HepMC3::GenParticle>(
        HepMC3::FourVector{0.0, 0.0, 7000.0, 7000.0},
        2212,
        4);
    auto beamMinus = std::make_shared<HepMC3::GenParticle>(
        HepMC3::FourVector{0.0, 0.0, -7000.0, 7000.0},
        2212,
        4);
    event.set_beam_particles(beamPlus, beamMinus);

    auto finalPion = std::make_shared<HepMC3::GenParticle>(
        HepMC3::FourVector{1000.0, 0.0, 0.0, 1100.0},
        211,
        1);
    auto hardVertex = std::make_shared<HepMC3::GenVertex>();
    hardVertex->add_particle_in(beamPlus);
    hardVertex->add_particle_in(beamMinus);
    hardVertex->add_particle_out(finalPion);
    event.add_vertex(hardVertex);

    const std::size_t particlesBefore = event.particles().size();
    const std::size_t verticesBefore = event.vertices().size();
    const HepMC3::FourVector finalPionMomentumBefore = finalPion->momentum();

    const std::vector<lowmc::ShapeUnweightedPhotonCandidate> photons{
        {
            lowmc::Photon{lowmc::FourVector{0.5, 0.3, 0.4, 0.0}, 1.0},
            0.5,
            0.0,
            0.9272952180016123,
        },
        {
            lowmc::Photon{lowmc::FourVector{1.0, 0.0, 0.6, 0.8}, 1.0},
            0.6,
            1.0986122886681098,
            1.5707963267948966,
        },
    };
    const lowmc::hepmc3::HepMC3LowPhotonOutputOptions options{
        1,
        "test-version",
        0.125,
        true,
        "ordinary",
        "ok",
        0.25,
        true,
    };

    const lowmc::hepmc3::HepMC3LowPhotonOutputResult result =
        lowmc::hepmc3::appendUnweightedLowPhotons(event, photons, options);

    require(result.photonsAdded == 2, "writer should report appended photon count");
    require(event.particles().size() == particlesBefore + 2, "writer should append particles");
    require(event.vertices().size() == verticesBefore, "writer should preserve vertex count");
    require(event.weights().size() == 1, "writer should preserve event weight count");
    requireNear(event.weights()[0], 2.5, 1.0e-12, "writer should preserve event weight value");
    require(
        event.attribute_as_string("LowMC:photon_mode") == "unweighted",
        "writer should record event photon mode");
    requireNear(
        requireAttribute(
            event.attribute<HepMC3::DoubleAttribute>("LowMC:enhancement_factor"),
            "writer should record event enhancement factor")
            ->value(),
        1.0,
        1.0e-12,
        "event enhancement factor value");
    require(
        !static_cast<bool>(event.attribute<HepMC3::DoubleAttribute>("LowMC:weight")),
        "unweighted physical event should not require an event weight attribute");
    requireNear(finalPion->momentum().px(), finalPionMomentumBefore.px(), 1.0e-12, "writer should not recoil original px");
    requireNear(finalPion->momentum().e(), finalPionMomentumBefore.e(), 1.0e-12, "writer should not recoil original energy");

    HepMC3::GenParticlePtr firstPhoton;
    for (const HepMC3::GenParticlePtr& particle : event.particles()) {
        if (particle->attribute_as_string("LowMC:origin") == "afterburner") {
            firstPhoton = particle;
            break;
        }
    }
    require(static_cast<bool>(firstPhoton), "writer should attach a LowMC photon");
    require(firstPhoton->pid() == 22, "LowMC photon should have PDG id 22");
    require(firstPhoton->status() == 1, "LowMC photon should have status 1");
    require(!firstPhoton->end_vertex(), "LowMC photon should be final-state");
    require(!firstPhoton->is_generated_mass_set(), "LowMC photon generated mass should be unset");
    requireNear(firstPhoton->momentum().px(), 300.0, 1.0e-12, "LowMC photon px should be converted to MeV");
    requireNear(firstPhoton->momentum().py(), 400.0, 1.0e-12, "LowMC photon py should be converted to MeV");
    requireNear(firstPhoton->momentum().pz(), 0.0, 1.0e-12, "LowMC photon pz should be converted to MeV");
    requireNear(firstPhoton->momentum().e(), 500.0, 1.0e-12, "LowMC photon E should be converted to MeV");
    require(
        firstPhoton->attribute_as_string("LowMC:origin") == "afterburner",
        "LowMC photon should record origin");
    require(
        firstPhoton->attribute_as_string("LowMC:model") == "leading-low",
        "LowMC photon should record model");
    require(
        firstPhoton->attribute_as_string("LowMC:recoil") == "none",
        "LowMC photon should record no recoil");
    require(
        firstPhoton->attribute_as_string("LowMC:photon_mode") == "unweighted",
        "LowMC photon should record photon mode");
    requireNear(
        requireAttribute(
            firstPhoton->attribute<HepMC3::DoubleAttribute>("LowMC:enhancement_factor"),
            "LowMC photon should record enhancement factor")
            ->value(),
        1.0,
        1.0e-12,
        "LowMC photon enhancement factor value");
    require(
        !static_cast<bool>(firstPhoton->attribute<HepMC3::DoubleAttribute>("LowMC:weight")),
        "unweighted physical LowMC photon should not require a weight attribute");
    require(
        firstPhoton->attribute_as_string("LowMC:version") == "test-version",
        "LowMC photon should record version when provided");
    require(
        firstPhoton->attribute_as_string("LowMC:validation_class") == "ordinary",
        "LowMC photon should record validation class");
    require(
        firstPhoton->attribute_as_string("LowMC:shape_status") == "ok",
        "LowMC photon should record shape status");
    requireNear(
        requireAttribute(
            firstPhoton->attribute<HepMC3::DoubleAttribute>("LowMC:lambda_event"),
            "LowMC photon should record lambda_event")
            ->value(),
        0.125,
        1.0e-12,
        "LowMC photon lambda_event value");
    requireNear(
        requireAttribute(
            firstPhoton->attribute<HepMC3::DoubleAttribute>("LowMC:max_trial_bound_ratio"),
            "LowMC photon should record max_trial_bound_ratio")
            ->value(),
        0.25,
        1.0e-12,
        "LowMC photon max_trial_bound_ratio value");
    require(
        requireAttribute(
            firstPhoton->attribute<HepMC3::IntAttribute>("LowMC:sample_index"),
            "LowMC photon should record sample index")
            ->value() == 0,
        "LowMC first photon sample index should be zero");
    requireNear(
        requireAttribute(
            firstPhoton->attribute<HepMC3::DoubleAttribute>("LowMC:kT"),
            "LowMC photon should record kT")
            ->value(),
        0.5,
        1.0e-12,
        "LowMC first photon kT attribute");
    requireNear(
        requireAttribute(
            firstPhoton->attribute<HepMC3::DoubleAttribute>("LowMC:eta"),
            "LowMC photon should record eta")
            ->value(),
        0.0,
        1.0e-12,
        "LowMC first photon eta attribute");
    requireNear(
        requireAttribute(
            firstPhoton->attribute<HepMC3::DoubleAttribute>("LowMC:phi"),
            "LowMC photon should record phi")
            ->value(),
        0.9272952180016123,
        1.0e-12,
        "LowMC first photon phi attribute");
}

void testAppendEnhancedLowPhotons() {
    HepMC3::GenEvent event(HepMC3::Units::GEV, HepMC3::Units::MM);
    const std::vector<lowmc::ShapeUnweightedPhotonCandidate> photons{
        {
            lowmc::Photon{lowmc::FourVector{1.0, 1.0, 0.0, 0.0}, 1.0},
            1.0,
            0.0,
            0.0,
        },
    };
    lowmc::hepmc3::HepMC3LowPhotonOutputOptions options;
    options.enhancementFactor = 5.0;

    const lowmc::hepmc3::HepMC3LowPhotonOutputResult result =
        lowmc::hepmc3::appendUnweightedLowPhotons(event, photons, options);

    require(result.photonsAdded == 1, "enhanced writer should report appended photon count");
    require(event.particles().size() == 1, "enhanced writer should append one particle");
    require(
        event.attribute_as_string("LowMC:photon_mode") == "enhanced",
        "enhanced writer should record event photon mode");
    requireNear(
        requireAttribute(
            event.attribute<HepMC3::DoubleAttribute>("LowMC:enhancement_factor"),
            "enhanced writer should record event enhancement factor")
            ->value(),
        5.0,
        1.0e-12,
        "enhanced event enhancement factor value");
    requireNear(
        requireAttribute(
            event.attribute<HepMC3::DoubleAttribute>("LowMC:weight"),
            "enhanced writer should record event analysis weight")
            ->value(),
        0.2,
        1.0e-12,
        "enhanced event weight should be inverse enhancement factor");
    const HepMC3::GenParticlePtr photon = event.particles().front();
    require(
        photon->attribute_as_string("LowMC:photon_mode") == "enhanced",
        "enhanced LowMC photon should record enhanced photon mode");
    requireNear(
        requireAttribute(
            photon->attribute<HepMC3::DoubleAttribute>("LowMC:enhancement_factor"),
            "enhanced LowMC photon should record enhancement factor")
            ->value(),
        5.0,
        1.0e-12,
        "enhanced LowMC photon enhancement factor value");
    requireNear(
        requireAttribute(
            photon->attribute<HepMC3::DoubleAttribute>("LowMC:weight"),
            "enhanced LowMC photon should record analysis weight")
            ->value(),
        0.2,
        1.0e-12,
        "enhanced LowMC photon weight should be inverse enhancement factor");
}

void testAppendRejectsInvalidEnhancementFactor() {
    HepMC3::GenEvent event(HepMC3::Units::GEV, HepMC3::Units::MM);
    const std::vector<lowmc::ShapeUnweightedPhotonCandidate> photons{
        {
            lowmc::Photon{lowmc::FourVector{1.0, 1.0, 0.0, 0.0}, 1.0},
            1.0,
            0.0,
            0.0,
        },
    };
    lowmc::hepmc3::HepMC3LowPhotonOutputOptions options;
    options.enhancementFactor = 0.0;

    bool threw = false;
    try {
        (void)lowmc::hepmc3::appendUnweightedLowPhotons(event, photons, options);
    } catch (const std::invalid_argument&) {
        threw = true;
    }

    require(threw, "HepMC3 writer should reject invalid enhancement factors");
    require(event.particles().empty(), "enhancement rejection should not modify particles");
    require(event.vertices().empty(), "enhancement rejection should not modify vertices");
}

}  // namespace

int main() {
    testCandidateCopiesPDGAndClassification();
    testCandidateConvertsMeVToGeV();
    testCollectCandidatesFromHandBuiltEvent();
    testAppendUnweightedLowPhotons();
    testAppendEnhancedLowPhotons();
    testAppendRejectsInvalidEnhancementFactor();

    std::cout << "All LowMC HepMC3 adapter tests passed\n";
    return EXIT_SUCCESS;
}
