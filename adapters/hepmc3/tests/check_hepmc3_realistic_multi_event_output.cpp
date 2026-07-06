#include "lowmc/hepmc3/HepMC3Adapter.h"

#include <HepMC3/Attribute.h>
#include <HepMC3/GenEvent.h>
#include <HepMC3/GenParticle.h>
#include <HepMC3/ReaderAscii.h>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory>
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

std::vector<HepMC3::ConstGenParticlePtr> lowmcPhotons(const HepMC3::GenEvent& event) {
    std::vector<HepMC3::ConstGenParticlePtr> photons;
    for (const HepMC3::ConstGenParticlePtr& particle : event.particles()) {
        if (!particle || particle->pid() != 22 || particle->status() != 1) {
            continue;
        }
        if (particle->attribute_as_string("LowMC:origin") == "afterburner") {
            photons.push_back(particle);
        }
    }
    return photons;
}

void requireEnhancedMetadata(
    const HepMC3::GenEvent& event,
    const std::vector<HepMC3::ConstGenParticlePtr>& photons,
    int expectedPhotons,
    double expectedWeight) {
    require(
        event.attribute_as_string("LowMC:photon_mode") == "enhanced",
        "generated event should record enhanced LowMC photon mode");
    requireNear(
        requireAttribute(
            event.attribute<HepMC3::DoubleAttribute>("LowMC:enhancement_factor"),
            "generated event should record LowMC enhancement factor")
            ->value(),
        10.0,
        1.0e-12,
        "event LowMC enhancement factor should match expectation");
    requireNear(
        requireAttribute(
            event.attribute<HepMC3::DoubleAttribute>("LowMC:weight"),
            "generated event should record LowMC weight")
            ->value(),
        expectedWeight,
        1.0e-12,
        "event LowMC weight should match expectation");
    require(
        static_cast<int>(photons.size()) == expectedPhotons,
        "generated event should contain expected LowMC photon count");

    int expectedSampleIndex = 0;
    for (const HepMC3::ConstGenParticlePtr& photon : photons) {
        require(!photon->end_vertex(), "LowMC photons should be attributed root photons");
        require(
            photon->attribute_as_string("LowMC:photon_mode") == "enhanced",
            "LowMC photons should record enhanced photon mode");
        require(
            photon->attribute_as_string("LowMC:shape_status") == "ok",
            "LowMC photons should record successful shape status");
        requireNear(
            requireAttribute(
                photon->attribute<HepMC3::DoubleAttribute>("LowMC:enhancement_factor"),
                "LowMC photons should record enhancement factor")
                ->value(),
            10.0,
            1.0e-12,
            "LowMC photon enhancement factor should match expectation");
        requireNear(
            requireAttribute(
                photon->attribute<HepMC3::DoubleAttribute>("LowMC:weight"),
                "LowMC photons should record weight")
                ->value(),
            expectedWeight,
            1.0e-12,
            "LowMC photon weight should match expectation");
        require(
            requireAttribute(
                photon->attribute<HepMC3::IntAttribute>("LowMC:sample_index"),
                "LowMC photons should record sample index")
                ->value() == expectedSampleIndex,
            "LowMC photon sample index should match order");
        require(
            requireAttribute(
                photon->attribute<HepMC3::DoubleAttribute>("LowMC:max_trial_bound_ratio"),
                "LowMC photons should record max trial/bound ratio")
                ->value() < 1.0,
            "LowMC max trial/bound ratio should be below one");
        ++expectedSampleIndex;
    }
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Usage: check_hepmc3_realistic_multi_event_output <output.hepmc3>\n";
        return EXIT_FAILURE;
    }

    HepMC3::ReaderAscii reader(argv[1]);
    require(!reader.failed(), "could not open output HepMC3 file");

    HepMC3::GenEvent firstEvent;
    HepMC3::GenEvent secondEvent;
    HepMC3::GenEvent thirdEvent;
    require(reader.read_event(firstEvent), "could not read first output event");
    require(reader.read_event(secondEvent), "could not read second output event");
    require(reader.read_event(thirdEvent), "could not read third output event");

    HepMC3::GenEvent extraEvent;
    const bool readExtraEvent = reader.read_event(extraEvent);
    require(
        !readExtraEvent ||
            (extraEvent.particles().empty() &&
             extraEvent.vertices().empty() &&
             extraEvent.beams().empty()),
        "output should contain exactly three non-empty events");

    require(firstEvent.vertices().size() == 1, "first event vertex topology should be preserved");
    require(secondEvent.vertices().size() == 1, "second event vertex topology should be preserved");
    require(thirdEvent.vertices().size() == 1, "third event vertex topology should be preserved");

    const std::vector<HepMC3::ConstGenParticlePtr> firstLowMCPhotons = lowmcPhotons(firstEvent);
    const std::vector<HepMC3::ConstGenParticlePtr> secondLowMCPhotons = lowmcPhotons(secondEvent);
    const std::vector<HepMC3::ConstGenParticlePtr> thirdLowMCPhotons = lowmcPhotons(thirdEvent);
    requireEnhancedMetadata(firstEvent, firstLowMCPhotons, 5, 0.1);
    requireEnhancedMetadata(secondEvent, secondLowMCPhotons, 2, 0.1);
    require(thirdLowMCPhotons.empty(), "skipped event should contain no LowMC photons");

    require(
        thirdEvent.attribute_as_string("LowMC:photon_mode").empty(),
        "skipped event should not record LowMC photon mode");
    require(
        !static_cast<bool>(thirdEvent.attribute<HepMC3::DoubleAttribute>("LowMC:enhancement_factor")),
        "skipped event should not record LowMC enhancement factor");
    require(
        !static_cast<bool>(thirdEvent.attribute<HepMC3::DoubleAttribute>("LowMC:weight")),
        "skipped event should not record LowMC weight");

    const lowmc::AdapterParticleSelectionResult firstSelection =
        lowmc::hepmc3::selectChargedParticles(firstEvent);
    require(
        firstSelection.diagnostic.selectedIncomingCharged == 2,
        "first event should preserve incoming charged particles");
    require(
        firstSelection.diagnostic.selectedOutgoingCharged == 4,
        "first event should preserve original outgoing charged particles");
    require(
        firstSelection.diagnostic.ignoredNeutralCandidates == 7,
        "first event should contain two original neutrals plus five LowMC photons");
    require(
        std::abs(firstSelection.diagnostic.signedChargeSum) < 1.0e-12,
        "first event should remain charge conserving");

    const lowmc::AdapterParticleSelectionResult secondSelection =
        lowmc::hepmc3::selectChargedParticles(secondEvent);
    require(
        secondSelection.diagnostic.selectedIncomingCharged == 2,
        "second event should preserve incoming charged particles");
    require(
        secondSelection.diagnostic.selectedOutgoingCharged == 4,
        "second event should preserve original outgoing charged particles");
    require(
        secondSelection.diagnostic.ignoredNeutralCandidates == 2,
        "second event should contain only added LowMC photons as neutrals");
    require(
        std::abs(secondSelection.diagnostic.signedChargeSum) < 1.0e-12,
        "second event should remain charge conserving");

    const lowmc::AdapterParticleSelectionResult thirdSelection =
        lowmc::hepmc3::selectChargedParticles(thirdEvent);
    require(
        thirdSelection.diagnostic.selectedIncomingCharged == 2,
        "third event should preserve incoming charged particles");
    require(
        thirdSelection.diagnostic.selectedOutgoingCharged == 1,
        "third event should preserve original outgoing charged particles");
    require(
        thirdSelection.diagnostic.ignoredNeutralCandidates == 2,
        "third event should preserve original neutral candidates only");
    requireNear(
        thirdSelection.diagnostic.signedChargeSum,
        -1.0,
        1.0e-12,
        "third event should remain the non-conserved skipped event");

    std::cout << "LowMC HepMC3 realistic multi-event output check passed\n";
    return EXIT_SUCCESS;
}
