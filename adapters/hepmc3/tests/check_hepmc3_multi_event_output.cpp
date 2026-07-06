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

}  // namespace

int main(int argc, char** argv) {
    if (argc != 6) {
        std::cerr << "Usage: check_hepmc3_multi_event_output <output.hepmc3> "
                     "<photon_mode> <enhancement_factor> <weight> <first_event_low_photons>\n";
        return EXIT_FAILURE;
    }

    const std::string expectedMode = argv[2];
    const double expectedEnhancementFactor = std::stod(argv[3]);
    const double expectedWeight = std::stod(argv[4]);
    const int expectedFirstEventLowPhotons = std::stoi(argv[5]);

    HepMC3::ReaderAscii reader(argv[1]);
    require(!reader.failed(), "could not open output HepMC3 file");

    HepMC3::GenEvent firstEvent;
    HepMC3::GenEvent secondEvent;
    require(reader.read_event(firstEvent), "could not read first output event");
    require(reader.read_event(secondEvent), "could not read second output event");

    HepMC3::GenEvent extraEvent;
    const bool readExtraEvent = reader.read_event(extraEvent);
    require(
        !readExtraEvent ||
            (extraEvent.particles().empty() &&
             extraEvent.vertices().empty() &&
             extraEvent.beams().empty()),
        "output should contain exactly two non-empty events");

    require(firstEvent.weights().empty(), "first event weight vector should remain empty");
    require(secondEvent.weights().empty(), "second event weight vector should remain empty");
    require(firstEvent.vertices().size() == 1, "first event vertex topology should be preserved");
    require(secondEvent.vertices().size() == 2, "second event vertex topology should be preserved");
    require(
        firstEvent.attribute_as_string("LowMC:photon_mode") == expectedMode,
        "first event should record expected LowMC photon mode");
    requireNear(
        requireAttribute(
            firstEvent.attribute<HepMC3::DoubleAttribute>("LowMC:enhancement_factor"),
            "first event should record LowMC enhancement factor")
            ->value(),
        expectedEnhancementFactor,
        1.0e-12,
        "first event LowMC enhancement factor should match expectation");
    if (expectedMode == "enhanced") {
        requireNear(
            requireAttribute(
                firstEvent.attribute<HepMC3::DoubleAttribute>("LowMC:weight"),
                "enhanced first event should record LowMC weight")
                ->value(),
            expectedWeight,
            1.0e-12,
            "enhanced first event LowMC weight should match expectation");
    } else {
        require(
            !static_cast<bool>(firstEvent.attribute<HepMC3::DoubleAttribute>("LowMC:weight")),
            "unweighted first event should not require a LowMC weight attribute");
    }
    require(
        secondEvent.attribute_as_string("LowMC:photon_mode").empty(),
        "skipped second event should not record LowMC photon mode");
    require(
        !static_cast<bool>(secondEvent.attribute<HepMC3::DoubleAttribute>("LowMC:enhancement_factor")),
        "skipped second event should not record LowMC enhancement factor");
    require(
        !static_cast<bool>(secondEvent.attribute<HepMC3::DoubleAttribute>("LowMC:weight")),
        "skipped second event should not record LowMC weight");

    const std::vector<HepMC3::ConstGenParticlePtr> firstLowMCPhotons = lowmcPhotons(firstEvent);
    const std::vector<HepMC3::ConstGenParticlePtr> secondLowMCPhotons = lowmcPhotons(secondEvent);
    require(
        static_cast<int>(firstLowMCPhotons.size()) == expectedFirstEventLowPhotons,
        "first event should contain expected LowMC photon count");
    require(secondLowMCPhotons.empty(), "skipped second event should contain no LowMC photons");

    const lowmc::AdapterParticleSelectionResult firstSelection =
        lowmc::hepmc3::selectChargedParticles(firstEvent);
    require(
        firstSelection.diagnostic.selectedIncomingCharged == 2,
        "first event should preserve incoming charged particles");
    require(
        firstSelection.diagnostic.selectedOutgoingCharged == 2,
        "first event should preserve original outgoing charged particles");
    require(
        firstSelection.diagnostic.ignoredNeutralCandidates ==
            static_cast<std::size_t>(1 + expectedFirstEventLowPhotons),
        "first event should contain original neutral plus LowMC photons");
    require(
        std::abs(firstSelection.diagnostic.signedChargeSum) < 1.0e-12,
        "first event should remain charge conserving");

    const lowmc::AdapterParticleSelectionResult secondSelection =
        lowmc::hepmc3::selectChargedParticles(secondEvent);
    require(
        secondSelection.diagnostic.selectedIncomingCharged == 2,
        "second event should preserve incoming charged particles");
    require(
        secondSelection.diagnostic.selectedOutgoingCharged == 1,
        "second event should preserve original outgoing charged particles");
    require(
        secondSelection.diagnostic.ignoredNeutralCandidates == 2,
        "second event should preserve original neutral candidates only");
    requireNear(
        secondSelection.diagnostic.signedChargeSum,
        -1.0,
        1.0e-12,
        "second event should remain the non-conserved skipped fixture");

    for (const HepMC3::ConstGenParticlePtr& photon : firstLowMCPhotons) {
        require(!photon->end_vertex(), "LowMC photons should be attributed root photons");
        require(
            photon->attribute_as_string("LowMC:origin") == "afterburner",
            "LowMC photons should record origin");
        require(
            photon->attribute_as_string("LowMC:photon_mode") == expectedMode,
            "LowMC photons should record expected photon mode");
        requireNear(
            requireAttribute(
                photon->attribute<HepMC3::DoubleAttribute>("LowMC:enhancement_factor"),
                "LowMC photons should record enhancement factor")
                ->value(),
            expectedEnhancementFactor,
            1.0e-12,
            "LowMC enhancement factor should match expectation");
        if (expectedMode == "enhanced") {
            requireNear(
                requireAttribute(
                    photon->attribute<HepMC3::DoubleAttribute>("LowMC:weight"),
                    "enhanced LowMC photons should record weight")
                    ->value(),
                expectedWeight,
                1.0e-12,
                "enhanced LowMC photon weight should match expectation");
        } else {
            require(
                !static_cast<bool>(photon->attribute<HepMC3::DoubleAttribute>("LowMC:weight")),
                "unweighted LowMC photons should not require a weight attribute");
        }
    }

    std::cout << "LowMC HepMC3 multi-event output check passed\n";
    return EXIT_SUCCESS;
}
