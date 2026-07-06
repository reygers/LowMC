#include "lowmc/hepmc3/HepMC3Adapter.h"

#include <HepMC3/Attribute.h>
#include <HepMC3/GenEvent.h>
#include <HepMC3/GenParticle.h>
#include <HepMC3/GenVertex.h>
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

}  // namespace

int main(int argc, char** argv) {
    if (argc != 6 && argc != 7) {
        std::cerr << "Usage: check_hepmc3_low_photon_output <output.hepmc3> "
                     "<photon_mode> <enhancement_factor> <weight> <low_photons> "
                     "[lambda_event]\n";
        return EXIT_FAILURE;
    }
    const std::string expectedMode = argv[2];
    const double expectedEnhancementFactor = std::stod(argv[3]);
    const double expectedWeight = std::stod(argv[4]);
    const int expectedLowPhotons = std::stoi(argv[5]);
    const double expectedLambdaEvent = argc == 7 ? std::stod(argv[6]) : 0.144148;

    HepMC3::ReaderAscii reader(argv[1]);
    require(!reader.failed(), "could not open output HepMC3 file");

    HepMC3::GenEvent event;
    require(reader.read_event(event), "could not read output HepMC3 event");

    const lowmc::AdapterParticleSelectionResult selection =
        lowmc::hepmc3::selectChargedParticles(event);
    require(
        selection.diagnostic.selectedIncomingCharged == 2,
        "output should preserve two incoming charged particles");
    require(
        selection.diagnostic.selectedOutgoingCharged == 2,
        "output should preserve two original outgoing charged particles");
    require(
        selection.diagnostic.ignoredNeutralCandidates ==
            static_cast<std::size_t>(1 + expectedLowPhotons),
        "output should contain original neutral plus expected LowMC photons");
    require(
        selection.diagnostic.rejectedChargedCount() == 0,
        "output should not add rejected charged candidates");
    require(
        event.attribute_as_string("LowMC:photon_mode") == expectedMode,
        "event should record expected LowMC photon mode");
    requireNear(
        requireAttribute(
            event.attribute<HepMC3::DoubleAttribute>("LowMC:enhancement_factor"),
            "event should record LowMC enhancement_factor")
            ->value(),
        expectedEnhancementFactor,
        1.0e-12,
        "event LowMC enhancement_factor should match expectation");
    if (expectedMode == "enhanced") {
        requireNear(
            requireAttribute(
                event.attribute<HepMC3::DoubleAttribute>("LowMC:weight"),
                "enhanced event should record LowMC weight")
                ->value(),
            expectedWeight,
            1.0e-12,
            "enhanced event LowMC weight should match expectation");
    } else {
        require(
            !static_cast<bool>(event.attribute<HepMC3::DoubleAttribute>("LowMC:weight")),
            "unweighted event should not require a LowMC weight attribute");
    }

    std::vector<HepMC3::GenParticlePtr> lowmcPhotons;
    for (const HepMC3::GenParticlePtr& particle : event.particles()) {
        if (!particle || particle->pid() != 22 || particle->status() != 1) {
            continue;
        }
        if (particle->attribute_as_string("LowMC:origin") == "afterburner") {
            lowmcPhotons.push_back(particle);
        }
    }

    require(
        static_cast<int>(lowmcPhotons.size()) == expectedLowPhotons,
        "output should contain expected number of LowMC photons");
    HepMC3::GenParticlePtr lowmcPhoton;
    if (!lowmcPhotons.empty()) {
        lowmcPhoton = lowmcPhotons.front();
    }
    require(static_cast<bool>(lowmcPhoton), "output should contain a LowMC photon");
    for (const HepMC3::GenParticlePtr& photon : lowmcPhotons) {
        require(!photon->end_vertex(), "LowMC photon should be final-state");
        require(
            photon->attribute_as_string("LowMC:photon_mode") == expectedMode,
            "LowMC photon should record expected photon mode");
        requireNear(
            requireAttribute(
                photon->attribute<HepMC3::DoubleAttribute>("LowMC:enhancement_factor"),
                "LowMC photon should record enhancement_factor")
                ->value(),
            expectedEnhancementFactor,
            1.0e-12,
            "LowMC enhancement_factor should match expectation");
        if (expectedMode == "enhanced") {
            requireNear(
                requireAttribute(
                    photon->attribute<HepMC3::DoubleAttribute>("LowMC:weight"),
                    "enhanced LowMC photon should record weight")
                    ->value(),
                expectedWeight,
                1.0e-12,
                "enhanced LowMC photon weight should match expectation");
        } else {
            require(
                !static_cast<bool>(photon->attribute<HepMC3::DoubleAttribute>("LowMC:weight")),
                "unweighted LowMC photon should not require a weight attribute");
        }
    }

    require(
        lowmcPhoton->attribute_as_string("LowMC:origin") == "afterburner",
        "LowMC photon should record origin");
    require(
        lowmcPhoton->attribute_as_string("LowMC:model") == "leading-low",
        "LowMC photon should record model");
    require(
        lowmcPhoton->attribute_as_string("LowMC:recoil") == "none",
        "LowMC photon should record no recoil");
    require(
        lowmcPhoton->attribute_as_string("LowMC:photon_mode") == expectedMode,
        "LowMC photon should record expected photon mode");
    require(
        lowmcPhoton->attribute_as_string("LowMC:validation_class") == "ordinary",
        "LowMC photon should record validation class");
    require(
        lowmcPhoton->attribute_as_string("LowMC:shape_status") == "ok",
        "LowMC photon should record successful shape status");

    requireNear(
        requireAttribute(
            lowmcPhoton->attribute<HepMC3::DoubleAttribute>("LowMC:lambda_event"),
            "LowMC photon should record lambda_event")
            ->value(),
        expectedLambdaEvent,
        1.0e-5,
        "LowMC lambda_event should match fixture yield");
    require(
        requireAttribute(
            lowmcPhoton->attribute<HepMC3::DoubleAttribute>("LowMC:max_trial_bound_ratio"),
            "LowMC photon should record max trial/bound ratio")
            ->value() < 1.0,
        "LowMC max trial/bound ratio should be below one");

    require(
        requireAttribute(
            lowmcPhoton->attribute<HepMC3::IntAttribute>("LowMC:sample_index"),
            "LowMC photon should record sample index")
            ->value() == 0,
        "LowMC photon sample index should be zero");
    require(
        static_cast<bool>(lowmcPhoton->attribute<HepMC3::DoubleAttribute>("LowMC:kT")),
        "LowMC photon should record kT");
    require(
        static_cast<bool>(lowmcPhoton->attribute<HepMC3::DoubleAttribute>("LowMC:eta")),
        "LowMC photon should record eta");
    require(
        static_cast<bool>(lowmcPhoton->attribute<HepMC3::DoubleAttribute>("LowMC:phi")),
        "LowMC photon should record phi");

    std::cout << "LowMC HepMC3 output file check passed\n";
    return EXIT_SUCCESS;
}
