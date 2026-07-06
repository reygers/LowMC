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

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "FAILED: " << message << '\n';
        std::exit(EXIT_FAILURE);
    }
}

[[nodiscard]] bool isEmptyTerminalEvent(const HepMC3::GenEvent& event) {
    return event.particles().empty() && event.vertices().empty() && event.beams().empty();
}

[[nodiscard]] bool isLowMCPhoton(const HepMC3::ConstGenParticlePtr& particle) {
    return particle &&
           particle->pid() == 22 &&
           particle->status() == 1 &&
           particle->attribute_as_string("LowMC:origin") == "afterburner";
}

[[nodiscard]] int lowmcPhotonCount(const HepMC3::GenEvent& event) {
    int count = 0;
    for (const HepMC3::ConstGenParticlePtr& particle : event.particles()) {
        if (isLowMCPhoton(particle)) {
            ++count;
            require(!particle->end_vertex(), "LowMC photons should be attributed root photons");
            require(
                particle->attribute_as_string("LowMC:photon_mode") == "enhanced",
                "LowMC photons should record enhanced mode");
            require(
                static_cast<bool>(particle->attribute<HepMC3::DoubleAttribute>("LowMC:weight")),
                "enhanced LowMC photons should record a weight");
        }
    }
    return count;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Usage: check_hepmc3_generated_sample_output <output.hepmc3>\n";
        return EXIT_FAILURE;
    }

    HepMC3::ReaderAscii reader(argv[1]);
    require(!reader.failed(), "could not open output HepMC3 file");

    int eventsRead = 0;
    int generatedEvents = 0;
    int skippedEvents = 0;
    int totalLowMCPhotons = 0;

    while (!reader.failed()) {
        HepMC3::GenEvent event;
        if (!reader.read_event(event) || isEmptyTerminalEvent(event)) {
            break;
        }
        ++eventsRead;
        require(event.vertices().size() == 1, "generated sample should preserve one vertex per event");

        const lowmc::AdapterParticleSelectionResult selection =
            lowmc::hepmc3::selectChargedParticles(event);
        require(
            selection.diagnostic.selectedIncomingCharged == 2,
            "generated sample should preserve two charged beams");
        require(
            selection.diagnostic.selectedOutgoingCharged >= 5,
            "generated sample should contain several stable charged final particles");

        const int photons = lowmcPhotonCount(event);
        totalLowMCPhotons += photons;
        if (photons > 0) {
            ++generatedEvents;
            require(
                event.attribute_as_string("LowMC:photon_mode") == "enhanced",
                "events with LowMC photons should record enhanced event metadata");
            require(
                std::abs(selection.diagnostic.signedChargeSum) < 1.0e-12,
                "events with LowMC photons should be charge conserving");
        } else {
            ++skippedEvents;
            require(
                event.attribute_as_string("LowMC:photon_mode").empty(),
                "events without LowMC photons should not record LowMC event metadata");
        }
    }

    require(eventsRead == 10, "generated sample output should contain ten events");
    require(generatedEvents >= 6, "most generated-sample events should receive LowMC photons");
    require(skippedEvents >= 1, "generated sample should include skipped events");
    require(totalLowMCPhotons >= 20, "generated sample should contain a useful LowMC photon sample");

    std::cout << "LowMC HepMC3 generated realistic sample output check passed"
              << " events=" << eventsRead
              << " generated_events=" << generatedEvents
              << " skipped_events=" << skippedEvents
              << " lowmc_photons=" << totalLowMCPhotons << '\n';
    return EXIT_SUCCESS;
}
