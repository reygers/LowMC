#include "lowmc/hepmc3/HepMC3Adapter.h"

#include <HepMC3/Attribute.h>
#include <HepMC3/GenEvent.h>
#include <HepMC3/GenParticle.h>
#include <HepMC3/ReaderAscii.h>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "FAILED: " << message << '\n';
        std::exit(EXIT_FAILURE);
    }
}

bool isEmptyTerminalEvent(const HepMC3::GenEvent& event) {
    return event.particles().empty() && event.vertices().empty() && event.beams().empty();
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Usage: check_hepmc3_charge_conserving_no_lowmc_output <output.hepmc3>\n";
        return EXIT_FAILURE;
    }

    HepMC3::ReaderAscii reader(argv[1]);
    require(!reader.failed(), "could not open output HepMC3 file");

    HepMC3::GenEvent event;
    require(reader.read_event(event), "could not read output HepMC3 event");
    HepMC3::GenEvent extraEvent;
    const bool readExtraEvent = reader.read_event(extraEvent);
    require(
        !readExtraEvent || isEmptyTerminalEvent(extraEvent),
        "output should contain exactly one non-empty event");

    require(event.particles().size() == 5, "output should preserve the five input particles");
    require(event.vertices().size() == 1, "output should preserve the original vertex topology");
    require(event.weights().empty(), "output should preserve empty event weight vector");
    require(event.attribute_as_string("LowMC:photon_mode").empty(), "event should have no LowMC mode");
    require(
        !static_cast<bool>(event.attribute<HepMC3::DoubleAttribute>("LowMC:enhancement_factor")),
        "event should have no LowMC enhancement factor");
    require(
        !static_cast<bool>(event.attribute<HepMC3::DoubleAttribute>("LowMC:weight")),
        "event should have no LowMC weight");

    int lowmcPhotonCount = 0;
    for (const HepMC3::GenParticlePtr& particle : event.particles()) {
        require(static_cast<bool>(particle), "all particles should be readable");
        if (particle->attribute_as_string("LowMC:origin") == "afterburner") {
            ++lowmcPhotonCount;
        }
    }
    require(lowmcPhotonCount == 0, "output should contain no LowMC photons");

    const lowmc::AdapterParticleSelectionResult selection =
        lowmc::hepmc3::selectChargedParticles(event);
    require(
        selection.diagnostic.selectedIncomingCharged == 2,
        "output should preserve two incoming charged particles");
    require(
        selection.diagnostic.selectedOutgoingCharged == 2,
        "output should preserve two outgoing charged particles");
    require(
        selection.diagnostic.ignoredNeutralCandidates == 1,
        "output should preserve only the original neutral candidate");
    require(
        selection.diagnostic.rejectedChargedCount() == 0,
        "output should not add rejected charged candidates");
    require(
        std::abs(selection.diagnostic.signedChargeSum) < 1.0e-12,
        "output should remain charge conserving");

    std::cout << "LowMC HepMC3 charge-conserving no-output check passed\n";
    return EXIT_SUCCESS;
}
