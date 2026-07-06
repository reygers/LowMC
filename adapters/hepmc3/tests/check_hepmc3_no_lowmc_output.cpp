#include <HepMC3/Attribute.h>
#include <HepMC3/GenEvent.h>
#include <HepMC3/GenParticle.h>
#include <HepMC3/ReaderAscii.h>

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
        std::cerr << "Usage: check_hepmc3_no_lowmc_output <output.hepmc3>\n";
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

    require(event.particles().size() == 1, "output should preserve one original neutral particle");
    require(event.vertices().empty(), "output should preserve no-vertex topology");
    require(event.weights().empty(), "output should preserve empty event weight vector");
    require(event.attribute_as_string("LowMC:photon_mode").empty(), "event should have no LowMC mode");
    require(
        !static_cast<bool>(event.attribute<HepMC3::DoubleAttribute>("LowMC:enhancement_factor")),
        "event should have no LowMC enhancement factor");
    require(
        !static_cast<bool>(event.attribute<HepMC3::DoubleAttribute>("LowMC:weight")),
        "event should have no LowMC weight");

    const HepMC3::ConstGenParticlePtr particle = event.particles().front();
    require(static_cast<bool>(particle), "neutral particle should be present");
    require(particle->pid() == 22, "neutral particle should remain a photon");
    require(particle->status() == 1, "neutral particle should remain final-state");
    require(
        particle->attribute_as_string("LowMC:origin").empty(),
        "neutral particle should not gain LowMC origin metadata");

    std::cout << "LowMC HepMC3 no-output check passed\n";
    return EXIT_SUCCESS;
}
