#include <HepMC3/FourVector.h>
#include <HepMC3/GenEvent.h>
#include <HepMC3/GenParticle.h>
#include <HepMC3/GenVertex.h>
#include <HepMC3/Units.h>
#include <HepMC3/WriterAscii.h>

#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct ParticleSpec {
    int pid;
    double px;
    double py;
    double pz;
    double mass;
};

[[nodiscard]] double energy(double px, double py, double pz, double mass) {
    return std::sqrt(px * px + py * py + pz * pz + mass * mass);
}

[[nodiscard]] HepMC3::GenParticlePtr makeParticle(const ParticleSpec& spec, int status) {
    return std::make_shared<HepMC3::GenParticle>(
        HepMC3::FourVector(
            spec.px,
            spec.py,
            spec.pz,
            energy(spec.px, spec.py, spec.pz, spec.mass)),
        spec.pid,
        status);
}

void addEvent(HepMC3::WriterAscii& writer, int eventNumber, const std::vector<ParticleSpec>& finalState) {
    HepMC3::GenEvent event(HepMC3::Units::MEV, HepMC3::Units::MM);
    event.set_event_number(eventNumber);

    const ParticleSpec forwardBeam{2212, 0.0, 0.0, 7.0e6, 938.2720813};
    const ParticleSpec backwardBeam{2212, 0.0, 0.0, -7.0e6, 938.2720813};
    HepMC3::GenParticlePtr beam1 = makeParticle(forwardBeam, 4);
    HepMC3::GenParticlePtr beam2 = makeParticle(backwardBeam, 4);
    auto vertex = std::make_shared<HepMC3::GenVertex>();
    vertex->add_particle_in(beam1);
    vertex->add_particle_in(beam2);

    for (const ParticleSpec& particle : finalState) {
        vertex->add_particle_out(makeParticle(particle, 1));
    }

    event.add_vertex(vertex);
    event.set_beam_particles(beam1, beam2);
    writer.write_event(event);
}

std::vector<ParticleSpec> shiftedEvent(int eventIndex) {
    const double scale = 1.0 + 0.06 * static_cast<double>(eventIndex);
    const double twist = 35.0 * static_cast<double>(eventIndex % 4);
    std::vector<ParticleSpec> particles{
        {211, 700.0 * scale + twist, 120.0, 900.0 + 80.0 * eventIndex, 139.57039},
        {-211, -520.0 * scale, -260.0 - twist, -650.0 - 50.0 * eventIndex, 139.57039},
        {321, -420.0, 310.0 * scale, 500.0 + 40.0 * eventIndex, 493.677},
        {-321, 300.0 + twist, -220.0 * scale, -420.0, 493.677},
        {2212, 380.0 * scale, 150.0, 1200.0 + 90.0 * eventIndex, 938.2720813},
        {2212, -340.0, -130.0 * scale, -950.0 - 60.0 * eventIndex, 938.2720813},
        {22, 160.0, -90.0, 250.0 + 20.0 * eventIndex, 0.0},
        {111, -180.0, 60.0, -300.0, 134.9768},
    };

    if (eventIndex % 3 == 1) {
        particles.push_back({211, 240.0, -150.0, 420.0, 139.57039});
        particles.push_back({-211, -210.0, 180.0, -380.0, 139.57039});
    }

    if (eventIndex == 4 || eventIndex == 9) {
        particles.pop_back();
    }
    return particles;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Usage: lowmc_hepmc3_generate_realistic_sample <output.hepmc3>\n";
        return EXIT_FAILURE;
    }

    try {
        HepMC3::WriterAscii writer(argv[1]);
        if (writer.failed()) {
            throw std::runtime_error("could not open output file");
        }

        constexpr int eventCount = 10;
        for (int eventIndex = 0; eventIndex < eventCount; ++eventIndex) {
            addEvent(writer, eventIndex, shiftedEvent(eventIndex));
        }

        std::cout << "generated_hepmc3_realistic_sample events=" << eventCount
                  << " output=" << argv[1] << '\n';
    } catch (const std::exception& error) {
        std::cerr << "FAILED: " << error.what() << '\n';
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
