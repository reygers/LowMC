#include "lowmc/LowMC.h"
#include "lowmc/pythia8/Pythia8Adapter.h"

#include <Pythia8/Pythia.h>

#include <cmath>
#include <cstdint>
#include <iostream>
#include <stdexcept>

#ifndef LOWMC_PYTHIA8_XMLDOC
#define LOWMC_PYTHIA8_XMLDOC "../share/Pythia8/xmldoc"
#endif

int main() {
    // Generate one hadronic Z decay with Pythia's own QED shower disabled so
    // that this example has one clearly identified photon source: LowMC.
    Pythia8::Pythia pythia(LOWMC_PYTHIA8_XMLDOC, false);
    pythia.readString("Print:quiet = on");
    pythia.readString("Random:setSeed = on");
    pythia.readString("Random:seed = 12345");
    pythia.readString("PDF:lepton = off");
    pythia.readString("WeakSingleBoson:ffbar2gmZ = on");
    pythia.readString("23:onMode = off");
    pythia.readString("23:onIfAny = 1 2 3 4 5");
    pythia.readString("Beams:idA = 11");
    pythia.readString("Beams:idB = -11");
    pythia.readString("Beams:eCM = 91.1876");
    pythia.readString("TimeShower:QEDshowerByQ = off");
    pythia.readString("TimeShower:QEDshowerByL = off");
    pythia.readString("SpaceShower:QEDshowerByQ = off");
    pythia.readString("SpaceShower:QEDshowerByL = off");
    pythia.readString("HadronLevel:QED = off");
    if (!pythia.init()) {
        throw std::runtime_error("Pythia initialization failed");
    }
    if (!pythia.next()) {
        throw std::runtime_error("Pythia did not generate the requested event");
    }

    const lowmc::AdapterParticleSelectionResult selection =
        lowmc::pythia8::selectChargedParticles(pythia.event);
    const bool chargeConserved =
        std::abs(selection.diagnostic.signedChargeSum) <= 1.0e-12;
    if (!chargeConserved) {
        throw std::runtime_error("selected Pythia event is not charge conserving");
    }

    const lowmc::WeightedEventGeneratorConfig config{
        lowmc::PhotonPhaseSpace{
            0.002,
            0.2,
            -5.0,
            5.0,
            0.0,
            2.0 * std::acos(-1.0),
        },
        100,
        std::uint64_t{24680},
        lowmc::tunedMassGridAdaptiveEtaPhiConfig(),
    };
    const lowmc::WeightedEventGenerationResult photons =
        lowmc::generateWeightedPhotons(selection.selectedParticles, config);

    std::cout << "Example: generate weighted Low photons from a Pythia event\n";
    std::cout << "selected_incoming_charged="
              << selection.diagnostic.selectedIncomingCharged << '\n';
    std::cout << "selected_outgoing_charged="
              << selection.diagnostic.selectedOutgoingCharged << '\n';
    std::cout << "charge_conserved=" << chargeConserved << '\n';
    std::cout << "weighted_photons=" << photons.samples.size() << '\n';
    std::cout << "mean_weight=" << photons.weightSummary.meanWeight << '\n';

    return photons.samples.empty() ? 1 : 0;
}
