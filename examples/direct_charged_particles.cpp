#include "lowmc/LowMC.h"

#include <cmath>
#include <iostream>
#include <vector>

int main() {
    const double pionMass = 0.13957039;
    const double momentum = 10.0;
    const double energy = std::sqrt(momentum * momentum + pionMass * pionMass);

    const std::vector<lowmc::ChargedParticle> selectedParticles = {
        {lowmc::FourVector{energy, 0.0, momentum, 0.0}, +1.0, lowmc::InOutSign::Outgoing},
        {lowmc::FourVector{energy, 0.0, -momentum, 0.0}, -1.0, lowmc::InOutSign::Outgoing},
    };

    const lowmc::WeightedEventGeneratorConfig config{
        lowmc::PhotonPhaseSpace{
            0.002,
            0.2,
            -5.0,
            5.0,
            0.0,
            2.0 * std::acos(-1.0),
        },
        1000,
        12345,
        lowmc::tunedMassGridAdaptiveEtaPhiConfig(),
    };

    const lowmc::WeightedEventGenerationResult result =
        lowmc::generateWeightedPhotons(selectedParticles, config);

    std::cout << "selected_particles=" << selectedParticles.size() << '\n';
    std::cout << "weighted_photons=" << result.samples.size() << '\n';
    std::cout << "mean_weight=" << result.weightSummary.meanWeight << '\n';
    std::cout << "effective_sample_size=" << result.weightSummary.effectiveSampleSize << '\n';
    std::cout << "signed_charge_sum=" << result.chargeConservation.netCharge << '\n';

    if (!result.samples.empty()) {
        const lowmc::WeightedEventPhotonSample& sample = result.samples.front();
        const lowmc::Photon& photon = sample.photon;
        std::cout << "first_photon_E=" << photon.momentum.E << '\n';
        std::cout << "first_photon_px=" << photon.momentum.px << '\n';
        std::cout << "first_photon_py=" << photon.momentum.py << '\n';
        std::cout << "first_photon_pz=" << photon.momentum.pz << '\n';
        std::cout << "first_photon_kT=" << sample.kT << '\n';
        std::cout << "first_photon_eta=" << sample.eta << '\n';
        std::cout << "first_photon_phi=" << sample.phi << '\n';
        std::cout << "first_photon_weight=" << photon.weight << '\n';
    }

    return 0;
}
