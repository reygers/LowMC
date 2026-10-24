#include "lowmc/LowMC.h"

#include <cmath>
#include <iostream>
#include <vector>

int main() {
    const double pionMass = 0.13957039;
    const double momentum = 10.0;
    const double energy = std::sqrt(momentum * momentum + pionMass * pionMass);
    const std::vector<lowmc::ChargedParticle> selectedParticles = {
        {lowmc::FourVector{energy, 0.0, momentum, 0.0},
         +1.0,
         lowmc::InOutSign::Outgoing},
        {lowmc::FourVector{energy, 0.0, -momentum, 0.0},
         -1.0,
         lowmc::InOutSign::Outgoing},
    };

    const lowmc::WeightedDielectronGeneratorConfig config{
        lowmc::DielectronPhaseSpace{
            0.002,                    // M_min [GeV]
            0.3,                      // M_max [GeV]
            0.002,                    // qT_min [GeV]
            0.2,                      // qT_max [GeV]
            -1.0,                     // y_min
            1.0,                      // y_max
            0.0,                      // phi_min
            2.0 * std::acos(-1.0),    // phi_max
        },
        1000,
        12345,
        1,
    };

    const lowmc::WeightedDielectronGenerationResult result =
        lowmc::generateWeightedDielectrons(selectedParticles, config);

    std::cout << "Example: generate weighted soft dielectrons\n";
    std::cout << "selected_particles=" << selectedParticles.size() << '\n';
    std::cout << "weighted_dielectrons=" << result.samples.size() << '\n';
    std::cout << "mean_weight=" << result.weightSummary.meanWeight << '\n';
    std::cout << "effective_sample_size="
              << result.weightSummary.effectiveSampleSize << '\n';
    std::cout << "signed_charge_sum="
              << result.chargeConservation.netCharge << '\n';

    if (!result.samples.empty()) {
        const lowmc::WeightedDielectronSample& sample = result.samples.front();
        std::cout << "first_pair_mass=" << sample.mass << '\n';
        std::cout << "first_pair_qT=" << sample.qT << '\n';
        std::cout << "first_pair_rapidity=" << sample.rapidity << '\n';
        std::cout << "first_electron_E=" << sample.dielectron.electron.E << '\n';
        std::cout << "first_positron_E=" << sample.dielectron.positron.E << '\n';
        std::cout << "first_pair_weight=" << sample.dielectron.weight << '\n';
    }

    return 0;
}
