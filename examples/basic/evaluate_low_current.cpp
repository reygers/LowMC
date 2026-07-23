#include "lowmc/LowMC.h"

#include <cmath>
#include <iostream>
#include <vector>

int main() {
    // A charge-neutral pair of outgoing pions.
    const double pionMass = 0.13957039;
    const double momentum = 1.0;
    const double energy = std::sqrt(momentum * momentum + pionMass * pionMass);
    const std::vector<lowmc::ChargedParticle> particles = {
        {lowmc::FourVector{energy, momentum, 0.0, 0.0}, +1.0, lowmc::InOutSign::Outgoing},
        {lowmc::FourVector{energy, -momentum, 0.0, 0.0}, -1.0, lowmc::InOutSign::Outgoing},
    };

    // A 20 MeV photon pointing along +y.
    const lowmc::FourVector photon{0.020, 0.0, 0.020, 0.0};
    const lowmc::LowCurrent lowCurrent(particles);
    const lowmc::ChargeConservationDiagnostic charge =
        lowCurrent.chargeConservation();
    const lowmc::FourVector current = lowCurrent.current(photon);
    const double intensity = lowCurrent.intensity(photon);
    const double density = lowmc::lowPhotonDensity(lowCurrent, photon);

    std::cout << "Example: evaluate the Low current at one photon momentum\n";
    std::cout << "charge_conserved=" << charge.conserved << '\n';
    std::cout << "signed_charge_sum=" << charge.netCharge << '\n';
    std::cout << "current=(" << current.E << ", " << current.px << ", "
              << current.py << ", " << current.pz << ")\n";
    std::cout << "intensity=" << intensity << '\n';
    std::cout << "low_photon_density=" << density << '\n';

    return charge.conserved && std::isfinite(density) && density >= 0.0 ? 0 : 1;
}
