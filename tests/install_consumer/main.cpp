#include <lowmc/LowMC.h>

#include <cmath>
#include <vector>

int main() {
    const double mass = 0.13957039;
    const double momentum = 1.0;
    const double energy = std::sqrt(momentum * momentum + mass * mass);
    const std::vector<lowmc::ChargedParticle> particles = {
        {lowmc::FourVector{energy, momentum, 0.0, 0.0}, +1.0, lowmc::InOutSign::Outgoing},
        {lowmc::FourVector{energy, -momentum, 0.0, 0.0}, -1.0, lowmc::InOutSign::Outgoing},
    };
    const lowmc::LowCurrent current(particles);
    const lowmc::FourVector photon{0.01, 0.0, 0.01, 0.0};
    const lowmc::ChargeConservationDiagnostic charge = current.chargeConservation();
    const double density = lowmc::lowPhotonDensity(current, photon);
    return charge.conserved && std::isfinite(density) && density >= 0.0 ? 0 : 1;
}
