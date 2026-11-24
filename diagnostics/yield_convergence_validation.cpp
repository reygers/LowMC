#include "lowmc/LowMC.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr double ChargedPionMassGeV = 0.13957039;
constexpr double KTMinGeV = 0.002;
constexpr double KTMaxGeV = 0.2;
constexpr double EtaMax = 5.0;

double pi() {
    return std::acos(-1.0);
}

std::vector<lowmc::ChargedParticle> makeTransverseDipole(double momentum) {
    const double energy = std::hypot(ChargedPionMassGeV, momentum);
    return {
        {lowmc::FourVector{energy, 0.0, momentum, 0.0},
         +1.0, lowmc::InOutSign::Outgoing},
        {lowmc::FourVector{energy, 0.0, -momentum, 0.0},
         -1.0, lowmc::InOutSign::Outgoing},
    };
}

double angularKernel(double velocity, double directionAlongDipole) {
    const double direction2 = directionAlongDipole * directionAlongDipole;
    const double denominator = 1.0 - velocity * velocity * direction2;
    return 4.0 * velocity * velocity * (1.0 - direction2) /
        (denominator * denominator);
}

double exactFiniteEtaYield(double momentum) {
    const double energy = std::hypot(ChargedPionMassGeV, momentum);
    const double velocity = momentum / energy;
    const double bracket =
        (1.0 + velocity * velocity) / (2.0 * velocity) *
            std::log((1.0 + velocity) / (1.0 - velocity)) -
        1.0;
    const double fullAngularIntegral = 8.0 * pi() * bracket;

    // The dipole is transverse to the eta axis.  Subtract the two smooth
    // polar caps excluded by |eta| < EtaMax; the collinear peaks remain in
    // the accepted region, so a modest midpoint grid suffices for the caps.
    constexpr int ZBins = 512;
    constexpr int PhiBins = 1024;
    const double zMin = std::tanh(EtaMax);
    const double zStep = (1.0 - zMin) / static_cast<double>(ZBins);
    const double phiStep = 2.0 * pi() / static_cast<double>(PhiBins);
    double oneCap = 0.0;
    for (int zIndex = 0; zIndex < ZBins; ++zIndex) {
        const double z = zMin + (static_cast<double>(zIndex) + 0.5) * zStep;
        const double transverse = std::sqrt(std::max(0.0, 1.0 - z * z));
        for (int phiIndex = 0; phiIndex < PhiBins; ++phiIndex) {
            const double phi = (static_cast<double>(phiIndex) + 0.5) * phiStep;
            oneCap += angularKernel(velocity, transverse * std::sin(phi));
        }
    }
    oneCap *= zStep * phiStep;
    const double acceptedAngularIntegral = fullAngularIntegral - 2.0 * oneCap;
    return lowmc::FineStructureConstant / (4.0 * pi() * pi()) *
        std::log(KTMaxGeV / KTMinGeV) * acceptedAngularIntegral;
}

lowmc::EventYieldEstimate estimate(double momentum, int level) {
    lowmc::GridAdaptiveEtaPhiProposalConfig grid =
        lowmc::tunedMassGridAdaptiveEtaPhiConfig();
    grid.etaBins = 32;
    grid.phiBins = 32;
    grid.refinementLevels = level;
    grid.refineFraction = 0.1;
    const lowmc::AdaptiveEventYieldIntegrationConfig config{16, 2, 2, grid};
    const lowmc::PhotonPhaseSpace phaseSpace{
        KTMinGeV, KTMaxGeV, -EtaMax, EtaMax, 0.0, 2.0 * pi()};
    return lowmc::estimateFactorizedAdaptiveEventYield(
        makeTransverseDipole(momentum), phaseSpace, config);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const std::string outputPath = argc > 1
            ? argv[1]
            : "yield_convergence_validation.csv";
        if (argc > 2) {
            throw std::invalid_argument(
                "usage: lowmc_yield_convergence_validation [output.csv]");
        }
        std::ofstream output(outputPath);
        if (!output) {
            throw std::runtime_error("could not open output CSV: " + outputPath);
        }
        output << std::setprecision(17);
        output << "momentum_gev,analytic_yield,level,lambda_event,"
                  "true_relative_error,coarse_fine_relative_difference,"
                  "evaluations,grid_cells\n";

        const std::vector<double> momenta{1.0, 10.0, 100.0, 1000.0};
        const std::vector<int> levels{3, 5, 7, 9};
        double maximumLevel9Error = 0.0;
        for (double momentum : momenta) {
            const double analytic = exactFiniteEtaYield(momentum);
            for (int level : levels) {
                const lowmc::EventYieldEstimate result = estimate(momentum, level);
                const double trueError = result.lambdaEvent / analytic - 1.0;
                output << momentum << ',' << analytic << ',' << level << ','
                       << result.lambdaEvent << ',' << trueError << ','
                       << result.relativeDifference << ',' << result.evaluations
                       << ',' << result.gridCells << '\n';
                if (level == 9) {
                    maximumLevel9Error = std::max(
                        maximumLevel9Error, std::abs(trueError));
                }
            }
        }
        if (maximumLevel9Error > 1.0e-2) {
            throw std::runtime_error("level-9 yield error exceeds validation gate");
        }
        std::cout << "yield_convergence_validation=ok"
                  << " max_level9_relative_error=" << maximumLevel9Error
                  << " output=" << outputPath << '\n';
    } catch (const std::exception& error) {
        std::cerr << "FAILED: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
