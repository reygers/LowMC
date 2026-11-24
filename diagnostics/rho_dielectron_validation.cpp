#include "lowmc/LowMC.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr double RhoMassGeV = 0.77526;
constexpr double ChargedPionMassGeV = 0.13957039;

double pi() {
    return std::acos(-1.0);
}

double conversionDensity(double mass) {
    return lowmc::krollWadaConversionDensityDMass2(mass);
}

lowmc::LowCurrent rhoBornCurrent() {
    const double pionEnergy = 0.5 * RhoMassGeV;
    const double pionMomentum = std::sqrt(
        pionEnergy * pionEnergy - ChargedPionMassGeV * ChargedPionMassGeV);
    return lowmc::LowCurrent({
        {lowmc::FourVector{pionEnergy, 0.0, 0.0, pionMomentum},
         +1.0, lowmc::InOutSign::Outgoing},
        {lowmc::FourVector{pionEnergy, 0.0, 0.0, -pionMomentum},
         -1.0, lowmc::InOutSign::Outgoing},
    });
}

double numericalDensity(
    const lowmc::LowCurrent& current,
    double mass,
    double momentum,
    lowmc::VirtualCurrentModel model) {
    constexpr int CosThetaBins = 20000;
    const double energy = std::hypot(mass, momentum);
    double angularIntegral = 0.0;
    for (int index = 0; index < CosThetaBins; ++index) {
        const double cosTheta = -1.0 +
            (2.0 * (static_cast<double>(index) + 0.5)) /
                static_cast<double>(CosThetaBins);
        const double sinTheta = std::sqrt(std::max(0.0, 1.0 - cosTheta * cosTheta));
        const lowmc::FourVector virtualPhoton{
            energy,
            momentum * sinTheta,
            0.0,
            momentum * cosTheta,
        };
        angularIntegral += current.intensity(virtualPhoton, model);
    }
    angularIntegral *= 4.0 * pi() / static_cast<double>(CosThetaBins);
    return lowmc::FineStructureConstant / (4.0 * pi() * pi()) *
        momentum * momentum / energy * angularIntegral *
        conversionDensity(mass);
}

double lichardAnalyticDensity(double mass, double momentum) {
    const double rho2 = RhoMassGeV * RhoMassGeV;
    const double mass2 = mass * mass;
    const double energy = std::hypot(mass, momentum);
    const double velocity = std::sqrt(
        1.0 - 4.0 * ChargedPionMassGeV * ChargedPionMassGeV / rho2);
    const double shiftedEnergy = energy + mass2 / RhoMassGeV;
    const double logarithm = std::log(
        (shiftedEnergy + velocity * momentum) /
        (shiftedEnergy - velocity * momentum));
    const double first =
        (1.0 + velocity * velocity - mass2 / rho2) /
        (2.0 * velocity * shiftedEnergy) * logarithm;
    const double denominator =
        energy * energy - velocity * velocity * momentum * momentum +
        (2.0 * energy + mass2 / RhoMassGeV) * mass2 / RhoMassGeV;
    const double second =
        (1.0 - velocity * velocity - mass2 / rho2) * momentum / denominator;
    return 2.0 * lowmc::FineStructureConstant / pi() *
        momentum / energy * conversionDensity(mass) * (first - second);
}

double bppAnalyticDensity(double mass, double momentum) {
    const double rho2 = RhoMassGeV * RhoMassGeV;
    const double energy = std::hypot(mass, momentum);
    const double velocity = std::sqrt(
        1.0 - 4.0 * ChargedPionMassGeV * ChargedPionMassGeV / rho2);
    const double logarithm = std::log(
        (energy + velocity * momentum) /
        (energy - velocity * momentum));
    const double first =
        (1.0 + velocity * velocity) / (2.0 * velocity) * logarithm;
    const double second =
        (1.0 - velocity * velocity) * energy * momentum /
        (energy * energy - velocity * velocity * momentum * momentum);
    return 2.0 * lowmc::FineStructureConstant / pi() *
        momentum / (energy * energy) * conversionDensity(mass) *
        (first - second);
}

double scalarQEDAnalyticDensity(double mass, double momentum) {
    const double rho2 = RhoMassGeV * RhoMassGeV;
    const double pion2 = ChargedPionMassGeV * ChargedPionMassGeV;
    const double mass2 = mass * mass;
    const double energy = std::hypot(mass, momentum);
    const double maximumEnergy = (rho2 + mass2 - 4.0 * pion2) /
        (2.0 * RhoMassGeV);
    if (!(energy < maximumEnergy)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const double xi2 =
        2.0 * RhoMassGeV * (maximumEnergy - energy) /
        (rho2 + mass2 - 2.0 * RhoMassGeV * energy);
    const double xi = std::sqrt(std::max(0.0, xi2));
    const double logarithm = std::log(
        (energy + xi * momentum) / (energy - xi * momentum));
    const double firstCoefficient =
        maximumEnergy * (rho2 - 2.0 * pion2 - 2.0 * energy * RhoMassGeV) +
        mass2 * (maximumEnergy + 0.25 * RhoMassGeV);
    const double secondBracket =
        (RhoMassGeV - 2.0 * maximumEnergy) *
            (2.0 * RhoMassGeV * maximumEnergy - mass2) /
            (2.0 * (energy * energy - xi * xi * momentum * momentum)) -
        2.0 * RhoMassGeV;
    const double braces = firstCoefficient * logarithm -
        energy * momentum * xi * secondBracket;
    const double virtualPhotonDensity =
        4.0 * lowmc::FineStructureConstant / pi() *
        momentum / (energy * energy) /
        std::pow(rho2 - 4.0 * pion2, 1.5) * braces;
    return conversionDensity(mass) * virtualPhotonDensity;
}

double relativeDifference(double value, double reference) {
    return std::abs(value / reference - 1.0);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const std::string outputPath = argc > 1
            ? argv[1]
            : "rho_dielectron_validation.csv";
        if (argc > 2) {
            throw std::invalid_argument(
                "usage: lowmc_rho_dielectron_validation [output.csv]");
        }

        std::ofstream output(outputPath);
        if (!output) {
            throw std::runtime_error("could not open output CSV: " + outputPath);
        }
        output << std::setprecision(17);
        output << "mass_gev,qstar_gev,leading_numeric,bpp_analytic,"
                  "leading_over_bpp,lichard_numeric,lichard_analytic,"
                  "lichard_over_analytic,scalar_qed_analytic,"
                  "leading_over_scalar_qed,lichard_over_scalar_qed\n";

        const lowmc::LowCurrent current = rhoBornCurrent();
        const std::vector<double> masses{0.01, 0.05, 0.10, 0.20};
        const std::vector<double> momenta{0.02, 0.05, 0.10, 0.20};
        double maximumLeadingDifference = 0.0;
        double maximumLichardDifference = 0.0;
        int rows = 0;
        for (double mass : masses) {
            for (double momentum : momenta) {
                const double exact = scalarQEDAnalyticDensity(mass, momentum);
                if (!std::isfinite(exact)) {
                    continue;
                }
                const double leadingNumeric = numericalDensity(
                    current, mass, momentum,
                    lowmc::VirtualCurrentModel::LeadingLow);
                const double leadingAnalytic = bppAnalyticDensity(mass, momentum);
                const double lichardNumeric = numericalDensity(
                    current, mass, momentum,
                    lowmc::VirtualCurrentModel::LichardEq2_8);
                const double lichardAnalytic =
                    lichardAnalyticDensity(mass, momentum);
                maximumLeadingDifference = std::max(
                    maximumLeadingDifference,
                    relativeDifference(leadingNumeric, leadingAnalytic));
                maximumLichardDifference = std::max(
                    maximumLichardDifference,
                    relativeDifference(lichardNumeric, lichardAnalytic));
                output
                    << mass << ',' << momentum << ','
                    << leadingNumeric << ',' << leadingAnalytic << ','
                    << leadingNumeric / leadingAnalytic << ','
                    << lichardNumeric << ',' << lichardAnalytic << ','
                    << lichardNumeric / lichardAnalytic << ','
                    << exact << ','
                    << leadingNumeric / exact << ','
                    << lichardNumeric / exact << '\n';
                ++rows;
            }
        }

        constexpr double Tolerance = 2.0e-7;
        if (maximumLeadingDifference > Tolerance ||
            maximumLichardDifference > Tolerance) {
            throw std::runtime_error("analytic rho validation tolerance exceeded");
        }
        std::cout << "rho_dielectron_validation=ok"
                  << " rows=" << rows
                  << " max_leading_relative_difference="
                  << maximumLeadingDifference
                  << " max_lichard_relative_difference="
                  << maximumLichardDifference
                  << " output=" << outputPath << '\n';
    } catch (const std::exception& error) {
        std::cerr << "FAILED: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
