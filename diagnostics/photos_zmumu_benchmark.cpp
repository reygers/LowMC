#include "lowmc/LowMC.h"

#include <HepMC3/FourVector.h>
#include <HepMC3/GenEvent.h>
#include <HepMC3/GenParticle.h>
#include <HepMC3/GenVertex.h>
#include <HepMC3/Units.h>

#include <Photos/Photos.h>
#include <Photos/PhotosHepMC3Event.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr double ZMassGeV = 91.1876;
constexpr double MuonMassGeV = 0.1056583755;
constexpr double Pi = 3.14159265358979323846;

struct Options {
    std::string photonOutput = "photos_zmumu_photons.csv";
    std::string summaryOutput = "photos_zmumu_summary.csv";
    std::string lowmcOutput = "lowmc_zmumu_photons.csv";
    int events = 10000;
    int lowmcSamples = 100000;
    int lowmcRefinementLevels = 9;
    std::uint64_t lowmcSeed = 12345;
    int seed1 = 1802;
    int seed2 = 9373;
    double infraredCutoff = 3.0e-4;
    double kTMin = 0.05;
    double kTMax = 0.2;
    double etaMax = 2.0;
    double energyMax = 0.5 * ZMassGeV;
};

struct BornKinematics {
    double energy = 0.0;
    double momentum = 0.0;
    lowmc::FourVector muPlus;
    lowmc::FourVector muMinus;
    lowmc::LowCurrent current;

    BornKinematics(double energyIn, double momentumIn)
        : energy(energyIn),
          momentum(momentumIn),
          muPlus{energy, momentum, 0.0, 0.0},
          muMinus{energy, -momentum, 0.0, 0.0},
          current(std::vector<lowmc::ChargedParticle>{
              {muPlus, +1.0, lowmc::InOutSign::Outgoing},
              {muMinus, -1.0, lowmc::InOutSign::Outgoing},
          }) {}
};

struct PhotonRecord {
    int eventNumber = 0;
    double energy = 0.0;
    double px = 0.0;
    double py = 0.0;
    double pz = 0.0;
    double kT = 0.0;
    double eta = 0.0;
    double phi = 0.0;
    double thetaPlusBorn = 0.0;
    double thetaMinusBorn = 0.0;
    double thetaPlusPost = 0.0;
    double thetaMinusPost = 0.0;
    double lowmcIntensity = 0.0;
    double analyticIntensity = 0.0;
    double relativeIntensityDifference = 0.0;
    bool accepted = false;
};

void printUsage(std::ostream& output, const char* program) {
    output
        << "Usage: " << program << " [photons.csv] [summary.csv] [options]\n\n"
        << "Run a fixed Z -> mu+ mu- PHOTOS 3.64 pilot and compare the accepted\n"
        << "real-photon rate with the LowMC leading-soft current in a common\n"
        << "rectangular (kT, eta, phi) acceptance. PHOTOS recoil is retained in\n"
        << "the event, while all primary angles use the common Born muons.\n\n"
        << "Options:\n"
        << "  --events N                 Born decays to process (default: 10000)\n"
        << "  --lowmc-output PATH        Weighted LowMC shape CSV\n"
        << "  --lowmc-samples N          Weighted LowMC samples (default: 100000)\n"
        << "  --lowmc-refinement-levels N  Shape-map refinement levels (default: 9)\n"
        << "  --lowmc-seed N             LowMC sampler seed (default: 12345)\n"
        << "  --seed1 N                  PHOTOS first seed, 0..31327 (default: 1802)\n"
        << "  --seed2 N                  PHOTOS second seed, 0..30080 (default: 9373)\n"
        << "  --infrared-cutoff X        Fraction of parent mass (default: 3e-4)\n"
        << "  --kt-min X                 Analysis kT minimum in GeV (default: 0.05)\n"
        << "  --kt-max X                 Analysis kT maximum in GeV (default: 0.2)\n"
        << "  --eta-max X                Symmetric |eta| limit (default: 2)\n"
        << "  --energy-max X             Maximum photon energy in GeV (default: mZ/2)\n";
}

std::string optionValue(int& index, int argc, char** argv, const std::string& name) {
    if (index + 1 >= argc) {
        throw std::invalid_argument(name + " requires a value");
    }
    return argv[++index];
}

int parseInt(const std::string& text, const std::string& name) {
    std::size_t parsed = 0;
    const int value = std::stoi(text, &parsed);
    if (parsed != text.size()) {
        throw std::invalid_argument(name + " must be an integer");
    }
    return value;
}

double parseDouble(const std::string& text, const std::string& name) {
    std::size_t parsed = 0;
    const double value = std::stod(text, &parsed);
    if (parsed != text.size() || !std::isfinite(value)) {
        throw std::invalid_argument(name + " must be a finite number");
    }
    return value;
}

Options parseOptions(int argc, char** argv) {
    Options options;
    int positional = 0;
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument == "--help" || argument == "-h") {
            printUsage(std::cout, argv[0]);
            std::exit(EXIT_SUCCESS);
        }
        if (argument == "--events") {
            options.events = parseInt(optionValue(index, argc, argv, argument), argument);
        } else if (argument == "--lowmc-output") {
            options.lowmcOutput = optionValue(index, argc, argv, argument);
        } else if (argument == "--lowmc-samples") {
            options.lowmcSamples = parseInt(optionValue(index, argc, argv, argument), argument);
        } else if (argument == "--lowmc-refinement-levels") {
            options.lowmcRefinementLevels =
                parseInt(optionValue(index, argc, argv, argument), argument);
        } else if (argument == "--lowmc-seed") {
            options.lowmcSeed = static_cast<std::uint64_t>(
                std::stoull(optionValue(index, argc, argv, argument)));
        } else if (argument == "--seed1") {
            options.seed1 = parseInt(optionValue(index, argc, argv, argument), argument);
        } else if (argument == "--seed2") {
            options.seed2 = parseInt(optionValue(index, argc, argv, argument), argument);
        } else if (argument == "--infrared-cutoff") {
            options.infraredCutoff =
                parseDouble(optionValue(index, argc, argv, argument), argument);
        } else if (argument == "--kt-min") {
            options.kTMin = parseDouble(optionValue(index, argc, argv, argument), argument);
        } else if (argument == "--kt-max") {
            options.kTMax = parseDouble(optionValue(index, argc, argv, argument), argument);
        } else if (argument == "--eta-max") {
            options.etaMax = parseDouble(optionValue(index, argc, argv, argument), argument);
        } else if (argument == "--energy-max") {
            options.energyMax = parseDouble(optionValue(index, argc, argv, argument), argument);
        } else if (!argument.empty() && argument.front() == '-') {
            throw std::invalid_argument("unknown option: " + argument);
        } else if (positional == 0) {
            options.photonOutput = argument;
            ++positional;
        } else if (positional == 1) {
            options.summaryOutput = argument;
            ++positional;
        } else {
            throw std::invalid_argument("too many positional arguments");
        }
    }

    if (options.events <= 0 || options.lowmcSamples <= 0) {
        throw std::invalid_argument("--events and --lowmc-samples must be positive");
    }
    if (options.lowmcRefinementLevels < 0 || options.lowmcRefinementLevels > 12) {
        throw std::invalid_argument("--lowmc-refinement-levels must lie in [0, 12]");
    }
    if (options.seed1 < 0 || options.seed1 > 31327 ||
        options.seed2 < 0 || options.seed2 > 30080) {
        throw std::invalid_argument("PHOTOS seeds are outside their supported ranges");
    }
    if (!(options.infraredCutoff > 0.0 && options.infraredCutoff < 1.0)) {
        throw std::invalid_argument("--infrared-cutoff must lie between zero and one");
    }
    if (!(options.kTMin > ZMassGeV * options.infraredCutoff)) {
        throw std::invalid_argument("--kt-min must exceed the PHOTOS generation cutoff");
    }
    if (!(options.kTMax > options.kTMin) || !(options.etaMax > 0.0)) {
        throw std::invalid_argument("invalid common analysis acceptance");
    }
    if (!(options.energyMax > options.kTMin)) {
        throw std::invalid_argument("--energy-max must exceed --kt-min");
    }
    return options;
}

BornKinematics makeBornKinematics() {
    const double energy = 0.5 * ZMassGeV;
    const double momentum = std::sqrt(energy * energy - MuonMassGeV * MuonMassGeV);
    return {energy, momentum};
}

double integrateEnergyWindow(
    const BornKinematics& born,
    const lowmc::PhotonPhaseSpace& phaseSpace,
    double energyMax,
    const lowmc::AdaptiveEventYieldIntegrationConfig& config,
    int& evaluations,
    int& gridCells) {
    const double referenceKT = std::sqrt(phaseSpace.kTMin * phaseSpace.kTMax);
    const double referenceJacobian = referenceKT * referenceKT;
    const lowmc::GridAdaptiveEtaPhiProposal proposal =
        lowmc::makeGridAdaptiveEtaPhiProposal(
            born.current,
            phaseSpace,
            referenceKT,
            config.gridConfig);
    gridCells = static_cast<int>(proposal.cells.size());

    double integral = 0.0;
    for (const lowmc::GridAdaptiveEtaPhiCell& cell : proposal.cells) {
        const double dEta =
            (cell.etaMax - cell.etaMin) / static_cast<double>(config.etaSubBins);
        const double dPhi =
            (cell.phiMax - cell.phiMin) / static_cast<double>(config.phiSubBins);
        for (int etaIndex = 0; etaIndex < config.etaSubBins; ++etaIndex) {
            const double eta =
                cell.etaMin + (static_cast<double>(etaIndex) + 0.5) * dEta;
            const double allowedKTMax =
                std::min(phaseSpace.kTMax, energyMax / std::cosh(eta));
            if (allowedKTMax <= phaseSpace.kTMin) {
                continue;
            }
            const double logKTInterval = std::log(allowedKTMax / phaseSpace.kTMin);
            for (int phiIndex = 0; phiIndex < config.phiSubBins; ++phiIndex) {
                const double phi =
                    cell.phiMin + (static_cast<double>(phiIndex) + 0.5) * dPhi;
                const lowmc::FourVector photon =
                    lowmc::photonMomentumFromKTEtaPhi(referenceKT, eta, phi);
                integral += lowmc::lowPhotonDensity(born.current, photon) *
                            referenceJacobian * logKTInterval * dEta * dPhi;
                ++evaluations;
            }
        }
    }
    return integral;
}

lowmc::EventYieldEstimate estimateEnergyWindowYield(
    const BornKinematics& born,
    const lowmc::PhotonPhaseSpace& phaseSpace,
    double energyMax) {
    lowmc::AdaptiveEventYieldIntegrationConfig fine;
    fine.gridConfig = lowmc::tunedMassGridAdaptiveEtaPhiConfig();
    fine.gridConfig.refinementLevels = 9;
    lowmc::AdaptiveEventYieldIntegrationConfig coarse = fine;
    coarse.gridConfig.refinementLevels = 8;

    lowmc::EventYieldEstimate estimate;
    estimate.method = "energy-window-factorized-adaptive-eta-phi";
    estimate.logKTBins = 0;
    estimate.etaBins = fine.etaSubBins;
    estimate.phiBins = fine.phiSubBins;
    estimate.lambdaEvent = integrateEnergyWindow(
        born,
        phaseSpace,
        energyMax,
        fine,
        estimate.evaluations,
        estimate.gridCells);
    int coarseEvaluations = 0;
    int coarseGridCells = 0;
    estimate.coarseLambdaEvent = integrateEnergyWindow(
        born,
        phaseSpace,
        energyMax,
        coarse,
        coarseEvaluations,
        coarseGridCells);
    estimate.absoluteDifference =
        std::abs(estimate.lambdaEvent - estimate.coarseLambdaEvent);
    estimate.relativeDifference =
        estimate.lambdaEvent > 0.0 ? estimate.absoluteDifference / estimate.lambdaEvent : 0.0;
    estimate.status =
        std::isfinite(estimate.lambdaEvent) && estimate.lambdaEvent >= 0.0 ? "ok" : "invalid";
    return estimate;
}

HepMC3::GenEvent makeBornEvent(const BornKinematics& born, int eventNumber) {
    HepMC3::GenEvent event(HepMC3::Units::GEV, HepMC3::Units::MM);
    event.set_event_number(eventNumber);

    auto z = std::make_shared<HepMC3::GenParticle>(
        HepMC3::FourVector{0.0, 0.0, 0.0, ZMassGeV}, 23, 2);
    z->set_generated_mass(ZMassGeV);
    auto muPlus = std::make_shared<HepMC3::GenParticle>(
        HepMC3::FourVector{born.momentum, 0.0, 0.0, born.energy}, -13, 1);
    muPlus->set_generated_mass(MuonMassGeV);
    auto muMinus = std::make_shared<HepMC3::GenParticle>(
        HepMC3::FourVector{-born.momentum, 0.0, 0.0, born.energy}, 13, 1);
    muMinus->set_generated_mass(MuonMassGeV);

    auto decay = std::make_shared<HepMC3::GenVertex>();
    decay->add_particle_in(z);
    decay->add_particle_out(muPlus);
    decay->add_particle_out(muMinus);
    event.add_vertex(decay);
    return event;
}

double magnitude(double x, double y, double z) {
    return std::sqrt(x * x + y * y + z * z);
}

double openingAngle(
    double px,
    double py,
    double pz,
    const lowmc::FourVector& direction) {
    const double left = magnitude(px, py, pz);
    const double right = magnitude(direction.px, direction.py, direction.pz);
    const double cosine = std::clamp(
        (px * direction.px + py * direction.py + pz * direction.pz) / (left * right),
        -1.0,
        1.0);
    return std::acos(cosine);
}

double analyticTwoMuonIntensity(
    const BornKinematics& born,
    const lowmc::FourVector& photon) {
    const double plusDenominator = born.muPlus.dot(photon);
    const double minusDenominator = born.muMinus.dot(photon);
    const double cross = born.muPlus.dot(born.muMinus);
    return 2.0 * cross / (plusDenominator * minusDenominator) -
           MuonMassGeV * MuonMassGeV /
               (plusDenominator * plusDenominator) -
           MuonMassGeV * MuonMassGeV /
               (minusDenominator * minusDenominator);
}

PhotonRecord makePhotonRecord(
    int eventNumber,
    const HepMC3::FourVector& momentum,
    const BornKinematics& born,
    const HepMC3::FourVector& postMuPlus,
    const HepMC3::FourVector& postMuMinus,
    const Options& options) {
    PhotonRecord record;
    record.eventNumber = eventNumber;
    record.energy = momentum.e();
    record.px = momentum.px();
    record.py = momentum.py();
    record.pz = momentum.pz();
    record.kT = std::hypot(record.px, record.py);
    record.eta = record.kT > 0.0 ? std::asinh(record.pz / record.kT)
                                 : std::copysign(std::numeric_limits<double>::infinity(), record.pz);
    record.phi = std::atan2(record.py, record.px);
    record.thetaPlusBorn = openingAngle(record.px, record.py, record.pz, born.muPlus);
    record.thetaMinusBorn = openingAngle(record.px, record.py, record.pz, born.muMinus);
    const lowmc::FourVector postPlus{
        postMuPlus.e(), postMuPlus.px(), postMuPlus.py(), postMuPlus.pz()};
    const lowmc::FourVector postMinus{
        postMuMinus.e(), postMuMinus.px(), postMuMinus.py(), postMuMinus.pz()};
    record.thetaPlusPost = openingAngle(record.px, record.py, record.pz, postPlus);
    record.thetaMinusPost = openingAngle(record.px, record.py, record.pz, postMinus);

    const lowmc::FourVector lowPhoton{record.energy, record.px, record.py, record.pz};
    record.lowmcIntensity = born.current.intensity(lowPhoton);
    record.analyticIntensity = analyticTwoMuonIntensity(born, lowPhoton);
    const double scale = std::max(std::abs(record.analyticIntensity), 1.0e-300);
    record.relativeIntensityDifference =
        std::abs(record.lowmcIntensity - record.analyticIntensity) / scale;
    record.accepted =
        record.kT >= options.kTMin && record.kT < options.kTMax &&
        std::abs(record.eta) < options.etaMax &&
        record.energy < options.energyMax &&
        record.phi >= -Pi && record.phi < Pi;
    return record;
}

double stableMomentumResidual(const HepMC3::GenEvent& event) {
    double energy = -ZMassGeV;
    double px = 0.0;
    double py = 0.0;
    double pz = 0.0;
    for (const HepMC3::ConstGenParticlePtr& particle : event.particles()) {
        if (particle->status() != 1) {
            continue;
        }
        energy += particle->momentum().e();
        px += particle->momentum().px();
        py += particle->momentum().py();
        pz += particle->momentum().pz();
    }
    return std::max({std::abs(energy), std::abs(px), std::abs(py), std::abs(pz)});
}

void writePhotonHeader(std::ofstream& output) {
    output
        << "event,accepted,E,px,py,pz,kT,eta,phi,theta_plus_born,theta_minus_born,"
        << "theta_plus_post,theta_minus_post,"
        << "lowmc_intensity,analytic_intensity,relative_intensity_difference\n";
}

void writePhoton(std::ofstream& output, const PhotonRecord& record) {
    output << record.eventNumber << ',' << (record.accepted ? 1 : 0) << ','
           << record.energy << ',' << record.px << ',' << record.py << ',' << record.pz << ','
           << record.kT << ',' << record.eta << ',' << record.phi << ','
           << record.thetaPlusBorn << ',' << record.thetaMinusBorn << ','
           << record.thetaPlusPost << ',' << record.thetaMinusPost << ','
           << record.lowmcIntensity << ',' << record.analyticIntensity << ','
           << record.relativeIntensityDifference << '\n';
}

void configurePhotos(const Options& options) {
    Photospp::Photos::setMomentumUnit(Photospp::Photos::GEV);
    Photospp::Photos::setSeed(options.seed1, options.seed2);
    Photospp::Photos::setExponentiation(false);
    Photospp::Photos::setDoubleBrem(false);
    Photospp::Photos::setQuatroBrem(false);
    Photospp::Photos::setPairEmission(false);
    Photospp::Photos::setPhotonEmission(true);
    Photospp::Photos::setInterference(true);
    Photospp::Photos::setCorrectionWtForW(false);
    Photospp::Photos::setMeCorrectionWtForScalar(false);
    Photospp::Photos::setMeCorrectionWtForW(false);
    Photospp::Photos::setMeCorrectionWtForZ(false);
    Photospp::Photos::setAlphaQED(lowmc::FineStructureConstant);
    Photospp::Photos::setInfraredCutOff(options.infraredCutoff);
    Photospp::Photos::forceMassFrom4Vector(true);
    Photospp::Photos::forceMassFromEventRecord(13);
    Photospp::Photos::initialize();
}

lowmc::WeightedEventGenerationResult writeLowMCSamples(
    const Options& options,
    const BornKinematics& born,
    const lowmc::PhotonPhaseSpace& phaseSpace) {
    lowmc::WeightedEventGeneratorConfig config;
    config.phaseSpace = phaseSpace;
    config.samples = options.lowmcSamples;
    config.seed = options.lowmcSeed;
    config.gridConfig = lowmc::tunedMassGridAdaptiveEtaPhiConfig();
    config.gridConfig.refinementLevels = options.lowmcRefinementLevels;

    const lowmc::WeightedEventGenerationResult result =
        lowmc::generateWeightedPhotons(born.current.particles(), config);
    std::ofstream output(options.lowmcOutput);
    if (!output) {
        throw std::runtime_error("could not open LowMC CSV: " + options.lowmcOutput);
    }
    output << std::setprecision(17)
           << "sample,E,px,py,pz,kT,eta,phi,theta_plus_born,theta_minus_born,"
           << "theta_nearest_born,weight,proposal_density\n";
    for (std::size_t index = 0; index < result.samples.size(); ++index) {
        const lowmc::WeightedEventPhotonSample& sample = result.samples[index];
        const lowmc::FourVector& photon = sample.photon.momentum;
        const double thetaPlus =
            openingAngle(photon.px, photon.py, photon.pz, born.muPlus);
        const double thetaMinus =
            openingAngle(photon.px, photon.py, photon.pz, born.muMinus);
        output << index << ',' << photon.E << ',' << photon.px << ',' << photon.py << ','
               << photon.pz << ',' << sample.kT << ',' << sample.eta << ',' << sample.phi << ','
               << thetaPlus << ',' << thetaMinus << ',' << std::min(thetaPlus, thetaMinus) << ','
               << sample.photon.weight << ',' << sample.proposalDensity << '\n';
    }
    return result;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Options options = parseOptions(argc, argv);
        const BornKinematics born = makeBornKinematics();
        const lowmc::ChargeConservationDiagnostic charge = born.current.chargeConservation();
        if (!charge.conserved) {
            throw std::runtime_error("fixed Born current is not charge conserving");
        }

        const lowmc::PhotonPhaseSpace phaseSpace{
            options.kTMin,
            options.kTMax,
            -options.etaMax,
            options.etaMax,
            -Pi,
            Pi,
        };
        const lowmc::EventYieldEstimate lowmcYield =
            estimateEnergyWindowYield(born, phaseSpace, options.energyMax);
        const lowmc::WeightedEventGenerationResult lowmcShapes =
            writeLowMCSamples(options, born, phaseSpace);

        std::ofstream photons(options.photonOutput);
        if (!photons) {
            throw std::runtime_error("could not open photon CSV: " + options.photonOutput);
        }
        photons << std::setprecision(17);
        writePhotonHeader(photons);

        configurePhotos(options);

        int zeroPhotonEvents = 0;
        int onePhotonEvents = 0;
        int multiphotonEvents = 0;
        int totalPhotons = 0;
        int acceptedPhotons = 0;
        double maxMomentumResidual = 0.0;
        double maxIntensityDifference = 0.0;

        for (int eventNumber = 0; eventNumber < options.events; ++eventNumber) {
            HepMC3::GenEvent event = makeBornEvent(born, eventNumber);
            Photospp::Photos::setEventNo(eventNumber);
            Photospp::PhotosHepMC3Event photosEvent(&event);
            photosEvent.process();

            HepMC3::FourVector postMuPlus;
            HepMC3::FourVector postMuMinus;
            bool foundPostMuPlus = false;
            bool foundPostMuMinus = false;
            for (const HepMC3::GenParticlePtr& particle : event.particles()) {
                if (particle->status() == 1 && particle->pid() == -13) {
                    postMuPlus = particle->momentum();
                    foundPostMuPlus = true;
                } else if (particle->status() == 1 && particle->pid() == 13) {
                    postMuMinus = particle->momentum();
                    foundPostMuMinus = true;
                }
            }
            if (!foundPostMuPlus || !foundPostMuMinus) {
                throw std::runtime_error("PHOTOS event has no final muon pair");
            }

            int eventPhotons = 0;
            for (const HepMC3::GenParticlePtr& particle : event.particles()) {
                if (particle->pid() != 22 || particle->status() != 1) {
                    continue;
                }
                ++eventPhotons;
                const PhotonRecord record =
                    makePhotonRecord(
                        eventNumber,
                        particle->momentum(),
                        born,
                        postMuPlus,
                        postMuMinus,
                        options);
                writePhoton(photons, record);
                if (record.accepted) {
                    ++acceptedPhotons;
                }
                maxIntensityDifference =
                    std::max(maxIntensityDifference, record.relativeIntensityDifference);
            }
            totalPhotons += eventPhotons;
            if (eventPhotons == 0) {
                ++zeroPhotonEvents;
            } else if (eventPhotons == 1) {
                ++onePhotonEvents;
            } else {
                ++multiphotonEvents;
            }
            maxMomentumResidual =
                std::max(maxMomentumResidual, stableMomentumResidual(event));
        }

        const double photosRate =
            static_cast<double>(acceptedPhotons) / static_cast<double>(options.events);
        const double photosRateError =
            std::sqrt(photosRate / static_cast<double>(options.events));
        const double rateRatio = lowmcYield.lambdaEvent > 0.0
                                     ? photosRate / lowmcYield.lambdaEvent
                                     : std::numeric_limits<double>::quiet_NaN();

        std::ofstream summary(options.summaryOutput);
        if (!summary) {
            throw std::runtime_error("could not open summary CSV: " + options.summaryOutput);
        }
        summary << std::setprecision(17);
        summary
            << "photos_version,z_mass_gev,muon_mass_gev,born_muon_energy_gev,"
            << "born_muon_momentum_gev,alpha,events,seed1,seed2,lowmc_shape_samples,"
            << "lowmc_shape_seed,lowmc_shape_refinement_levels,lowmc_shape_mean_weight,"
            << "lowmc_shape_ess_fraction,interference,"
            << "exponentiation,double_brem,quatro_brem,pair_emission,me_correction_z,"
            << "infrared_cutoff_fraction,infrared_cutoff_gev,kt_min,kt_max,eta_min,eta_max,"
            << "energy_max_gev,"
            << "zero_photon_events,one_photon_events,multiphoton_events,total_photons,"
            << "accepted_photons,photos_accepted_rate,photos_accepted_rate_error,"
            << "lowmc_lambda_event,photos_over_lowmc,lowmc_over_photos,"
            << "lowmc_relative_integration_difference,max_momentum_residual,"
            << "max_current_relative_difference\n";
        summary
            << "3.64," << ZMassGeV << ',' << MuonMassGeV << ',' << born.energy << ','
            << born.momentum << ',' << lowmc::FineStructureConstant << ',' << options.events << ','
            << options.seed1 << ',' << options.seed2 << ',' << options.lowmcSamples << ','
            << options.lowmcSeed << ',' << options.lowmcRefinementLevels << ','
            << lowmcShapes.weightSummary.meanWeight << ','
            << lowmcShapes.weightSummary.effectiveSampleSize /
                   static_cast<double>(lowmcShapes.weightSummary.count)
            << ",1,0,0,0,0,0,"
            << options.infraredCutoff << ','
            << options.infraredCutoff * ZMassGeV << ',' << options.kTMin << ','
            << options.kTMax << ',' << -options.etaMax << ',' << options.etaMax << ','
            << options.energyMax << ','
            << zeroPhotonEvents << ',' << onePhotonEvents << ',' << multiphotonEvents << ','
            << totalPhotons << ',' << acceptedPhotons << ',' << photosRate << ','
            << photosRateError << ',' << lowmcYield.lambdaEvent << ',' << rateRatio << ','
            << (photosRate > 0.0 ? lowmcYield.lambdaEvent / photosRate
                                 : std::numeric_limits<double>::quiet_NaN()) << ','
            << lowmcYield.relativeDifference << ',' << maxMomentumResidual << ','
            << maxIntensityDifference << '\n';

        const bool ok =
            multiphotonEvents == 0 && maxMomentumResidual < 1.0e-6 &&
            maxIntensityDifference < 1.0e-8 && lowmcYield.status == std::string("ok");
        std::cout << std::setprecision(8)
                  << "status=" << (ok ? "ok" : "failed")
                  << " events=" << options.events
                  << " zero_photon_events=" << zeroPhotonEvents
                  << " one_photon_events=" << onePhotonEvents
                  << " multiphoton_events=" << multiphotonEvents
                  << " accepted_photons=" << acceptedPhotons
                  << " photos_rate=" << photosRate
                  << " lowmc_lambda_event=" << lowmcYield.lambdaEvent
                  << " photos_over_lowmc=" << rateRatio
                  << " max_momentum_residual=" << maxMomentumResidual
                  << " max_current_relative_difference=" << maxIntensityDifference
                  << " photons_csv=" << options.photonOutput
                  << " lowmc_csv=" << options.lowmcOutput
                  << " summary_csv=" << options.summaryOutput << '\n';
        return ok ? EXIT_SUCCESS : EXIT_FAILURE;
    } catch (const std::exception& error) {
        std::cerr << "FAILED: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
