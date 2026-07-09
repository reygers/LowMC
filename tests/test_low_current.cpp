#include "lowmc/LowMC.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "FAILED: " << message << '\n';
        std::exit(EXIT_FAILURE);
    }
}

void requireNear(double actual, double expected, double tolerance, const std::string& message) {
    if (std::abs(actual - expected) > tolerance) {
        std::cerr << "FAILED: " << message << "\n  actual: " << actual
                  << "\n  expected: " << expected << "\n  tolerance: " << tolerance << '\n';
        std::exit(EXIT_FAILURE);
    }
}

void requireThrowsInvalidArgument(void (*action)(), const std::string& message) {
    try {
        action();
    } catch (const std::invalid_argument&) {
        return;
    } catch (...) {
        std::cerr << "FAILED: " << message << "\n  threw unexpected exception type\n";
        std::exit(EXIT_FAILURE);
    }

    std::cerr << "FAILED: " << message << "\n  did not throw\n";
    std::exit(EXIT_FAILURE);
}

lowmc::FourVector masslessPhotonFromEnergyAndDirection(
    double energy, double nx, double ny, double nz) {
    return lowmc::FourVector{energy, energy * nx, energy * ny, energy * nz};
}

void testChargeNeutralTwoParticleCurrentIsTransverse() {
    const std::vector<lowmc::ChargedParticle> particles = {
        {lowmc::FourVector{2.0, 0.0, 0.0, 1.0}, +1.0, lowmc::InOutSign::Outgoing},
        {lowmc::FourVector{2.0, 0.0, 0.0, -1.0}, -1.0, lowmc::InOutSign::Outgoing},
    };

    const lowmc::LowCurrent current(particles);
    const lowmc::FourVector photon = masslessPhotonFromEnergyAndDirection(
        0.1, 1.0 / std::sqrt(2.0), 0.0, 1.0 / std::sqrt(2.0));

    const lowmc::ChargeConservationDiagnostic diagnostic = current.chargeConservation();
    require(diagnostic.conserved, "opposite outgoing charges should conserve selected charge");
    requireNear(diagnostic.netCharge, 0.0, 1.0e-12, "net selected charge");

    const lowmc::FourVector J = current.current(photon);
    requireNear(photon.dot(J), 0.0, 1.0e-12, "charge-conserving current should satisfy k.J=0");
}

void testFiniteMassCollinearRegularization() {
    const double mass = 0.5;
    const double pz = 10.0;
    const double energy = std::sqrt(pz * pz + mass * mass);
    const lowmc::FourVector particle{energy, 0.0, 0.0, pz};
    const lowmc::FourVector collinearPhoton{0.1, 0.0, 0.0, 0.1};

    const double denominator = particle.dot(collinearPhoton);
    require(denominator > 0.0, "massive collinear p.k should remain positive");
    require(std::isfinite(denominator), "massive collinear p.k should be finite");

    const lowmc::LowCurrent current({
        {particle, +1.0, lowmc::InOutSign::Outgoing},
    });
    const double intensity = current.intensity(collinearPhoton);
    require(std::isfinite(intensity), "massive collinear intensity should be finite");
}

void testCurrentSquaredMatchesPairSum() {
    const std::vector<lowmc::ChargedParticle> particles = {
        {lowmc::FourVector{5.0, 0.0, 0.0, 5.0}, +1.0, lowmc::InOutSign::Incoming},
        {lowmc::FourVector{5.0, 0.0, 0.0, -5.0}, -1.0, lowmc::InOutSign::Incoming},
        {lowmc::FourVector{3.0, 1.0, 0.4, 2.0}, +1.0, lowmc::InOutSign::Outgoing},
        {lowmc::FourVector{2.5, -0.7, -0.2, -1.2}, -1.0, lowmc::InOutSign::Outgoing},
    };

    const lowmc::LowCurrent current(particles);
    const lowmc::FourVector photon = masslessPhotonFromEnergyAndDirection(
        0.2, 0.3, 0.4, std::sqrt(1.0 - 0.3 * 0.3 - 0.4 * 0.4));

    requireNear(
        current.intensity(photon),
        current.intensityPairSum(photon),
        1.0e-12,
        "current-squared intensity should match explicit pair sum");
}

void testLowPhotonDensityUsesLowNormalization() {
    const std::vector<lowmc::ChargedParticle> particles = {
        {lowmc::FourVector{2.0, 0.0, 0.0, 1.0}, +1.0, lowmc::InOutSign::Outgoing},
        {lowmc::FourVector{2.0, 0.0, 0.0, -1.0}, -1.0, lowmc::InOutSign::Outgoing},
    };

    const lowmc::LowCurrent current(particles);
    const lowmc::FourVector photon = masslessPhotonFromEnergyAndDirection(
        0.1, 1.0 / std::sqrt(2.0), 0.0, 1.0 / std::sqrt(2.0));
    const double alpha = 0.01;
    const double twoPi = 2.0 * std::acos(-1.0);

    requireNear(
        lowmc::lowPhotonDensity(current, photon, alpha),
        alpha / (twoPi * twoPi) * current.intensity(photon),
        1.0e-12,
        "low photon density should be alpha/(2pi)^2 times S(k)");
}

void testLowPhotonDensityDefaultAlphaIsPositive() {
    const lowmc::LowCurrent current({
        {lowmc::FourVector{2.0, 0.0, 0.0, 1.0}, +1.0, lowmc::InOutSign::Outgoing},
        {lowmc::FourVector{2.0, 0.0, 0.0, -1.0}, -1.0, lowmc::InOutSign::Outgoing},
    });
    const lowmc::FourVector photon = masslessPhotonFromEnergyAndDirection(
        0.1, 1.0 / std::sqrt(2.0), 0.0, 1.0 / std::sqrt(2.0));

    const double density = lowmc::lowPhotonDensity(current, photon);
    require(density > 0.0, "default-alpha low photon density should be positive");
    require(std::isfinite(density), "default-alpha low photon density should be finite");
}

lowmc::PhotonPhaseSpace testPhaseSpace() {
    return lowmc::PhotonPhaseSpace{
        0.02,
        0.2,
        -1.5,
        1.5,
        0.0,
        2.0 * std::acos(-1.0),
    };
}

void testWeightedPhotonSampleUsesProposalDensityConvention() {
    const lowmc::LowCurrent current({
        {lowmc::FourVector{2.0, 0.0, 0.0, 1.0}, +1.0, lowmc::InOutSign::Outgoing},
        {lowmc::FourVector{2.0, 0.0, 0.0, -1.0}, -1.0, lowmc::InOutSign::Outgoing},
    });
    const lowmc::FourVector photon = masslessPhotonFromEnergyAndDirection(
        0.1, 1.0 / std::sqrt(2.0), 0.0, 1.0 / std::sqrt(2.0));
    const double proposalDensity = 3.0;
    const double alpha = 0.01;

    const lowmc::WeightedPhotonSample sample = lowmc::makeWeightedPhotonSample(
        current,
        testPhaseSpace(),
        lowmc::WeightedPhotonProposal{photon, proposalDensity},
        alpha);

    requireNear(
        sample.photon.weight,
        lowmc::lowPhotonDensity(current, photon, alpha) / proposalDensity,
        1.0e-12,
        "weighted photon should use density divided by q_mu");
    requireNear(
        sample.proposalDensity,
        proposalDensity,
        1.0e-12,
        "weighted photon sample should retain proposal density");
    requireNear(
        sample.photon.momentum.E,
        photon.E,
        1.0e-12,
        "weighted photon sample should retain photon momentum");
}

void makeSampleWithZeroProposalDensity() {
    const lowmc::LowCurrent current({
        {lowmc::FourVector{2.0, 0.0, 0.0, 1.0}, +1.0, lowmc::InOutSign::Outgoing},
        {lowmc::FourVector{2.0, 0.0, 0.0, -1.0}, -1.0, lowmc::InOutSign::Outgoing},
    });
    const lowmc::FourVector photon = masslessPhotonFromEnergyAndDirection(
        0.1, 1.0 / std::sqrt(2.0), 0.0, 1.0 / std::sqrt(2.0));

    static_cast<void>(lowmc::makeWeightedPhotonSample(
        current,
        testPhaseSpace(),
        lowmc::WeightedPhotonProposal{photon, 0.0}));
}

void makeSampleWithInvalidPhaseSpace() {
    const lowmc::LowCurrent current({
        {lowmc::FourVector{2.0, 0.0, 0.0, 1.0}, +1.0, lowmc::InOutSign::Outgoing},
        {lowmc::FourVector{2.0, 0.0, 0.0, -1.0}, -1.0, lowmc::InOutSign::Outgoing},
    });
    const lowmc::FourVector photon = masslessPhotonFromEnergyAndDirection(
        0.1, 1.0 / std::sqrt(2.0), 0.0, 1.0 / std::sqrt(2.0));
    lowmc::PhotonPhaseSpace phaseSpace = testPhaseSpace();
    phaseSpace.kTMax = phaseSpace.kTMin;

    static_cast<void>(lowmc::makeWeightedPhotonSample(
        current,
        phaseSpace,
        lowmc::WeightedPhotonProposal{photon, 1.0}));
}

void testLogUniformKTSamplingEndpoints() {
    const lowmc::PhotonPhaseSpace phaseSpace = testPhaseSpace();

    const lowmc::LogUniformKTProposalPoint lower =
        lowmc::sampleLogUniformKT(phaseSpace, 0.0);
    requireNear(
        lower.kT,
        phaseSpace.kTMin,
        1.0e-15,
        "log-uniform kT sampler should map u=0 to kTMin");

    const lowmc::LogUniformKTProposalPoint middle =
        lowmc::sampleLogUniformKT(phaseSpace, 0.5);
    requireNear(
        middle.kT,
        std::sqrt(phaseSpace.kTMin * phaseSpace.kTMax),
        1.0e-15,
        "log-uniform kT sampler should map u=0.5 to geometric mean");

    const double nearOne = std::nextafter(1.0, 0.0);
    const lowmc::LogUniformKTProposalPoint upper =
        lowmc::sampleLogUniformKT(phaseSpace, nearOne);
    require(
        upper.kT < phaseSpace.kTMax,
        "log-uniform kT sampler should keep u in [0,1) below kTMax");
    require(
        upper.kT > phaseSpace.kTMin,
        "log-uniform kT sampler near u=1 should remain inside range");
}

void testLogUniformKTDensityFormulas() {
    const lowmc::PhotonPhaseSpace phaseSpace = testPhaseSpace();
    const double kT = 0.05;
    const double normalization = std::log(phaseSpace.kTMax / phaseSpace.kTMin);

    requireNear(
        lowmc::logUniformKTDensityDKT(phaseSpace, kT),
        1.0 / (kT * normalization),
        1.0e-12,
        "log-uniform kT density with respect to dkT");
    requireNear(
        lowmc::logUniformKTDensityPhysicalMeasure(phaseSpace, kT),
        1.0 / (kT * kT * normalization),
        1.0e-12,
        "log-uniform kT density with respect to kT dkT");

    const lowmc::LogUniformKTProposalPoint proposal =
        lowmc::sampleLogUniformKT(phaseSpace, 0.25);
    requireNear(
        proposal.densityDKT,
        lowmc::logUniformKTDensityDKT(phaseSpace, proposal.kT),
        1.0e-12,
        "sampled log-uniform kT point should include dkT density");
    requireNear(
        proposal.densityPhysicalMeasure,
        lowmc::logUniformKTDensityPhysicalMeasure(phaseSpace, proposal.kT),
        1.0e-12,
        "sampled log-uniform kT point should include physical-measure density");
}

void sampleLogUniformKTWithUnitOne() {
    static_cast<void>(lowmc::sampleLogUniformKT(testPhaseSpace(), 1.0));
}

void densityWithKTOutsideRange() {
    const lowmc::PhotonPhaseSpace phaseSpace = testPhaseSpace();
    static_cast<void>(lowmc::logUniformKTDensityDKT(phaseSpace, 0.5));
}

void sampleLogUniformKTWithInvalidPhaseSpace() {
    lowmc::PhotonPhaseSpace phaseSpace = testPhaseSpace();
    phaseSpace.kTMax = phaseSpace.kTMin;
    static_cast<void>(lowmc::sampleLogUniformKT(phaseSpace, 0.5));
}

void testLogUniformKTRejectsInvalidInputs() {
    requireThrowsInvalidArgument(
        sampleLogUniformKTWithUnitOne,
        "log-uniform kT sampler should reject u=1");
    requireThrowsInvalidArgument(
        densityWithKTOutsideRange,
        "log-uniform kT density should reject kT outside range");
    requireThrowsInvalidArgument(
        sampleLogUniformKTWithInvalidPhaseSpace,
        "log-uniform kT sampler should reject invalid phase space");
}

void testUniformEtaPhiSamplingEndpoints() {
    const lowmc::PhotonPhaseSpace phaseSpace = testPhaseSpace();

    const lowmc::UniformEtaPhiProposalPoint lower =
        lowmc::sampleUniformEtaPhi(phaseSpace, 0.0, 0.0);
    requireNear(
        lower.eta,
        phaseSpace.etaMin,
        1.0e-15,
        "uniform eta-phi sampler should map eta unit 0 to etaMin");
    requireNear(
        lower.phi,
        phaseSpace.phiMin,
        1.0e-15,
        "uniform eta-phi sampler should map phi unit 0 to phiMin");

    const lowmc::UniformEtaPhiProposalPoint middle =
        lowmc::sampleUniformEtaPhi(phaseSpace, 0.5, 0.5);
    requireNear(
        middle.eta,
        0.5 * (phaseSpace.etaMin + phaseSpace.etaMax),
        1.0e-15,
        "uniform eta-phi sampler should map eta unit 0.5 to midpoint");
    requireNear(
        middle.phi,
        0.5 * (phaseSpace.phiMin + phaseSpace.phiMax),
        1.0e-15,
        "uniform eta-phi sampler should map phi unit 0.5 to midpoint");

    const double nearOne = std::nextafter(1.0, 0.0);
    const lowmc::UniformEtaPhiProposalPoint upper =
        lowmc::sampleUniformEtaPhi(phaseSpace, nearOne, nearOne);
    require(
        upper.eta < phaseSpace.etaMax,
        "uniform eta-phi sampler should keep eta below etaMax for u in [0,1)");
    require(
        upper.phi < phaseSpace.phiMax,
        "uniform eta-phi sampler should keep phi below phiMax for u in [0,1)");
}

void testUniformEtaPhiDensityFormula() {
    const lowmc::PhotonPhaseSpace phaseSpace = testPhaseSpace();
    const double area =
        (phaseSpace.etaMax - phaseSpace.etaMin) *
        (phaseSpace.phiMax - phaseSpace.phiMin);

    requireNear(
        lowmc::uniformEtaPhiDensity(phaseSpace),
        1.0 / area,
        1.0e-12,
        "uniform eta-phi density should be inverse accepted area");

    const lowmc::UniformEtaPhiProposalPoint proposal =
        lowmc::sampleUniformEtaPhi(phaseSpace, 0.25, 0.75);
    requireNear(
        proposal.density,
        lowmc::uniformEtaPhiDensity(phaseSpace),
        1.0e-12,
        "sampled eta-phi proposal should include uniform density");
}

void sampleUniformEtaPhiWithUnitOne() {
    static_cast<void>(lowmc::sampleUniformEtaPhi(testPhaseSpace(), 1.0, 0.5));
}

void sampleUniformEtaPhiWithInvalidPhaseSpace() {
    lowmc::PhotonPhaseSpace phaseSpace = testPhaseSpace();
    phaseSpace.phiMax = phaseSpace.phiMin;
    static_cast<void>(lowmc::sampleUniformEtaPhi(phaseSpace, 0.5, 0.5));
}

void testUniformEtaPhiRejectsInvalidInputs() {
    requireThrowsInvalidArgument(
        sampleUniformEtaPhiWithUnitOne,
        "uniform eta-phi sampler should reject unit value 1");
    requireThrowsInvalidArgument(
        sampleUniformEtaPhiWithInvalidPhaseSpace,
        "uniform eta-phi sampler should reject invalid phase space");
}

void testPhotonMomentumFromKTEtaPhiIsLightlikeAndHasRequestedKT() {
    const double kT = 0.08;
    const double eta = 0.7;
    const double phi = 1.2;
    const lowmc::FourVector photon = lowmc::photonMomentumFromKTEtaPhi(kT, eta, phi);

    requireNear(
        photon.pt(),
        kT,
        1.0e-15,
        "photon momentum from kT eta phi should preserve kT");
    requireNear(
        photon.mass2(),
        0.0,
        1.0e-15,
        "photon momentum from kT eta phi should be lightlike");
}

void testFixedKTUniformEtaPhiWeightedSampleUsesSliceProposalDensity() {
    const lowmc::LowCurrent current({
        {lowmc::FourVector{2.0, 0.0, 0.0, 1.0}, +1.0, lowmc::InOutSign::Outgoing},
        {lowmc::FourVector{2.0, 0.0, 0.0, -1.0}, -1.0, lowmc::InOutSign::Outgoing},
    });
    const lowmc::PhotonPhaseSpace phaseSpace = testPhaseSpace();
    const double kT = 0.08;
    const double alpha = 0.01;

    const lowmc::FixedKTWeightedSample sample =
        lowmc::makeFixedKTUniformEtaPhiWeightedSample(
            current,
            phaseSpace,
            kT,
            0.25,
            0.75,
            alpha);

    requireNear(
        sample.kT,
        kT,
        1.0e-15,
        "fixed-kT weighted sample should retain requested kT");
    requireNear(
        sample.etaPhiProposal.density,
        lowmc::uniformEtaPhiDensity(phaseSpace),
        1.0e-12,
        "fixed-kT weighted sample should use conditional eta-phi density");
    requireNear(
        sample.weightedPhoton.proposalDensity,
        sample.etaPhiProposal.density,
        1.0e-12,
        "fixed-kT weighted sample should pass q_eta_phi as proposal density");
    requireNear(
        sample.weightedPhoton.photon.momentum.pt(),
        kT,
        1.0e-15,
        "fixed-kT weighted sample photon should have requested kT");
    requireNear(
        sample.weightedPhoton.photon.weight,
        lowmc::lowPhotonDensity(current, sample.weightedPhoton.photon.momentum, alpha) /
            sample.etaPhiProposal.density,
        1.0e-12,
        "fixed-kT weighted sample should weight by density divided by q_eta_phi");
}

void makeFixedKTSampleOutsideRange() {
    const lowmc::LowCurrent current({
        {lowmc::FourVector{2.0, 0.0, 0.0, 1.0}, +1.0, lowmc::InOutSign::Outgoing},
        {lowmc::FourVector{2.0, 0.0, 0.0, -1.0}, -1.0, lowmc::InOutSign::Outgoing},
    });

    static_cast<void>(lowmc::makeFixedKTUniformEtaPhiWeightedSample(
        current,
        testPhaseSpace(),
        1.0,
        0.5,
        0.5));
}

void testFixedKTUniformEtaPhiWeightedSampleRejectsInvalidInputs() {
    requireThrowsInvalidArgument(
        makeFixedKTSampleOutsideRange,
        "fixed-kT weighted sample should reject kT outside phase-space range");
}

void testAdaptiveEtaPhiProposalSamplesLocalParticlePatch() {
    const double mass = 0.13957039;
    const double momentum = 10.0;
    const double energy = std::sqrt(momentum * momentum + mass * mass);
    const lowmc::PhotonPhaseSpace phaseSpace = testPhaseSpace();
    const lowmc::AdaptiveEtaPhiProposalConfig config{
        0.2,
        5.0,
        0.1,
        0.1,
    };
    const std::vector<lowmc::ChargedParticle> particles = {
        {lowmc::FourVector{energy, 0.0, momentum, 0.0}, +1.0, lowmc::InOutSign::Outgoing},
        {lowmc::FourVector{energy, 0.0, -momentum, 0.0}, -1.0, lowmc::InOutSign::Outgoing},
    };

    const lowmc::AdaptiveEtaPhiProposal proposal =
        lowmc::makeAdaptiveEtaPhiProposal(phaseSpace, particles, config);
    require(proposal.components.size() == 2, "adaptive proposal should create two patches");

    const lowmc::AdaptiveEtaPhiProposalPoint sample =
        lowmc::sampleAdaptiveEtaPhi(proposal, 0.3, 0.5, 0.5);
    requireNear(sample.eta, 0.0, 1.0e-12, "adaptive patch sample eta center");
    requireNear(
        sample.phi,
        0.5 * std::acos(-1.0),
        1.0e-12,
        "adaptive patch sample phi center");

    const double area =
        (phaseSpace.etaMax - phaseSpace.etaMin) *
        (phaseSpace.phiMax - phaseSpace.phiMin);
    const double patchArea = 0.2 * 0.2;
    const double expectedDensity = config.uniformFraction / area + 0.4 / patchArea;
    requireNear(
        sample.density,
        expectedDensity,
        1.0e-12,
        "adaptive proposal density should include uniform floor and selected patch");
}

void testFixedKTAdaptiveEtaPhiWeightedSampleUsesAdaptiveDensity() {
    const double mass = 0.13957039;
    const double momentum = 10.0;
    const double energy = std::sqrt(momentum * momentum + mass * mass);
    const lowmc::LowCurrent current({
        {lowmc::FourVector{energy, 0.0, momentum, 0.0}, +1.0, lowmc::InOutSign::Outgoing},
        {lowmc::FourVector{energy, 0.0, -momentum, 0.0}, -1.0, lowmc::InOutSign::Outgoing},
    });
    const lowmc::PhotonPhaseSpace phaseSpace = testPhaseSpace();
    const lowmc::AdaptiveEtaPhiProposalConfig config{
        0.2,
        5.0,
        0.1,
        0.1,
    };
    const double kT = 0.08;
    const double alpha = 0.01;

    const lowmc::FixedKTAdaptiveWeightedSample sample =
        lowmc::makeFixedKTAdaptiveEtaPhiWeightedSample(
            current,
            phaseSpace,
            kT,
            0.3,
            0.5,
            0.5,
            config,
            alpha);

    requireNear(
        sample.weightedPhoton.proposalDensity,
        sample.etaPhiProposal.density,
        1.0e-12,
        "adaptive fixed-kT sample should pass adaptive q_eta_phi");
    requireNear(
        sample.weightedPhoton.photon.weight,
        lowmc::lowPhotonDensity(current, sample.weightedPhoton.photon.momentum, alpha) /
            sample.etaPhiProposal.density,
        1.0e-12,
        "adaptive fixed-kT sample should weight by exact density divided by adaptive q");
}

void testSmoothEtaPhiProposalSamplesLocalKernel() {
    const double mass = 0.13957039;
    const double momentum = 10.0;
    const double energy = std::sqrt(momentum * momentum + mass * mass);
    const lowmc::PhotonPhaseSpace phaseSpace = testPhaseSpace();
    const lowmc::SmoothEtaPhiProposalConfig config{
        0.2,
        5.0,
        0.1,
        0.1,
    };
    const std::vector<lowmc::ChargedParticle> particles = {
        {lowmc::FourVector{energy, 0.0, momentum, 0.0}, +1.0, lowmc::InOutSign::Outgoing},
        {lowmc::FourVector{energy, 0.0, -momentum, 0.0}, -1.0, lowmc::InOutSign::Outgoing},
    };

    const lowmc::SmoothEtaPhiProposal proposal =
        lowmc::makeSmoothEtaPhiProposal(phaseSpace, particles, config);
    require(proposal.components.size() == 2, "smooth proposal should create two kernels");

    const lowmc::SmoothEtaPhiProposalPoint sample =
        lowmc::sampleSmoothEtaPhi(proposal, 0.3, 0.5, 0.5);
    requireNear(sample.eta, 0.0, 1.0e-12, "smooth kernel sample eta center");
    requireNear(
        sample.phi,
        0.5 * std::acos(-1.0),
        1.0e-12,
        "smooth kernel sample phi center");

    const double area =
        (phaseSpace.etaMax - phaseSpace.etaMin) *
        (phaseSpace.phiMax - phaseSpace.phiMin);
    const double expectedDensity = config.uniformFraction / area + 0.4 * 100.0;
    requireNear(
        sample.density,
        expectedDensity,
        1.0e-12,
        "smooth proposal density should include uniform floor and selected kernel");
}

void testFixedKTSmoothEtaPhiWeightedSampleUsesSmoothDensity() {
    const double mass = 0.13957039;
    const double momentum = 10.0;
    const double energy = std::sqrt(momentum * momentum + mass * mass);
    const lowmc::LowCurrent current({
        {lowmc::FourVector{energy, 0.0, momentum, 0.0}, +1.0, lowmc::InOutSign::Outgoing},
        {lowmc::FourVector{energy, 0.0, -momentum, 0.0}, -1.0, lowmc::InOutSign::Outgoing},
    });
    const lowmc::PhotonPhaseSpace phaseSpace = testPhaseSpace();
    const lowmc::SmoothEtaPhiProposalConfig config{
        0.2,
        5.0,
        0.1,
        0.1,
    };
    const double kT = 0.08;
    const double alpha = 0.01;

    const lowmc::FixedKTSmoothWeightedSample sample =
        lowmc::makeFixedKTSmoothEtaPhiWeightedSample(
            current,
            phaseSpace,
            kT,
            0.3,
            0.5,
            0.5,
            config,
            alpha);

    requireNear(
        sample.weightedPhoton.proposalDensity,
        sample.etaPhiProposal.density,
        1.0e-12,
        "smooth fixed-kT sample should pass smooth q_eta_phi");
    requireNear(
        sample.weightedPhoton.photon.weight,
        lowmc::lowPhotonDensity(current, sample.weightedPhoton.photon.momentum, alpha) /
            sample.etaPhiProposal.density,
        1.0e-12,
        "smooth fixed-kT sample should weight by exact density divided by smooth q");
}

void testGridAdaptiveEtaPhiProposalBuildsRefinedCells() {
    const lowmc::LowCurrent current({
        {lowmc::FourVector{2.0, 0.0, 0.0, 1.0}, +1.0, lowmc::InOutSign::Outgoing},
        {lowmc::FourVector{2.0, 0.0, 0.0, -1.0}, -1.0, lowmc::InOutSign::Outgoing},
    });
    const lowmc::PhotonPhaseSpace phaseSpace = testPhaseSpace();
    const lowmc::GridAdaptiveEtaPhiProposalConfig config{
        2,
        2,
        1,
        0.25,
        0.2,
    };

    const lowmc::GridAdaptiveEtaPhiProposal proposal =
        lowmc::makeGridAdaptiveEtaPhiProposal(current, phaseSpace, 0.08, config);
    require(
        proposal.cells.size() == 7,
        "grid-adaptive proposal should replace one refined cell by four children");

    double mixtureWeight = 0.0;
    for (const lowmc::GridAdaptiveEtaPhiCell& cell : proposal.cells) {
        require(cell.targetMass >= 0.0, "grid-adaptive target masses should be non-negative");
        mixtureWeight += cell.mixtureWeight;
    }
    requireNear(
        mixtureWeight,
        1.0 - config.uniformFraction,
        1.0e-12,
        "grid-adaptive cell mixture weights should sum to adaptive fraction");
}

void testGridAdaptiveEtaPhiSampleUsesPositiveDensity() {
    const lowmc::LowCurrent current({
        {lowmc::FourVector{2.0, 0.0, 0.0, 1.0}, +1.0, lowmc::InOutSign::Outgoing},
        {lowmc::FourVector{2.0, 0.0, 0.0, -1.0}, -1.0, lowmc::InOutSign::Outgoing},
    });
    const lowmc::GridAdaptiveEtaPhiProposalConfig config{
        4,
        4,
        1,
        0.25,
        0.2,
    };
    const lowmc::GridAdaptiveEtaPhiProposal proposal =
        lowmc::makeGridAdaptiveEtaPhiProposal(current, testPhaseSpace(), 0.08, config);
    const lowmc::GridAdaptiveEtaPhiProposalPoint sample =
        lowmc::sampleGridAdaptiveEtaPhi(proposal, 0.5, 0.25, 0.5, 0.5);

    require(std::isfinite(sample.density), "grid-adaptive sample density should be finite");
    require(sample.density > 0.0, "grid-adaptive sample density should be positive");
    requireNear(
        sample.density,
        lowmc::gridAdaptiveEtaPhiDensity(proposal, sample.eta, sample.phi),
        1.0e-12,
        "grid-adaptive sampled density should match density evaluator");
}

void testGridAdaptiveVarianceRefinementScoreBuildsProposal() {
    const lowmc::LowCurrent current({
        {lowmc::FourVector{2.0, 0.0, 0.0, 1.0}, +1.0, lowmc::InOutSign::Outgoing},
        {lowmc::FourVector{2.0, 0.0, 0.0, -1.0}, -1.0, lowmc::InOutSign::Outgoing},
    });
    lowmc::GridAdaptiveEtaPhiProposalConfig config;
    config.etaBins = 4;
    config.phiBins = 4;
    config.refinementLevels = 1;
    config.refineFraction = 0.25;
    config.uniformFraction = 0.2;
    config.refinementScore = lowmc::GridAdaptiveRefinementScore::Variance;
    config.varianceProbeBins = 3;

    const lowmc::GridAdaptiveEtaPhiProposal proposal =
        lowmc::makeGridAdaptiveEtaPhiProposal(current, testPhaseSpace(), 0.08, config);

    require(!proposal.cells.empty(), "variance-refined proposal should build cells");
    for (const lowmc::GridAdaptiveEtaPhiCell& cell : proposal.cells) {
        require(
            std::isfinite(cell.refinementScore),
            "variance-refined cell score should be finite");
        require(
            cell.refinementScore >= 0.0,
            "variance-refined cell score should be non-negative");
    }
}

void testFixedKTGridAdaptiveEtaPhiWeightedSampleUsesGridDensity() {
    const lowmc::LowCurrent current({
        {lowmc::FourVector{2.0, 0.0, 0.0, 1.0}, +1.0, lowmc::InOutSign::Outgoing},
        {lowmc::FourVector{2.0, 0.0, 0.0, -1.0}, -1.0, lowmc::InOutSign::Outgoing},
    });
    const lowmc::GridAdaptiveEtaPhiProposalConfig config{
        4,
        4,
        1,
        0.25,
        0.2,
    };
    const double kT = 0.08;
    const double alpha = 0.01;

    const lowmc::FixedKTGridAdaptiveWeightedSample sample =
        lowmc::makeFixedKTGridAdaptiveEtaPhiWeightedSample(
            current,
            testPhaseSpace(),
            kT,
            0.5,
            0.25,
            0.5,
            0.5,
            config,
            alpha);

    requireNear(
        sample.weightedPhoton.proposalDensity,
        sample.etaPhiProposal.density,
        1.0e-12,
        "grid-adaptive fixed-kT sample should pass grid q_eta_phi");
    requireNear(
        sample.weightedPhoton.photon.weight,
        lowmc::lowPhotonDensity(current, sample.weightedPhoton.photon.momentum, alpha) /
            sample.etaPhiProposal.density,
        1.0e-12,
        "grid-adaptive fixed-kT sample should weight by exact density divided by grid q");
}

void testMinimalWeightedPhotonSampleUsesFullProposalDensity() {
    const lowmc::LowCurrent current({
        {lowmc::FourVector{2.0, 0.0, 0.0, 1.0}, +1.0, lowmc::InOutSign::Outgoing},
        {lowmc::FourVector{2.0, 0.0, 0.0, -1.0}, -1.0, lowmc::InOutSign::Outgoing},
    });
    const lowmc::PhotonPhaseSpace phaseSpace = testPhaseSpace();
    const double alpha = 0.01;

    const lowmc::MinimalWeightedPhotonSample sample = lowmc::makeMinimalWeightedPhotonSample(
        current,
        phaseSpace,
        0.5,
        0.25,
        0.75,
        alpha);

    const double expectedProposalDensity =
        sample.kTProposal.densityPhysicalMeasure * sample.etaPhiProposal.density;
    requireNear(
        sample.proposalDensity,
        expectedProposalDensity,
        1.0e-12,
        "minimal weighted sample should combine q_kT_mu and q_eta_phi");
    requireNear(
        sample.weightedPhoton.proposalDensity,
        expectedProposalDensity,
        1.0e-12,
        "minimal weighted photon should receive full q_mu");
    requireNear(
        sample.weightedPhoton.photon.momentum.pt(),
        sample.kTProposal.kT,
        1.0e-15,
        "minimal weighted sample photon should use sampled kT");
    requireNear(
        sample.weightedPhoton.photon.weight,
        lowmc::lowPhotonDensity(current, sample.weightedPhoton.photon.momentum, alpha) /
            expectedProposalDensity,
        1.0e-12,
        "minimal weighted sample should weight by density divided by full q_mu");
}

void testGridAdaptiveWeightedPhotonSampleUsesFullProposalDensity() {
    const lowmc::LowCurrent current({
        {lowmc::FourVector{2.0, 0.0, 0.0, 1.0}, +1.0, lowmc::InOutSign::Outgoing},
        {lowmc::FourVector{2.0, 0.0, 0.0, -1.0}, -1.0, lowmc::InOutSign::Outgoing},
    });
    const lowmc::PhotonPhaseSpace phaseSpace = testPhaseSpace();
    const lowmc::GridAdaptiveEtaPhiProposalConfig config{
        4,
        4,
        1,
        0.25,
        0.2,
    };
    const lowmc::GridAdaptiveEtaPhiProposal proposal =
        lowmc::makeGridAdaptiveEtaPhiProposal(current, phaseSpace, 0.08, config);
    const double alpha = 0.01;

    const lowmc::GridAdaptiveWeightedPhotonSample sample =
        lowmc::makeGridAdaptiveWeightedPhotonSample(
            current,
            phaseSpace,
            proposal,
            0.5,
            0.5,
            0.25,
            0.5,
            0.5,
            alpha);

    const double expectedProposalDensity =
        sample.kTProposal.densityPhysicalMeasure * sample.etaPhiProposal.density;
    requireNear(
        sample.proposalDensity,
        expectedProposalDensity,
        1.0e-12,
        "grid full weighted sample should combine q_kT_mu and q_eta_phi");
    requireNear(
        sample.weightedPhoton.proposalDensity,
        expectedProposalDensity,
        1.0e-12,
        "grid full weighted photon should receive full q_mu");
    requireNear(
        sample.weightedPhoton.photon.weight,
        lowmc::lowPhotonDensity(current, sample.weightedPhoton.photon.momentum, alpha) /
            expectedProposalDensity,
        1.0e-12,
        "grid full weighted sample should weight by density divided by full q_mu");
}

void makeMinimalWeightedSampleWithInvalidUnitKT() {
    const lowmc::LowCurrent current({
        {lowmc::FourVector{2.0, 0.0, 0.0, 1.0}, +1.0, lowmc::InOutSign::Outgoing},
        {lowmc::FourVector{2.0, 0.0, 0.0, -1.0}, -1.0, lowmc::InOutSign::Outgoing},
    });

    static_cast<void>(lowmc::makeMinimalWeightedPhotonSample(
        current,
        testPhaseSpace(),
        1.0,
        0.5,
        0.5));
}

void testMinimalWeightedPhotonSampleRejectsInvalidInputs() {
    requireThrowsInvalidArgument(
        makeMinimalWeightedSampleWithInvalidUnitKT,
        "minimal weighted sample should reject invalid kT unit value");
}

void testWeightSummaryComputesBasicDiagnostics() {
    const lowmc::WeightSummary summary = lowmc::summarizeWeights({1.0, 2.0, 3.0});

    require(summary.count == 3, "weight summary should report count");
    requireNear(summary.sumWeight, 6.0, 1.0e-12, "weight summary sum");
    requireNear(summary.minWeight, 1.0, 1.0e-12, "weight summary minimum");
    requireNear(summary.maxWeight, 3.0, 1.0e-12, "weight summary maximum");
    requireNear(summary.meanWeight, 2.0, 1.0e-12, "weight summary mean");
    requireNear(
        summary.meanSquaredWeight,
        14.0 / 3.0,
        1.0e-12,
        "weight summary mean squared weight");
    requireNear(
        summary.varianceWeight,
        2.0 / 3.0,
        1.0e-12,
        "weight summary variance");
    requireNear(
        summary.effectiveSampleSize,
        36.0 / 14.0,
        1.0e-12,
        "weight summary effective sample size");
}

void testWeightSummaryHandlesEmptyAndZeroWeights() {
    const lowmc::WeightSummary empty = lowmc::summarizeWeights({});
    require(empty.count == 0, "empty weight summary should report zero count");
    requireNear(empty.effectiveSampleSize, 0.0, 1.0e-12, "empty effective sample size");

    const lowmc::WeightSummary zeros = lowmc::summarizeWeights({0.0, 0.0});
    require(zeros.count == 2, "zero weight summary should preserve count");
    requireNear(zeros.sumWeight, 0.0, 1.0e-12, "zero weight summary sum");
    requireNear(zeros.effectiveSampleSize, 0.0, 1.0e-12, "zero effective sample size");
}

void summarizeNegativeWeight() {
    static_cast<void>(lowmc::summarizeWeights({1.0, -1.0}));
}

void testWeightSummaryRejectsInvalidWeights() {
    requireThrowsInvalidArgument(
        summarizeNegativeWeight,
        "weight summary should reject negative weights");
}

void testWeightSummaryCanSummarizePhotonsAndWeightedSamples() {
    const std::vector<lowmc::Photon> photons = {
        {lowmc::FourVector{}, 1.0},
        {lowmc::FourVector{}, 3.0},
    };
    const lowmc::WeightSummary photonSummary = lowmc::summarizePhotons(photons);
    requireNear(photonSummary.sumWeight, 4.0, 1.0e-12, "photon weight summary sum");

    const std::vector<lowmc::WeightedPhotonSample> samples = {
        {lowmc::Photon{lowmc::FourVector{}, 2.0}, 0.5},
        {lowmc::Photon{lowmc::FourVector{}, 4.0}, 0.25},
    };
    const lowmc::WeightSummary sampleSummary = lowmc::summarizeWeightedPhotonSamples(samples);
    requireNear(sampleSummary.meanWeight, 3.0, 1.0e-12, "sample weight summary mean");
}

std::vector<lowmc::ChargedParticle> eventGeneratorToyParticles() {
    const double mass = 0.13957039;
    const double momentum = 10.0;
    const double energy = std::sqrt(momentum * momentum + mass * mass);
    return {
        {lowmc::FourVector{energy, 0.0, momentum, 0.0}, +1.0, lowmc::InOutSign::Outgoing},
        {lowmc::FourVector{energy, 0.0, -momentum, 0.0}, -1.0, lowmc::InOutSign::Outgoing},
    };
}

lowmc::WeightedEventGeneratorConfig eventGeneratorTestConfig(std::uint64_t seed) {
    lowmc::GridAdaptiveEtaPhiProposalConfig gridConfig;
    gridConfig.etaBins = 4;
    gridConfig.phiBins = 4;
    gridConfig.refinementLevels = 1;
    gridConfig.refineFraction = 0.25;
    gridConfig.uniformFraction = 0.2;
    return lowmc::WeightedEventGeneratorConfig{
        testPhaseSpace(),
        6,
        seed,
        gridConfig,
    };
}

lowmc::ShapeUnweightingConfig shapeUnweightingTestConfig(std::uint64_t seed) {
    lowmc::GridAdaptiveEtaPhiProposalConfig gridConfig;
    gridConfig.etaBins = 4;
    gridConfig.phiBins = 4;
    gridConfig.refinementLevels = 1;
    gridConfig.refineFraction = 0.25;
    gridConfig.uniformFraction = 0.2;
    return lowmc::ShapeUnweightingConfig{
        testPhaseSpace(),
        12,
        seed,
        gridConfig,
        12,
        2,
        lowmc::ShapeUnweightingProbeRule::Interior,
        lowmc::ShapeUnweightingProbePolicy::Fixed,
        4,
        100.0,
        1.0,
        0,
        0.0,
        100.0,
        lowmc::EnvelopeFailurePolicy::Strict,
    };
}

lowmc::EventYieldIntegrationConfig eventYieldTestIntegrationConfig() {
    return lowmc::EventYieldIntegrationConfig{
        4,
        8,
        10,
    };
}

void testEventYieldEstimateIsPositiveAndDeterministic() {
    const lowmc::EventYieldEstimate first =
        lowmc::estimateEventYield(
            eventGeneratorToyParticles(),
            testPhaseSpace(),
            eventYieldTestIntegrationConfig());
    const lowmc::EventYieldEstimate second =
        lowmc::estimateEventYield(
            eventGeneratorToyParticles(),
            testPhaseSpace(),
            eventYieldTestIntegrationConfig());

    require(first.lambdaEvent > 0.0, "event yield estimate should be positive");
    require(std::isfinite(first.lambdaEvent), "event yield estimate should be finite");
    require(first.coarseLambdaEvent > 0.0, "event yield coarse estimate should be positive");
    require(
        std::string(first.method) == "log-kt-eta-phi-midpoint",
        "event yield estimate should report integration method");
    require(std::string(first.status) == "ok", "event yield estimate should report ok status");
    require(
        first.evaluations == first.logKTBins * first.etaBins * first.phiBins,
        "event yield estimate should report evaluation count");
    requireNear(
        first.lambdaEvent,
        second.lambdaEvent,
        1.0e-15,
        "event yield estimate should be deterministic");
}

void testEventYieldEstimateUsesLogKTMeasure() {
    lowmc::PhotonPhaseSpace shiftedPhaseSpace = testPhaseSpace();
    shiftedPhaseSpace.kTMin *= 10.0;
    shiftedPhaseSpace.kTMax *= 10.0;

    const lowmc::EventYieldIntegrationConfig config = eventYieldTestIntegrationConfig();
    const lowmc::EventYieldEstimate base =
        lowmc::estimateEventYield(eventGeneratorToyParticles(), testPhaseSpace(), config);
    const lowmc::EventYieldEstimate shifted =
        lowmc::estimateEventYield(eventGeneratorToyParticles(), shiftedPhaseSpace, config);

    requireNear(
        shifted.lambdaEvent,
        base.lambdaEvent,
        1.0e-12,
        "event yield should be invariant when the same log-kT interval is rescaled");
}

void estimateEventYieldWithInvalidGrid() {
    lowmc::EventYieldIntegrationConfig config = eventYieldTestIntegrationConfig();
    config.etaBins = 0;
    static_cast<void>(
        lowmc::estimateEventYield(eventGeneratorToyParticles(), testPhaseSpace(), config));
}

void estimateEventYieldWithInvalidPhaseSpace() {
    lowmc::PhotonPhaseSpace phaseSpace = testPhaseSpace();
    phaseSpace.kTMax = phaseSpace.kTMin;
    static_cast<void>(
        lowmc::estimateEventYield(
            eventGeneratorToyParticles(),
            phaseSpace,
            eventYieldTestIntegrationConfig()));
}

void testEventYieldEstimateRejectsInvalidInputs() {
    requireThrowsInvalidArgument(
        estimateEventYieldWithInvalidGrid,
        "event yield estimate should reject invalid grid");
    requireThrowsInvalidArgument(
        estimateEventYieldWithInvalidPhaseSpace,
        "event yield estimate should reject invalid phase space");
}

void testAdaptiveEventYieldEstimateIsPositiveAndReportsCells() {
    lowmc::AdaptiveEventYieldIntegrationConfig config;
    config.logKTBins = 2;
    config.etaSubBins = 1;
    config.phiSubBins = 1;
    config.gridConfig.etaBins = 4;
    config.gridConfig.phiBins = 4;
    config.gridConfig.refinementLevels = 1;
    config.gridConfig.refineFraction = 0.25;

    const lowmc::EventYieldEstimate estimate =
        lowmc::estimateAdaptiveEventYield(eventGeneratorToyParticles(), testPhaseSpace(), config);

    require(estimate.lambdaEvent > 0.0, "adaptive event yield should be positive");
    require(std::isfinite(estimate.lambdaEvent), "adaptive event yield should be finite");
    require(
        std::string(estimate.method) == "adaptive-log-kt-grid-eta-phi-midpoint",
        "adaptive event yield should report method name");
    require(estimate.gridCells > 0, "adaptive event yield should report grid cells");
    require(
        estimate.evaluations == estimate.gridCells,
        "one sub-bin adaptive yield should evaluate once per reported grid cell");
}

void estimateAdaptiveEventYieldWithInvalidConfig() {
    lowmc::AdaptiveEventYieldIntegrationConfig config;
    config.logKTBins = 0;
    static_cast<void>(
        lowmc::estimateAdaptiveEventYield(eventGeneratorToyParticles(), testPhaseSpace(), config));
}

void testAdaptiveEventYieldRejectsInvalidInputs() {
    requireThrowsInvalidArgument(
        estimateAdaptiveEventYieldWithInvalidConfig,
        "adaptive event yield should reject invalid config");
}

void testPhotonMultiplicityUsesPoissonMeanAndOversamplingWeight() {
    const double lambdaEvent = 0.25;
    const lowmc::PhotonMultiplicityConfig config{
        4.0,
        12345,
    };
    const lowmc::PhotonMultiplicityDiagnostic first =
        lowmc::samplePhotonMultiplicity(lambdaEvent, config);
    const lowmc::PhotonMultiplicityDiagnostic second =
        lowmc::samplePhotonMultiplicity(lambdaEvent, config);

    requireNear(first.lambdaEvent, lambdaEvent, 1.0e-15, "multiplicity lambda_event");
    requireNear(first.poissonMean, 1.0, 1.0e-15, "oversampled Poisson mean");
    requireNear(first.sampleWeight, 0.25, 1.0e-15, "oversampled analysis weight");
    requireNear(
        first.expectedWeightedCount,
        first.sampleWeight * static_cast<double>(first.sampledPhotons),
        1.0e-15,
        "weighted count should use diagnostic sample weight");
    require(
        first.sampledPhotons == second.sampledPhotons,
        "same-seed multiplicity diagnostic should be reproducible");
}

void samplePhotonMultiplicityWithNegativeLambda() {
    static_cast<void>(lowmc::samplePhotonMultiplicity(-1.0));
}

void samplePhotonMultiplicityWithInvalidOversampling() {
    static_cast<void>(
        lowmc::samplePhotonMultiplicity(1.0, lowmc::PhotonMultiplicityConfig{0.0, 12345}));
}

void testPhotonMultiplicityRejectsInvalidInputs() {
    requireThrowsInvalidArgument(
        samplePhotonMultiplicityWithNegativeLambda,
        "multiplicity diagnostic should reject negative lambda_event");
    requireThrowsInvalidArgument(
        samplePhotonMultiplicityWithInvalidOversampling,
        "multiplicity diagnostic should reject invalid oversampling factor");
}

void testParticleSelectionConfigDefaultsMatchAdapterContract() {
    const lowmc::ParticleSelectionConfig config;

    require(
        config.enabled(lowmc::ParticleSelection::IncludeIncomingBeams),
        "adapter selection should include incoming beams by default");
    require(
        config.enabled(lowmc::ParticleSelection::IncludeFinalStableCharged),
        "adapter selection should include final stable charged particles by default");
    require(
        !config.enabled(lowmc::ParticleSelection::IncludeUnstableIntermediateCharged),
        "adapter selection should exclude unstable intermediate charged particles by default");
    require(
        !config.enabled(lowmc::ParticleSelection::IncludeBeamRemnants),
        "adapter selection should exclude beam remnants by default");
}

void testParticleSelectionConfigCanToggleOptions() {
    lowmc::ParticleSelectionConfig config;

    config.set(lowmc::ParticleSelection::IncludeIncomingBeams, false);
    config.set(lowmc::ParticleSelection::IncludeBeamRemnants, true);

    require(
        !config.enabled(lowmc::ParticleSelection::IncludeIncomingBeams),
        "adapter selection should allow disabling incoming beams");
    require(
        config.enabled(lowmc::ParticleSelection::IncludeBeamRemnants),
        "adapter selection should allow explicitly enabling beam remnants");
}

void testAdapterSelectionDiagnosticSummarizesSelectedParticles() {
    const std::vector<lowmc::ChargedParticle> selectedParticles = {
        {lowmc::FourVector{10.0, 0.0, 0.0, 10.0}, +1.0, lowmc::InOutSign::Incoming},
        {lowmc::FourVector{10.0, 0.0, 0.0, -10.0}, -1.0, lowmc::InOutSign::Incoming},
        {lowmc::FourVector{2.0, 0.2, 0.0, 1.0}, +1.0, lowmc::InOutSign::Outgoing},
        {lowmc::FourVector{3.0, -0.2, 0.0, -1.0}, -1.0, lowmc::InOutSign::Outgoing},
    };
    std::vector<lowmc::RejectedChargedCandidate> rejected = {
        {lowmc::RejectedParticleReason::UnstableIntermediate, +1.0, "generator-index=42"},
        {lowmc::RejectedParticleReason::BeamRemnant, -1.0, "generator-index=43"},
    };

    const lowmc::AdapterSelectionDiagnostic diagnostic =
        lowmc::summarizeAdapterSelection(selectedParticles, rejected);

    require(
        diagnostic.selectedIncomingCharged == 2,
        "adapter diagnostic should count selected incoming charged particles");
    require(
        diagnostic.selectedOutgoingCharged == 2,
        "adapter diagnostic should count selected outgoing charged particles");
    require(
        diagnostic.selectedChargedCount() == selectedParticles.size(),
        "adapter diagnostic should report total selected charged count");
    require(
        diagnostic.rejectedChargedCount() == 2,
        "adapter diagnostic should report rejected charged candidate count");
    requireNear(
        diagnostic.signedChargeSum,
        0.0,
        1.0e-12,
        "adapter diagnostic should compute signed charge sum");
    requireNear(
        diagnostic.ordinaryChargeSum,
        0.0,
        1.0e-12,
        "adapter diagnostic should compute ordinary charge sum");
    require(
        diagnostic.rejectedChargedCandidates[0].reason ==
            lowmc::RejectedParticleReason::UnstableIntermediate,
        "adapter diagnostic should retain rejected charged candidate reasons");
    require(
        diagnostic.rejectedChargedCandidates[0].sourceLabel == "generator-index=42",
        "adapter diagnostic should retain rejected charged candidate source labels");
}

void testAdapterParticleSelectionUsesPDGChargeAndDefaultPolicy() {
    const std::vector<lowmc::AdapterParticleCandidate> candidates = {
        {lowmc::FourVector{100.0, 0.0, 0.0, 100.0}, 2212, true, false, false, "beam"},
        {lowmc::FourVector{2.0, 1.0, 0.0, 0.0}, -211, false, true, false, "pi-"},
        {lowmc::FourVector{1.0, 0.0, 0.0, 1.0}, 22, false, true, false, "gamma"},
        {lowmc::FourVector{3.0, 0.0, 1.0, 0.0}, 9'000'000, false, true, false, "private"},
        {lowmc::FourVector{4.0, 0.0, -1.0, 0.0}, 211, false, true, true, "remnant"},
        {lowmc::FourVector{5.0, 0.0, 0.0, 2.0}, 321, false, false, false, "intermediate"},
    };

    const lowmc::AdapterParticleSelectionResult result =
        lowmc::selectAdapterParticles(candidates);

    require(
        result.selectedParticles.size() == 2,
        "adapter particle selection should select incoming and final stable charged records");
    require(
        result.selectedParticles[0].inOutSign == lowmc::InOutSign::Incoming,
        "adapter particle selection should mark incoming beams");
    requireNear(
        result.selectedParticles[0].charge,
        +1.0,
        1.0e-12,
        "adapter particle selection should assign proton charge");
    require(
        result.selectedParticles[1].inOutSign == lowmc::InOutSign::Outgoing,
        "adapter particle selection should mark final stable records as outgoing");
    requireNear(
        result.selectedParticles[1].charge,
        -1.0,
        1.0e-12,
        "adapter particle selection should assign pion charge");

    require(
        result.diagnostic.selectedIncomingCharged == 1,
        "adapter particle selection diagnostic should count selected incoming charged records");
    require(
        result.diagnostic.selectedOutgoingCharged == 1,
        "adapter particle selection diagnostic should count selected outgoing charged records");
    require(
        result.diagnostic.ignoredNeutralCandidates == 1,
        "adapter particle selection diagnostic should count ignored neutral candidates");
    require(
        result.diagnostic.rejectedChargedCount() == 3,
        "adapter particle selection diagnostic should count rejected charged candidates");
    requireNear(
        result.diagnostic.signedChargeSum,
        -2.0,
        1.0e-12,
        "adapter particle selection diagnostic should compute signed charge sum");
    requireNear(
        result.diagnostic.ordinaryChargeSum,
        0.0,
        1.0e-12,
        "adapter particle selection diagnostic should compute ordinary charge sum");
    require(
        result.diagnostic.rejectedChargedCandidates[0].reason ==
            lowmc::RejectedParticleReason::UnsupportedCharge,
        "adapter particle selection should reject unsupported PDG ids explicitly");
    require(
        result.diagnostic.rejectedChargedCandidates[0].sourceLabel == "private",
        "adapter particle selection should keep source label for unsupported charge");
    require(
        result.diagnostic.rejectedChargedCandidates[1].reason ==
            lowmc::RejectedParticleReason::BeamRemnant,
        "adapter particle selection should reject beam remnants by default");
    require(
        result.diagnostic.rejectedChargedCandidates[2].reason ==
            lowmc::RejectedParticleReason::UnstableIntermediate,
        "adapter particle selection should reject unstable intermediates by default");
}

void testAdapterParticleSelectionCanDisableIncomingBeams() {
    lowmc::ParticleSelectionConfig config;
    config.includeIncomingBeams = false;

    const lowmc::AdapterParticleSelectionResult result = lowmc::selectAdapterParticles(
        {{lowmc::FourVector{100.0, 0.0, 0.0, 100.0}, 2212, true, false, false, "beam"}},
        config);

    require(
        result.selectedParticles.empty(),
        "adapter particle selection should obey disabled incoming-beam option");
    require(
        result.diagnostic.rejectedChargedCount() == 1,
        "adapter particle selection should report excluded incoming beam");
    require(
        result.diagnostic.rejectedChargedCandidates[0].reason ==
            lowmc::RejectedParticleReason::OutsideSelection,
        "adapter particle selection should use outside-selection reason for disabled incoming");
}

void testAdapterParticleSelectionCanIncludeUnstableIntermediates() {
    lowmc::ParticleSelectionConfig config;
    config.includeUnstableIntermediateCharged = true;

    const lowmc::AdapterParticleSelectionResult result = lowmc::selectAdapterParticles(
        {{lowmc::FourVector{5.0, 0.0, 0.0, 2.0}, 321, false, false, false, "intermediate"}},
        config);

    require(
        result.selectedParticles.size() == 1,
        "adapter particle selection should allow explicit unstable-intermediate inclusion");
    require(
        result.selectedParticles[0].inOutSign == lowmc::InOutSign::Outgoing,
        "included unstable intermediates should be treated as outgoing candidates");
    requireNear(
        result.selectedParticles[0].charge,
        +1.0,
        1.0e-12,
        "included unstable intermediates should use PDG charge");
}

void testAdapterParticleSelectionRejectsAmbiguousRecords() {
    const lowmc::AdapterParticleSelectionResult result = lowmc::selectAdapterParticles({
        {lowmc::FourVector{10.0, 0.0, 0.0, 10.0}, 2212, true, true, false, "ambiguous"},
    });

    require(
        result.selectedParticles.empty(),
        "adapter particle selection should not select ambiguous records");
    require(
        result.diagnostic.rejectedChargedCount() == 1,
        "adapter particle selection should report ambiguous charged records");
    require(
        result.diagnostic.rejectedChargedCandidates[0].reason ==
            lowmc::RejectedParticleReason::AmbiguousRecord,
        "adapter particle selection should use ambiguous-record reason");
}

void requirePDGCharge(
    int pdgId,
    int expectedChargeTimesThree,
    const std::string& message) {
    const lowmc::PDGChargeLookupResult result = lowmc::chargeTimesThreeFromPDG(pdgId);
    require(result.known, message + " should have known charge");
    require(
        result.chargeTimesThree == expectedChargeTimesThree,
        message + " should have expected charge");
}

void requireUnknownPDGCharge(int pdgId, const std::string& message) {
    const lowmc::PDGChargeLookupResult result = lowmc::chargeTimesThreeFromPDG(pdgId);
    require(!result.known, message + " should have unsupported charge");
}

void testPDGChargeLookupHandlesElementaryParticles() {
    requirePDGCharge(2, +2, "u quark");
    requirePDGCharge(-2, -2, "anti-u quark");
    requirePDGCharge(1, -1, "d quark");
    requirePDGCharge(-1, +1, "anti-d quark");
    requirePDGCharge(11, -3, "electron");
    requirePDGCharge(-11, +3, "positron");
    requirePDGCharge(12, 0, "electron neutrino");
    requirePDGCharge(22, 0, "photon");
    requirePDGCharge(21, 0, "gluon");
    requirePDGCharge(23, 0, "Z boson");
    requirePDGCharge(24, +3, "W plus");
    requirePDGCharge(-24, -3, "W minus");
    requirePDGCharge(25, 0, "Higgs boson");
}

void testPDGChargeLookupHandlesOrdinaryMesons() {
    requirePDGCharge(211, +3, "pi plus");
    requirePDGCharge(-211, -3, "pi minus");
    requirePDGCharge(111, 0, "pi zero");
    requirePDGCharge(321, +3, "K plus");
    requirePDGCharge(-321, -3, "K minus");
    requirePDGCharge(311, 0, "K zero");
    requirePDGCharge(411, +3, "D plus");
    requirePDGCharge(421, 0, "D zero");
    requirePDGCharge(521, +3, "B plus");
    requirePDGCharge(511, 0, "B zero");
}

void testPDGChargeLookupHandlesOrdinaryBaryons() {
    requirePDGCharge(2212, +3, "proton");
    requirePDGCharge(-2212, -3, "antiproton");
    requirePDGCharge(2112, 0, "neutron");
    requirePDGCharge(3122, 0, "Lambda zero");
    requirePDGCharge(3222, +3, "Sigma plus");
    requirePDGCharge(3112, -3, "Sigma minus");
    requirePDGCharge(3312, -3, "Xi minus");
    requirePDGCharge(3334, -3, "Omega minus");
    requirePDGCharge(-3334, +3, "anti-Omega plus");
}

void testPDGChargeLookupHandlesIonsAndUnknownIds() {
    requirePDGCharge(1'000'010'020, +3, "deuteron");
    requirePDGCharge(-1'000'010'020, -3, "anti-deuteron");
    requirePDGCharge(1'000'020'040, +6, "alpha particle");
    requireUnknownPDGCharge(0, "zero PDG id");
    requireUnknownPDGCharge(9'000'000, "generator-specific private id");
    requireUnknownPDGCharge(10'000'001, "unsupported extension id");
    requireUnknownPDGCharge(1103, "diquark id");
}

void testWeightedEventGeneratorIsReproducible() {
    const lowmc::WeightedEventGenerationResult first =
        lowmc::generateWeightedPhotons(
            eventGeneratorToyParticles(),
            eventGeneratorTestConfig(123));
    const lowmc::WeightedEventGenerationResult second =
        lowmc::generateWeightedPhotons(
            eventGeneratorToyParticles(),
            eventGeneratorTestConfig(123));

    require(first.photons.size() == 6, "event generator should return requested samples");
    require(first.samples.size() == 6, "event generator should return requested sample records");
    require(
        first.photons.size() == second.photons.size(),
        "same-seed event generator runs should have same photon count");
    for (std::size_t index = 0; index < first.photons.size(); ++index) {
        requireNear(
            first.photons[index].momentum.E,
            second.photons[index].momentum.E,
            1.0e-15,
            "same-seed photon energy");
        requireNear(
            first.photons[index].momentum.px,
            second.photons[index].momentum.px,
            1.0e-15,
            "same-seed photon px");
        requireNear(
            first.photons[index].weight,
            second.photons[index].weight,
            1.0e-15,
            "same-seed photon weight");
        requireNear(
            first.samples[index].proposalDensity,
            second.samples[index].proposalDensity,
            1.0e-15,
            "same-seed proposal density");
    }
}

void testWeightedEventGeneratorChangesWithSeed() {
    const lowmc::WeightedEventGenerationResult first =
        lowmc::generateWeightedPhotons(
            eventGeneratorToyParticles(),
            eventGeneratorTestConfig(123));
    const lowmc::WeightedEventGenerationResult second =
        lowmc::generateWeightedPhotons(
            eventGeneratorToyParticles(),
            eventGeneratorTestConfig(124));

    require(
        std::abs(first.photons.front().momentum.E - second.photons.front().momentum.E) >
            1.0e-15 ||
        std::abs(first.photons.front().weight - second.photons.front().weight) > 1.0e-15,
        "different seeds should change generated photons");
}

void testWeightedEventGeneratorReportsDiagnostics() {
    const lowmc::WeightedEventGenerationResult result =
        lowmc::generateWeightedPhotons(
            eventGeneratorToyParticles(),
            eventGeneratorTestConfig(123));

    require(result.chargeConservation.conserved, "two-pion event should conserve charge");
    require(
        result.samples.size() == result.photons.size(),
        "event generator sample records and photon convenience vector should match");
    for (std::size_t index = 0; index < result.samples.size(); ++index) {
        requireNear(
            result.samples[index].photon.weight,
            result.photons[index].weight,
            1.0e-15,
            "sample record photon should match photon vector");
        require(
            result.samples[index].proposalDensity > 0.0,
            "event generator sample should retain positive proposal density");
        requireNear(
            result.samples[index].photon.momentum.pt(),
            result.samples[index].kT,
            1.0e-15,
            "sample record kT should match photon momentum");
    }
    requireNear(
        result.gridBuildKT,
        std::sqrt(testPhaseSpace().kTMin * testPhaseSpace().kTMax),
        1.0e-15,
        "event generator should build grid at geometric mean kT");
    requireNear(
        result.weightSummary.meanWeight,
        lowmc::summarizePhotons(result.photons).meanWeight,
        1.0e-12,
        "event generator weight summary should summarize returned photons");
}

void testWeightedEventGeneratorAppliesOversamplingFactor() {
    lowmc::WeightedEventGeneratorConfig baseConfig = eventGeneratorTestConfig(123);
    lowmc::WeightedEventGeneratorConfig oversampledConfig = baseConfig;
    oversampledConfig.oversamplingFactor = 3;

    const lowmc::WeightedEventGenerationResult base =
        lowmc::generateWeightedPhotons(eventGeneratorToyParticles(), baseConfig);
    const lowmc::WeightedEventGenerationResult oversampled =
        lowmc::generateWeightedPhotons(eventGeneratorToyParticles(), oversampledConfig);

    require(
        oversampled.samples.size() == base.samples.size() * 3,
        "oversampling factor should multiply returned sample records");
    require(
        oversampled.photons.size() == base.photons.size() * 3,
        "oversampling factor should multiply returned photons");

    for (std::size_t index = 0; index < base.samples.size(); ++index) {
        requireNear(
            oversampled.samples[index].photon.momentum.E,
            base.samples[index].photon.momentum.E,
            1.0e-15,
            "oversampling should preserve sampled photon sequence");
        requireNear(
            oversampled.samples[index].proposalDensity,
            base.samples[index].proposalDensity,
            1.0e-15,
            "oversampling should preserve proposal density sequence");
        requireNear(
            oversampled.samples[index].photon.weight * 3.0,
            base.samples[index].photon.weight,
            1.0e-14,
            "oversampling should scale individual weights down");
    }

    require(
        oversampled.weightSummary.count == oversampled.photons.size(),
        "oversampled weight summary should summarize scaled returned photons");
}

void generateWeightedEventWithInvalidPhaseSpace() {
    lowmc::WeightedEventGeneratorConfig config = eventGeneratorTestConfig(123);
    config.phaseSpace.kTMax = config.phaseSpace.kTMin;
    static_cast<void>(lowmc::generateWeightedPhotons(eventGeneratorToyParticles(), config));
}

void generateWeightedEventWithInvalidOversamplingFactor() {
    lowmc::WeightedEventGeneratorConfig config = eventGeneratorTestConfig(123);
    config.oversamplingFactor = 0;
    static_cast<void>(lowmc::generateWeightedPhotons(eventGeneratorToyParticles(), config));
}

void generateWeightedEventWithTooManyOversampledPhotons() {
    lowmc::WeightedEventGeneratorConfig config = eventGeneratorTestConfig(123);
    config.samples = std::numeric_limits<int>::max();
    config.oversamplingFactor = 2;
    static_cast<void>(lowmc::generateWeightedPhotons(eventGeneratorToyParticles(), config));
}

void testWeightedEventGeneratorRejectsInvalidInputs() {
    requireThrowsInvalidArgument(
        generateWeightedEventWithInvalidPhaseSpace,
        "event generator should reject invalid phase space");
    requireThrowsInvalidArgument(
        generateWeightedEventWithInvalidOversamplingFactor,
        "event generator should reject invalid oversampling factor");
    requireThrowsInvalidArgument(
        generateWeightedEventWithTooManyOversampledPhotons,
        "event generator should reject oversampled sample-count overflow");
}

void testShapeUnweightingReturnsDiagnosticsAndUnitWeightCandidates() {
    const lowmc::ShapeUnweightingConfig config = shapeUnweightingTestConfig(123);
    const lowmc::ShapeUnweightingResult result =
        lowmc::generateShapeUnweightedPhotonCandidates(eventGeneratorToyParticles(), config);

    require(
        result.status == lowmc::ShapeUnweightingStatus::Ok,
        "shape unweighting should succeed with conservative safety factor");
    require(result.chargeConservation.conserved, "shape unweighting should report charge conservation");
    require(
        result.envelope.trialCount == config.trials,
        "shape unweighting diagnostics should report trial count");
    require(
        result.envelope.calibrationCount == config.calibrationTrials,
        "shape unweighting diagnostics should report calibration count");
    require(
        result.envelope.gridCellCount > 0,
        "shape unweighting diagnostics should report grid cell count");
    require(
        result.envelope.probeEvaluations ==
            result.envelope.gridCellCount *
                static_cast<std::size_t>(config.probeBins) *
                static_cast<std::size_t>(config.probeBins),
        "shape unweighting diagnostics should report probe evaluation count");
    require(
        result.envelope.promotedCellCount == 0,
        "fixed shape unweighting should not promote cells");
    require(
        result.envelope.acceptedCount == static_cast<int>(result.candidates.size()),
        "shape unweighting accepted count should match candidates");
    require(
        result.envelope.failedEnvelopeCount == 0,
        "conservative shape unweighting config should have no envelope failures");
    require(
        result.envelope.observedAcceptanceEfficiency >= 0.0 &&
            result.envelope.observedAcceptanceEfficiency <= 1.0,
        "shape unweighting observed acceptance should be a probability");
    require(
        result.envelope.envelopeIntegral > 0.0,
        "shape unweighting should report positive envelope integral");

    for (const lowmc::ShapeUnweightedPhotonCandidate& candidate : result.candidates) {
        requireNear(
            candidate.photon.weight,
            1.0,
            1.0e-15,
            "shape-unweighted candidates should carry unit photon weight");
        requireNear(
            candidate.photon.momentum.pt(),
            candidate.kT,
            1.0e-15,
            "shape-unweighted candidate kT should match photon momentum");
    }
}

void testShapeUnweightingIndexMatchesLinearReference() {
    lowmc::ShapeUnweightingConfig config = shapeUnweightingTestConfig(123);
    config.gridConfig.etaBins = 8;
    config.gridConfig.phiBins = 8;
    config.gridConfig.refinementLevels = 2;
    config.gridConfig.refineFraction = 0.25;
    config.targetedRefinementLevels = 2;
    config.targetedRefinementNeighborScale = 0.5;

    const lowmc::ShapeUnweightingIndexConsistencyDiagnostic diagnostic =
        lowmc::validateShapeUnweightingIndexConsistency(eventGeneratorToyParticles(), config);

    require(diagnostic.consistent, "shape-unweighting index should match linear reference");
    require(
        diagnostic.gridCellCount > 0,
        "shape-unweighting index consistency should build cells");
    require(
        diagnostic.pointLookupChecks >= 5 * diagnostic.gridCellCount,
        "shape-unweighting index consistency should check representative cell points");
    require(
        diagnostic.envelopeSamplingChecks > 0,
        "shape-unweighting index consistency should check envelope sampling");
}

void testShapeUnweightingBoostedCellPolicyPromotesCells() {
    lowmc::ShapeUnweightingConfig config = shapeUnweightingTestConfig(123);
    config.probePolicy = lowmc::ShapeUnweightingProbePolicy::BoostedCells;
    config.boostedProbeBins = 4;
    config.boostedCellSafetyFactor = config.safetyFactor;

    const lowmc::ShapeUnweightingResult result =
        lowmc::generateShapeUnweightedPhotonCandidates(eventGeneratorToyParticles(), config);

    require(
        result.envelope.promotedCellCount > 0,
        "boosted-cell probe policy should promote cells near charged particles");
    require(
        result.envelope.probeEvaluations >
            result.envelope.gridCellCount *
                static_cast<std::size_t>(config.probeBins) *
                static_cast<std::size_t>(config.probeBins),
        "boosted-cell probe policy should report extra probe evaluations");
}

void testShapeUnweightingStrictFailureClearsCandidates() {
    lowmc::ShapeUnweightingConfig config = shapeUnweightingTestConfig(123);
    config.safetyFactor = 1.0e-12;
    config.failurePolicy = lowmc::EnvelopeFailurePolicy::Strict;

    const lowmc::ShapeUnweightingResult result =
        lowmc::generateShapeUnweightedPhotonCandidates(eventGeneratorToyParticles(), config);

    require(
        result.status == lowmc::ShapeUnweightingStatus::EnvelopeFailed,
        "strict shape unweighting should report envelope failure");
    require(
        result.envelope.failedEnvelopeCount > 0,
        "strict shape unweighting should count envelope failures");
    require(
        result.candidates.empty(),
        "strict shape unweighting should clear candidates after envelope failure");
    require(
        result.envelope.acceptedCount == 0,
        "strict shape unweighting should clear accepted count after envelope failure");
}

void testShapeUnweightingDiagnosticFailureKeepsCandidates() {
    lowmc::ShapeUnweightingConfig config = shapeUnweightingTestConfig(123);
    config.safetyFactor = 1.0e-12;
    config.failurePolicy = lowmc::EnvelopeFailurePolicy::Diagnostic;

    const lowmc::ShapeUnweightingResult result =
        lowmc::generateShapeUnweightedPhotonCandidates(eventGeneratorToyParticles(), config);

    require(
        result.status == lowmc::ShapeUnweightingStatus::EnvelopeFailed,
        "diagnostic shape unweighting should report envelope failure");
    require(
        result.envelope.failedEnvelopeCount > 0,
        "diagnostic shape unweighting should count envelope failures");
    require(
        !result.candidates.empty(),
        "diagnostic shape unweighting should keep accepted candidates with failure evidence");
}

void generateShapeUnweightingWithInvalidProbeBins() {
    lowmc::ShapeUnweightingConfig config = shapeUnweightingTestConfig(123);
    config.probeBins = 0;
    static_cast<void>(
        lowmc::generateShapeUnweightedPhotonCandidates(eventGeneratorToyParticles(), config));
}

void generateShapeUnweightingWithInvalidSafetyFactor() {
    lowmc::ShapeUnweightingConfig config = shapeUnweightingTestConfig(123);
    config.safetyFactor = 0.0;
    static_cast<void>(
        lowmc::generateShapeUnweightedPhotonCandidates(eventGeneratorToyParticles(), config));
}

void testShapeUnweightingRejectsInvalidInputs() {
    requireThrowsInvalidArgument(
        generateShapeUnweightingWithInvalidProbeBins,
        "shape unweighting should reject invalid probe bins");
    requireThrowsInvalidArgument(
        generateShapeUnweightingWithInvalidSafetyFactor,
        "shape unweighting should reject invalid safety factor");
}

void testWeightedPhotonSampleRejectsInvalidInputs() {
    requireThrowsInvalidArgument(
        makeSampleWithZeroProposalDensity,
        "weighted photon sample should reject zero proposal density");
    requireThrowsInvalidArgument(
        makeSampleWithInvalidPhaseSpace,
        "weighted photon sample should reject invalid phase space");
}

void testNonConservingChargeSelectionIsReported() {
    const lowmc::LowCurrent current({
        {lowmc::FourVector{2.0, 0.0, 0.0, 1.0}, +1.0, lowmc::InOutSign::Outgoing},
    });

    const lowmc::ChargeConservationDiagnostic diagnostic = current.chargeConservation();
    require(!diagnostic.conserved, "single outgoing charged particle should not conserve charge");
    requireNear(diagnostic.netCharge, 1.0, 1.0e-12, "reported non-conserved charge");
}

void requirePhaseSpaceError(
    const lowmc::PhotonPhaseSpace& phaseSpace,
    lowmc::PhotonPhaseSpaceValidationError expectedError,
    const std::string& message) {
    const lowmc::PhotonPhaseSpaceValidation validation =
        lowmc::validatePhotonPhaseSpace(phaseSpace);
    require(!validation.valid(), message + " should be invalid");
    require(validation.error == expectedError, message + " should report the expected error");
}

void testValidPhotonPhaseSpaceIsAccepted() {
    const lowmc::PhotonPhaseSpace phaseSpace{
        0.02,
        0.2,
        -1.5,
        1.5,
        0.0,
        2.0 * std::acos(-1.0),
    };

    const lowmc::PhotonPhaseSpaceValidation validation =
        lowmc::validatePhotonPhaseSpace(phaseSpace);
    require(validation.valid(), "positive non-empty photon phase space should be valid");
}

void testInvalidPhotonPhaseSpacesAreRejected() {
    const lowmc::PhotonPhaseSpace valid{
        0.02,
        0.2,
        -1.5,
        1.5,
        0.0,
        2.0 * std::acos(-1.0),
    };

    lowmc::PhotonPhaseSpace phaseSpace = valid;
    phaseSpace.kTMin = 0.0;
    requirePhaseSpaceError(
        phaseSpace,
        lowmc::PhotonPhaseSpaceValidationError::NonPositiveKTMin,
        "zero kTMin");

    phaseSpace = valid;
    phaseSpace.kTMax = phaseSpace.kTMin;
    requirePhaseSpaceError(
        phaseSpace,
        lowmc::PhotonPhaseSpaceValidationError::EmptyKTRange,
        "empty kT range");

    phaseSpace = valid;
    phaseSpace.etaMax = phaseSpace.etaMin;
    requirePhaseSpaceError(
        phaseSpace,
        lowmc::PhotonPhaseSpaceValidationError::EmptyEtaRange,
        "empty eta range");

    phaseSpace = valid;
    phaseSpace.phiMax = phaseSpace.phiMin;
    requirePhaseSpaceError(
        phaseSpace,
        lowmc::PhotonPhaseSpaceValidationError::EmptyPhiRange,
        "empty phi range");

    phaseSpace = valid;
    phaseSpace.etaMax = std::numeric_limits<double>::infinity();
    requirePhaseSpaceError(
        phaseSpace,
        lowmc::PhotonPhaseSpaceValidationError::NonFiniteBound,
        "non-finite bound");
}

}  // namespace

int main() {
    testChargeNeutralTwoParticleCurrentIsTransverse();
    testFiniteMassCollinearRegularization();
    testCurrentSquaredMatchesPairSum();
    testLowPhotonDensityUsesLowNormalization();
    testLowPhotonDensityDefaultAlphaIsPositive();
    testWeightedPhotonSampleUsesProposalDensityConvention();
    testLogUniformKTSamplingEndpoints();
    testLogUniformKTDensityFormulas();
    testLogUniformKTRejectsInvalidInputs();
    testUniformEtaPhiSamplingEndpoints();
    testUniformEtaPhiDensityFormula();
    testUniformEtaPhiRejectsInvalidInputs();
    testPhotonMomentumFromKTEtaPhiIsLightlikeAndHasRequestedKT();
    testFixedKTUniformEtaPhiWeightedSampleUsesSliceProposalDensity();
    testFixedKTUniformEtaPhiWeightedSampleRejectsInvalidInputs();
    testAdaptiveEtaPhiProposalSamplesLocalParticlePatch();
    testFixedKTAdaptiveEtaPhiWeightedSampleUsesAdaptiveDensity();
    testSmoothEtaPhiProposalSamplesLocalKernel();
    testFixedKTSmoothEtaPhiWeightedSampleUsesSmoothDensity();
    testGridAdaptiveEtaPhiProposalBuildsRefinedCells();
    testGridAdaptiveEtaPhiSampleUsesPositiveDensity();
    testGridAdaptiveVarianceRefinementScoreBuildsProposal();
    testFixedKTGridAdaptiveEtaPhiWeightedSampleUsesGridDensity();
    testMinimalWeightedPhotonSampleUsesFullProposalDensity();
    testGridAdaptiveWeightedPhotonSampleUsesFullProposalDensity();
    testMinimalWeightedPhotonSampleRejectsInvalidInputs();
    testWeightSummaryComputesBasicDiagnostics();
    testWeightSummaryHandlesEmptyAndZeroWeights();
    testWeightSummaryRejectsInvalidWeights();
    testWeightSummaryCanSummarizePhotonsAndWeightedSamples();
    testEventYieldEstimateIsPositiveAndDeterministic();
    testEventYieldEstimateUsesLogKTMeasure();
    testEventYieldEstimateRejectsInvalidInputs();
    testAdaptiveEventYieldEstimateIsPositiveAndReportsCells();
    testAdaptiveEventYieldRejectsInvalidInputs();
    testPhotonMultiplicityUsesPoissonMeanAndOversamplingWeight();
    testPhotonMultiplicityRejectsInvalidInputs();
    testParticleSelectionConfigDefaultsMatchAdapterContract();
    testParticleSelectionConfigCanToggleOptions();
    testAdapterSelectionDiagnosticSummarizesSelectedParticles();
    testAdapterParticleSelectionUsesPDGChargeAndDefaultPolicy();
    testAdapterParticleSelectionCanDisableIncomingBeams();
    testAdapterParticleSelectionCanIncludeUnstableIntermediates();
    testAdapterParticleSelectionRejectsAmbiguousRecords();
    testPDGChargeLookupHandlesElementaryParticles();
    testPDGChargeLookupHandlesOrdinaryMesons();
    testPDGChargeLookupHandlesOrdinaryBaryons();
    testPDGChargeLookupHandlesIonsAndUnknownIds();
    testWeightedEventGeneratorIsReproducible();
    testWeightedEventGeneratorChangesWithSeed();
    testWeightedEventGeneratorReportsDiagnostics();
    testWeightedEventGeneratorAppliesOversamplingFactor();
    testWeightedEventGeneratorRejectsInvalidInputs();
    testShapeUnweightingReturnsDiagnosticsAndUnitWeightCandidates();
    testShapeUnweightingIndexMatchesLinearReference();
    testShapeUnweightingBoostedCellPolicyPromotesCells();
    testShapeUnweightingStrictFailureClearsCandidates();
    testShapeUnweightingDiagnosticFailureKeepsCandidates();
    testShapeUnweightingRejectsInvalidInputs();
    testWeightedPhotonSampleRejectsInvalidInputs();
    testNonConservingChargeSelectionIsReported();
    testValidPhotonPhaseSpaceIsAccepted();
    testInvalidPhotonPhaseSpacesAreRejected();

    std::cout << "All LowMC core tests passed\n";
    return EXIT_SUCCESS;
}
