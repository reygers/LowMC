#include "lowmc/hepmc3/HepMC3Adapter.h"

#include <HepMC3/Attribute.h>
#include <HepMC3/FourVector.h>
#include <HepMC3/GenEvent.h>
#include <HepMC3/GenParticle.h>
#include <HepMC3/GenVertex.h>
#include <HepMC3/Units.h>

#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

[[nodiscard]] std::string sourceLabel(const HepMC3::GenParticle& particle) {
    return "hepmc3:id=" + std::to_string(particle.id()) +
           ",pid=" + std::to_string(particle.pid()) +
           ",status=" + std::to_string(particle.status());
}

[[nodiscard]] bool containsParticlePointer(
    const std::vector<const HepMC3::GenParticle*>& particles,
    const HepMC3::GenParticle* candidate) {
    for (const HepMC3::GenParticle* particle : particles) {
        if (particle == candidate) {
            return true;
        }
    }
    return false;
}

void addStringAttribute(
    HepMC3::GenParticle& particle,
    const std::string& name,
    const std::string& value) {
    particle.add_attribute(name, std::make_shared<HepMC3::StringAttribute>(value));
}

void addIntAttribute(
    HepMC3::GenParticle& particle,
    const std::string& name,
    int value) {
    particle.add_attribute(name, std::make_shared<HepMC3::IntAttribute>(value));
}

void addDoubleAttribute(
    HepMC3::GenParticle& particle,
    const std::string& name,
    double value) {
    particle.add_attribute(name, std::make_shared<HepMC3::DoubleAttribute>(value));
}

void addStringAttribute(
    HepMC3::GenEvent& event,
    const std::string& name,
    const std::string& value) {
    event.add_attribute(name, std::make_shared<HepMC3::StringAttribute>(value));
}

void addDoubleAttribute(
    HepMC3::GenEvent& event,
    const std::string& name,
    double value) {
    event.add_attribute(name, std::make_shared<HepMC3::DoubleAttribute>(value));
}

[[nodiscard]] HepMC3::FourVector toHepMC3Momentum(
    const lowmc::FourVector& momentum,
    HepMC3::Units::MomentumUnit outputUnit) {
    HepMC3::FourVector hepmcMomentum{
        momentum.px,
        momentum.py,
        momentum.pz,
        momentum.E,
    };
    HepMC3::Units::convert(hepmcMomentum, HepMC3::Units::GEV, outputUnit);
    return hepmcMomentum;
}

}  // namespace

namespace lowmc::hepmc3 {

const char* hepmc3AdapterBuildMarker() {
    return "lowmc_hepmc3_adapter";
}

AdapterParticleCandidate makeAdapterParticleCandidate(
    const HepMC3::GenParticle& particle,
    HepMC3::Units::MomentumUnit momentumUnit,
    const HepMC3ParticleClassification& classification) {
    HepMC3::FourVector momentum = particle.momentum();
    HepMC3::Units::convert(momentum, momentumUnit, HepMC3::Units::GEV);

    return AdapterParticleCandidate{
        FourVector{momentum.e(), momentum.px(), momentum.py(), momentum.pz()},
        particle.pid(),
        classification.isIncomingBeam,
        classification.isFinalStable,
        classification.isBeamRemnant,
        classification.sourceLabel,
    };
}

std::vector<AdapterParticleCandidate> collectAdapterParticleCandidates(
    const HepMC3::GenEvent& event) {
    std::vector<AdapterParticleCandidate> candidates;
    std::vector<const HepMC3::GenParticle*> beamParticles;
    const HepMC3::Units::MomentumUnit momentumUnit = event.momentum_unit();

    for (const HepMC3::ConstGenParticlePtr& beam : event.beams()) {
        if (!beam) {
            continue;
        }
        beamParticles.push_back(beam.get());
        candidates.push_back(makeAdapterParticleCandidate(
            *beam,
            momentumUnit,
            HepMC3ParticleClassification{
                true,
                false,
                false,
                sourceLabel(*beam),
            }));
    }

    for (const HepMC3::ConstGenParticlePtr& particle : event.particles()) {
        if (!particle || containsParticlePointer(beamParticles, particle.get())) {
            continue;
        }
        if (particle->status() != 1 || particle->end_vertex()) {
            continue;
        }

        candidates.push_back(makeAdapterParticleCandidate(
            *particle,
            momentumUnit,
            HepMC3ParticleClassification{
                false,
                true,
                false,
                sourceLabel(*particle),
            }));
    }

    return candidates;
}

AdapterParticleSelectionResult selectChargedParticles(
    const HepMC3::GenEvent& event,
    const ParticleSelectionConfig& config) {
    return selectAdapterParticles(collectAdapterParticleCandidates(event), config);
}

HepMC3LowPhotonOutputResult appendUnweightedLowPhotons(
    HepMC3::GenEvent& event,
    const std::vector<ShapeUnweightedPhotonCandidate>& photons,
    const HepMC3LowPhotonOutputOptions& options) {
    if (!std::isfinite(options.enhancementFactor) || options.enhancementFactor <= 0.0) {
        throw std::invalid_argument(
            "HepMC3 Low-photon output requires a positive finite enhancementFactor");
    }

    HepMC3LowPhotonOutputResult result;
    if (photons.empty()) {
        return result;
    }

    const HepMC3::Units::MomentumUnit outputUnit = event.momentum_unit();
    const bool enhanced = options.enhancementFactor != 1.0;
    const std::string photonMode = enhanced ? "enhanced" : "unweighted";
    const double photonWeight = 1.0 / options.enhancementFactor;
    addStringAttribute(event, "LowMC:photon_mode", photonMode);
    addDoubleAttribute(event, "LowMC:enhancement_factor", options.enhancementFactor);
    if (enhanced) {
        addDoubleAttribute(event, "LowMC:weight", photonWeight);
    }

    int sampleIndex = 0;
    for (const ShapeUnweightedPhotonCandidate& photon : photons) {
        auto particle = std::make_shared<HepMC3::GenParticle>(
            toHepMC3Momentum(photon.photon.momentum, outputUnit),
            22,
            1);
        particle->unset_generated_mass();
        event.add_particle(particle);

        addStringAttribute(*particle, "LowMC:origin", "afterburner");
        addStringAttribute(*particle, "LowMC:model", "leading-low");
        addStringAttribute(*particle, "LowMC:recoil", "none");
        addStringAttribute(*particle, "LowMC:photon_mode", photonMode);
        addDoubleAttribute(*particle, "LowMC:enhancement_factor", options.enhancementFactor);
        if (enhanced) {
            addDoubleAttribute(*particle, "LowMC:weight", photonWeight);
        }
        if (!options.version.empty()) {
            addStringAttribute(*particle, "LowMC:version", options.version);
        }
        if (options.hasLambdaEvent) {
            addDoubleAttribute(*particle, "LowMC:lambda_event", options.lambdaEvent);
        }
        if (!options.validationClass.empty()) {
            addStringAttribute(*particle, "LowMC:validation_class", options.validationClass);
        }
        if (!options.shapeStatus.empty()) {
            addStringAttribute(*particle, "LowMC:shape_status", options.shapeStatus);
        }
        if (options.hasMaxTrialBoundRatio) {
            addDoubleAttribute(*particle, "LowMC:max_trial_bound_ratio", options.maxTrialBoundRatio);
        }
        addIntAttribute(*particle, "LowMC:sample_index", sampleIndex);
        addDoubleAttribute(*particle, "LowMC:kT", photon.kT);
        addDoubleAttribute(*particle, "LowMC:eta", photon.eta);
        addDoubleAttribute(*particle, "LowMC:phi", photon.phi);

        ++sampleIndex;
    }

    result.photonsAdded = sampleIndex;
    return result;
}

}  // namespace lowmc::hepmc3
