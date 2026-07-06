#include <HepMC3/Attribute.h>
#include <HepMC3/GenEvent.h>
#include <HepMC3/GenParticle.h>
#include <HepMC3/ReaderAscii.h>

#include <cstdlib>
#include <exception>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct Options {
    std::string inputPath;
    int maxEvents = std::numeric_limits<int>::max();
};

void printUsage(std::ostream& out) {
    out << "Usage: lowmc_hepmc3_output_summary_diagnostic [--max-events N] <input.hepmc3>\n"
        << "\n"
        << "Reads HepMC3 ASCII events and prints a compact read-only summary of\n"
        << "LowMC afterburner photons and LowMC event-level metadata.\n"
        << "No photons are generated and the input file is not modified.\n";
}

std::string readStringOption(int& index, int argc, char** argv, const std::string& optionName) {
    if (index + 1 >= argc) {
        throw std::invalid_argument(optionName + " requires a value");
    }
    return argv[++index];
}

int parsePositiveInt(const std::string& value, const std::string& optionName) {
    std::size_t parsed = 0;
    const long parsedValue = std::stol(value, &parsed);
    if (parsed != value.size() ||
        parsedValue <= 0 ||
        parsedValue > std::numeric_limits<int>::max()) {
        throw std::invalid_argument(optionName + " expects a positive integer");
    }
    return static_cast<int>(parsedValue);
}

Options parseOptions(int argc, char** argv) {
    Options options;

    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument == "--help" || argument == "-h") {
            printUsage(std::cout);
            std::exit(EXIT_SUCCESS);
        }
        if (argument == "--max-events") {
            options.maxEvents =
                parsePositiveInt(readStringOption(index, argc, argv, argument), argument);
            continue;
        }
        if (!argument.empty() && argument.front() == '-') {
            throw std::invalid_argument("unknown option: " + argument);
        }
        if (!options.inputPath.empty()) {
            throw std::invalid_argument("expected exactly one input file path");
        }
        options.inputPath = argument;
    }

    if (options.inputPath.empty()) {
        throw std::invalid_argument("missing input file path");
    }
    return options;
}

[[nodiscard]] bool isEmptyTerminalEvent(const HepMC3::GenEvent& event) {
    return event.particles().empty() && event.vertices().empty() && event.beams().empty();
}

[[nodiscard]] bool isLowMCPhoton(const HepMC3::ConstGenParticlePtr& particle) {
    return particle &&
           particle->pid() == 22 &&
           particle->status() == 1 &&
           particle->attribute_as_string("LowMC:origin") == "afterburner";
}

[[nodiscard]] double doubleAttributeValue(
    const std::shared_ptr<HepMC3::DoubleAttribute>& attribute,
    double defaultValue) {
    if (!attribute) {
        return defaultValue;
    }
    return attribute->value();
}

void printSummary(int eventIndex, const HepMC3::GenEvent& event) {
    int lowmcPhotonCount = 0;
    int photonWeightAttributeCount = 0;
    double photonWeightSum = 0.0;
    std::string photonMode = event.attribute_as_string("LowMC:photon_mode");
    double enhancementFactor = doubleAttributeValue(
        event.attribute<HepMC3::DoubleAttribute>("LowMC:enhancement_factor"),
        0.0);

    for (const HepMC3::ConstGenParticlePtr& particle : event.particles()) {
        if (!isLowMCPhoton(particle)) {
            continue;
        }
        ++lowmcPhotonCount;
        if (photonMode.empty()) {
            photonMode = particle->attribute_as_string("LowMC:photon_mode");
        }
        if (enhancementFactor == 0.0) {
            enhancementFactor = doubleAttributeValue(
                particle->attribute<HepMC3::DoubleAttribute>("LowMC:enhancement_factor"),
                0.0);
        }
        const std::shared_ptr<HepMC3::DoubleAttribute> photonWeight =
            particle->attribute<HepMC3::DoubleAttribute>("LowMC:weight");
        if (photonWeight) {
            ++photonWeightAttributeCount;
        }
        photonWeightSum += doubleAttributeValue(photonWeight, 1.0);
    }

    const bool hasEventMetadata =
        !event.attribute_as_string("LowMC:photon_mode").empty() ||
        static_cast<bool>(event.attribute<HepMC3::DoubleAttribute>("LowMC:enhancement_factor")) ||
        static_cast<bool>(event.attribute<HepMC3::DoubleAttribute>("LowMC:weight"));

    std::cout << "event_index=" << eventIndex
              << " lowmc_photons=" << lowmcPhotonCount
              << " photon_mode=" << (photonMode.empty() ? "none" : photonMode)
              << " enhancement_factor=" << enhancementFactor
              << " photon_weight_sum=" << photonWeightSum
              << " photon_weight_attributes=" << photonWeightAttributeCount
              << " event_metadata=" << (hasEventMetadata ? 1 : 0)
              << '\n';
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Options options = parseOptions(argc, argv);
        HepMC3::ReaderAscii reader(options.inputPath);
        if (reader.failed()) {
            std::cerr << "FAILED: could not open HepMC3 input file: "
                      << options.inputPath << '\n';
            return EXIT_FAILURE;
        }

        int eventsRead = 0;
        while (eventsRead < options.maxEvents && !reader.failed()) {
            HepMC3::GenEvent event;
            if (!reader.read_event(event)) {
                break;
            }
            if (isEmptyTerminalEvent(event)) {
                break;
            }
            printSummary(eventsRead, event);
            ++eventsRead;
        }

        if (eventsRead == 0) {
            std::cerr << "FAILED: no events read from HepMC3 input file: "
                      << options.inputPath << '\n';
            return EXIT_FAILURE;
        }
    } catch (const std::exception& error) {
        std::cerr << "FAILED: " << error.what() << '\n';
        printUsage(std::cerr);
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
