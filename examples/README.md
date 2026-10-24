# LowMC Examples

The examples form a short learning path. Filenames describe the task; the
`basic` and `advanced` directories indicate the level of complexity.

## Basic

These examples use only the generator-independent core and are built by
default.

| Example | Main idea | Executable |
|---|---|---|
| `basic/evaluate_low_current.cpp` | Evaluate the coherent current, intensity, and density at one photon momentum | `lowmc_example_low_current` |
| `basic/generate_weighted_photons.cpp` | Generate weighted photons from a hand-constructed charged event | `lowmc_example_weighted_photons` |
| `basic/generate_weighted_dielectrons.cpp` | Generate weighted virtual photons and polarization-aware electron-positron decays | `lowmc_example_weighted_dielectrons` |
| `basic/sample_event_multiplicity.cpp` | Integrate `lambda_event` and make one physical Poisson multiplicity draw | `lowmc_example_event_multiplicity` |

## Advanced

These examples introduce shape unweighting or optional event interfaces.

| Example | Main idea | Requirement |
|---|---|---|
| `advanced/shape_unweighting.cpp` | Generate unweighted Low-shape candidates and inspect the envelope diagnostics | Core only |
| `advanced/hepmc3_afterburner.py` | Run the provisional HepMC3 afterburner and summarize its output | `LOWMC_ENABLE_HEPMC3=ON` |
| `advanced/pythia8_weighted_photons.cpp` | Select a generated Pythia event in memory and generate weighted Low photons | `LOWMC_ENABLE_PYTHIA8=ON` |

## Suggested Order

1. `evaluate_low_current.cpp`
2. `generate_weighted_photons.cpp`
3. `generate_weighted_dielectrons.cpp`
4. `sample_event_multiplicity.cpp`
5. `shape_unweighting.cpp`
6. `hepmc3_afterburner.py`
7. `pythia8_weighted_photons.cpp`

The examples use fixed, visible configurations and intentionally avoid general
command-line parsing, CSV scans, plotting, and tuning. Programs for those tasks
belong under `diagnostics/` or `scripts/`.

The shape-unweighting example returns candidates distributed according to the
Low-photon shape. It does not sample the physical event multiplicity. The
HepMC3 example demonstrates the composed yield, multiplicity, shape, and
event-record workflow, which remains provisional in LowMC 0.1.0.
