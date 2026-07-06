# REFERENCES

## Core Physics

- F. E. Low, "Bremsstrahlung of very low-energy quanta in elementary particle
  collisions."
- T. H. Burnett and N. M. Kroll, "Extension of the Low Soft-Photon Theorem."
- V. Balek, N. Pisutova, and J. Pisut, "A Search of a Mechanism Responsible
  for Bremsstrahlung Enhancement in Hadronic Reactions. III. Low Mass Dilepton
  Production", Acta Physica Slovaca 41, no. 4 (1991).
  Local working note:
  `literature/Balek-Pisutova-Pisut-1991-soft-dileptons/balek-pisutova-pisut-soft-dileptons-summary.md`

## Related Concepts

- Soft-photon theorem
- Coherent electromagnetic radiation from external charged particles
- Infrared factorization
- Finite-mass collinear regularization
- Soft-photon exponentiation
- Kroll-Wada conversion for virtual photons and dielectrons
- Adaptive Monte Carlo integration
- Importance sampling
- Cellular adaptive event generation

## Tools And Interfaces To Consider Later

- HepMC3 event records
- Pythia8 event generation
- PHOTOS and other QED-radiation tools

## Related Software References

- VEGAS adaptive Monte Carlo integration
- FOAM cellular Monte Carlo event generator
- HepMC3 manual, Doxygen, and source headers for event-record concepts,
  particle status, momentum units, and beam-particle access.
  - Upstream mirror headers checked: `GenEvent.h`, `GenParticle.h`,
    `FourVector.h`, `Units.h`, and `Data/GenParticleData.h` from
    `https://github.com/hep-mirrors/hepmc3`.
- PDG Monte Carlo particle-numbering scheme for PDG id structure, quark-content
  digits, special codes, private-code ranges, and ion code `10LZZZAAAI`.
  - `https://pdg.lbl.gov/2023/reviews/rpp2023-rev-monte-carlo-numbering.pdf`
- Pythia8 manual and UserHooks documentation
- Pythia8 event-record manual for particle status helpers, momenta, charge, and
  beam/documentation entries.
- Pythia8 HepMC3 output interface
- PHOTOS Monte Carlo documentation and validation papers

## Open Bibliography Questions

- Full bibliographic data for Low's paper should be added.
- Full bibliographic data for the Burnett-Kroll paper should be added.
- Full bibliographic data for the Kroll-Wada formula should be added.
- Which experimental soft-photon measurements should be used for validation?
- Which references should define the preferred exponentiation formalism?
- Which generator-interface references are needed for HepMC3 and Pythia8
  adapters?
