# LowMC

*A Monte Carlo afterburner for coherent soft-photon emission based on Low's theorem.*

## Motivation

Soft-photon production has been measured in hadronic collisions for decades.
Experimental analyses usually compare the observed photon yield to the
event-by-event prediction obtained from Low's soft-photon theorem.

LowMC turns this prediction into a generator-independent library plus optional
adapter and diagnostic tools.

Given an arbitrary event record, LowMC computes the coherent electromagnetic
current

Jμ(k)

of all charged particles and generates soft photons according to the leading
Low formula.

## Features

- event-wide coherent radiation
- generator-independent core
- HepMC3 support
- Pythia8 interface
- weighted generation and provisional unweighted HepMC3 output
- future support for virtual photons and dielectrons

## Status

Early development.

The project is currently implementing the leading-order soft-photon generator.

## Example

The first user-facing workflow is analysis-oriented:

```text
input event or selected charged particles
-> LowMC core
-> weighted Low-photon samples
-> user histograms and diagnostics
```

The provisional HepMC3 diagnostic also provides an event-record afterburner
workflow while production validation continues:

```text
input.hepmc3 -> LowMC afterburner -> output.hepmc3
```

## Documentation

- docs/DESIGN.md
- docs/USAGE.md
- docs/CODEX.md
- docs/STATUS.md

## References

Low, F.E.
Bremsstrahlung of very low-energy quanta in elementary particle collisions.

...

## License

MIT
