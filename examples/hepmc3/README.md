# HepMC3 Quickstart Example

This directory is the installed, self-contained HepMC3 example for LowMC.

Run it from an installed LowMC prefix with:

```sh
python3 hepmc3_quickstart.py --output-dir quickstart-output
```

The script uses the installed LowMC HepMC3 tools from the same prefix by
default, reads the local `charge_conserving_event.hepmc3` fixture, writes a new
HepMC3 file with Low photons, and summarizes the result.

The quickstart uses atomic output by default. Pass `--no-atomic-output` only
when explicitly testing the lower-level diagnostic streaming behavior.
