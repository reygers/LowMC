# HepMC3 Afterburner Example

The installed copy of this example is self-contained and includes a small
charge-conserving HepMC3 fixture.

Run it from an installed LowMC prefix with:

```sh
python3 hepmc3_afterburner.py --output-dir example-output
```

The script uses the installed LowMC HepMC3 tools from the same prefix by
default, reads the local `charge_conserving_event.hepmc3` fixture, writes a new
HepMC3 file with Low photons, and summarizes the result.

The example uses atomic output by default. Pass `--no-atomic-output` only
when explicitly testing the lower-level diagnostic streaming behavior.
