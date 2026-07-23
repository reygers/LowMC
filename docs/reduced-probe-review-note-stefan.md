# Reduced-Probe Review Note For Stefan

## Context

LowMC's current HepMC3 afterburner diagnostic uses a full `11 x 11` base-probe
lattice as the boosted-cell reference. Reduced base probing is being studied
only as a diagnostic optimization. It is not a default, and it does not change
the Low density. The compensation factor multiplies reduced-cell envelope
estimates to recover safety margin when the reduced lattice misses narrow
structure.

The current question is whether `5 x 5` reduced base probing with compensation
factor `1.50` is a reasonable conservative candidate for independent review.
Factor `1.30` is kept as the speed/efficiency comparison, not as the preferred
safety candidate.

## Recommendation

Present `5 x 5`, factor `1.50` as the only current conservative reduced-probe
candidate.

Keep `5 x 5`, factor `1.30` in the comparison table because it is useful for
understanding the speed/efficiency tradeoff, but do not present it as having a
clean safety record.

Do not recommend factor `1.35`: it failed the per-seed full-reference criterion
in the first thirty-seed block.

Reduced probing should remain diagnostic-only until the safety/efficiency
criterion is reviewed and at least one further validation scope is agreed.

## Evidence

All runs below used one Pythia-produced event per seed, enhancement factor `5`,
multiplicity seed `13583`, shape seed `97535`, the default boosted-cell path,
and the same photon accounting checks.

| Seed block | Configuration | Time | Max trial/bound | Photons | Envelope failures | Per-seed overshoots |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| `24681..24710` | full `11 x 11` | `79.87 s` | `0.384364` | `174/174` | `0` | reference |
| `24681..24710` | `5 x 5`, factor `1.30` | `36.93 s` | `0.361652` | `174/174` | `0` | `0` |
| `24681..24710` | `5 x 5`, factor `1.35` | `37.12 s` | `0.373041` | `174/174` | `0` | `1` |
| `24681..24710` | `5 x 5`, factor `1.50` | `36.76 s` | `0.321687` | `174/174` | `0` | `0` |
| `24711..24740` | full `11 x 11` | `151.36 s` | `0.395003` | `313/313` | `0` | reference |
| `24711..24740` | `5 x 5`, factor `1.30` | `67.06 s` | `0.405221` | `313/313` | `0` | `4` |
| `24711..24740` | `5 x 5`, factor `1.50` | `67.38 s` | `0.351192` | `313/313` | `0` | `0` |

The fresh held-out block `24711..24740` was harder than the first block:
selected outgoing charged multiplicity reached `274`, and requested photons
summed to `313`.

In that fresh block, factor `1.30` overshot the matching full-reference ratio
on seeds `24711`, `24713`, `24725`, and `24729`. The largest overshoot was
seed `24713`, where factor `1.30` gave `0.405221` versus full-reference
`0.395003`.

Factor `1.50` had no per-seed overshoots in either block. In the fresh block,
its narrowest margin relative to the matching full reference was `0.040705` at
seed `24711`.

## Efficiency Tradeoff

The main cost of factor `1.50` is not wall time in these runs. It is lower
accepted-trial efficiency due to a more inflated envelope.

For `24711..24740`, all three configurations used `662000` trial points:

| Configuration | Accepted trials |
| --- | ---: |
| full `11 x 11` | `166475` |
| `5 x 5`, factor `1.30` | `138274` |
| `5 x 5`, factor `1.50` | `120172` |

Thus factor `1.50` keeps the reduced-probe wall time but gives fewer accepted
candidates. This is the expected price of a more conservative envelope.

## Interpretation

The compensation factor is an envelope-safety device. It inflates reduced-cell
bound estimates; it does not alter the Low-current density, photon
normalization, multiplicity convention, or output metadata.

The relevant validation criterion is therefore not just "no envelope failures"
in the reduced run. The stronger criterion used here is:

- same requested/retained photon counts as the full reference
- zero envelope failures and zero summary mismatches
- no per-seed overshoot relative to the matching full `11 x 11` reference
- acceptable accepted-trial efficiency
- explicit diagnostics for any event close to the bound

On this criterion, factor `1.50` passes the current two-block evidence; factor
`1.30` does not.

## Questions For Review

1. Is the "no per-seed overshoot relative to full `11 x 11`" criterion the
   right current safety criterion, or should the review use a different margin?
2. Is the accepted-trial penalty of factor `1.50` acceptable for a non-default
   diagnostic option?
3. Before exposing reduced probing as a documented non-default diagnostic
   option, should the next validation use more fresh Pythia seeds, different
   Pythia processes, or targeted topology classes?
4. Are there simpler diagnostics we should add so future failures identify
   whether the issue is calibration, boosted-cell promotion, or residual
   reduced-probe under-bounding?

## Current Proposed Decision

Keep the full `11 x 11` boosted-cell reference as the default.

Keep reduced probing diagnostic-only.

Use `5 x 5`, factor `1.50` as the conservative reduced-probe candidate for
independent review.

Keep `5 x 5`, factor `1.30` only as the performance comparison.
