# Reference data

`alice_fig4_hepdata.csv` contains the ordinary data points from Tables 4 and
5 of HEPData record `ins1798514` (version 1):

- <https://doi.org/10.17182/hepdata.111331.v1/t4>
- <https://doi.org/10.17182/hepdata.111331.v1/t5>

The table preserves the separate statistical and systematic uncertainties.
The bins published only as upper limits are not present in the HEPData YAML
tables and are therefore not invented here. The additional 5% luminosity
uncertainty quoted by ALICE is global and is handled as a plot annotation.

`alice_fig4_published_bremsstrahlung.csv` contains the green hadronic-
bremsstrahlung curves exactly as drawn in Figure 4 of the ALICE publication.
The vector paths were extracted from the publication PDF with `pdftocairo`
and transformed from page coordinates using the two plotted axes. These rows
are a trace of the published line, not a reconstruction of the underlying
calculation and not an independent set of theory values supplied by its
authors. This distinction matters because the published transverse-momentum
curve is not reproduced exactly by the BPP Eq. (9a) dipole term alone.
