# DECISIONS

## Accepted

### Generator-independent core

The physics core must not depend on Pythia, HepMC, or another event generator.
Adapters can be added around the core later.

HepMC3 is the preferred first event-record interface because it keeps LowMC
portable across generators and stored event samples. A Pythia8 integration
should remain a thin adapter around the same core.

The first HepMC3 implementation should be optional at build time and live
outside `lowmc_core`. The core must keep building without HepMC3. Exact HepMC3
status-code, beam-particle, charge, and unit mappings must be documented before
generator-specific code is added.

HepMC3 particle charge should be inferred from PDG id by default. The lookup
should use integer `3Q` internally, combine a small explicit table for
elementary and special particles with algorithmic decoding for ordinary mesons
and baryons, and reject unsupported ids with diagnostics rather than guessing.
LowMC should not require Pythia8, ROOT, HepPDT, or a generated particle table
for the first HepMC3 adapter. Optional user overrides or generated tables can
be added later if needed.

The first direct Pythia8 adapter should follow the same rule: Pythia-specific
code remains outside `lowmc_core`, translates `Pythia8::Event` records into
generator-independent `AdapterParticleCandidate` objects, and then reuses
`selectAdapterParticles(...)`. Charge lookup should use the LowMC PDG helper by
default so HepMC3, Pythia8, and direct/manual inputs share one charge
convention. Pythia particle-data charge may be useful as an optional diagnostic
cross-check or future override, but it should not be the default source of
truth for the first adapter boundary.

### Leading Low theorem first

The initial model is the leading-order Low theorem for coherent soft-photon
emission. Higher-order corrections are outside the first implementation.

### Event-wide coherent current

LowMC computes one current from all selected charged particles in the event.
Interference terms are therefore included automatically through `S(k) = -J^2`.

### Physical four-momenta

Charged particles are represented by their physical four-momenta. Finite masses
regularize the collinear region, so no additional angular regulator should be
introduced in the physics model.

### Compute `J` before `S`

The implementation should compute `J^mu(k)` first and then form `S(k) = -J^2`.
This preserves the intended `O(N)` evaluation cost per photon direction.

### Default particle selection

The default coherent current should include incoming beam particles and final
stable charged particles.

Incoming particles use `eta_i = -1`; outgoing final-state particles use
`eta_i = +1`.

The symbol `eta_i` is used in formulas. Code and public APIs should call the
same quantity `inOutSign` to avoid confusion with pseudorapidity.

Particle charges `Q_i` are measured in units of the positron charge `e`.
Four-momenta use positive energies for both incoming and outgoing particles; the
incoming sign is represented only by `inOutSign`.

For HepMC-style event records, the adapter should treat stable charged
final-state particles as the default outgoing selection, using event-record
finality such as no end vertex and, where appropriate, conventional stable
status codes. The exact mapping belongs in the adapter, not in the physics core.

The first HepMC3-to-weighted-photon diagnostic should not write photons back
into the HepMC3 event record. It should report adapter-selection diagnostics
and weighted-photon summaries as text. For events with nonzero signed charge
sum, it should print `photon_status=skipped_nonconserved_charge` and continue
with later events, rather than silently repairing the selection or failing the
whole file.

The first production HepMC3 output writer should wait for unweighted Low
photons. Weighted photon samples should not be written as ordinary final-state
particles by default, because downstream analyses could mistake sampling
records for physical photons. If a weighted HepMC3 writer is added for
diagnostics, it must be explicitly marked as diagnostic output and attach the
sample weight and sampling coordinates as attributes.

The first HepMC3 output-writer contract is settled for unweighted output.
Added Low photons should be appended to the input event as stable PDG `22`,
status-`1` particles. HepMC3 should assign normal particle ids/barcodes; LowMC
sample indices are provenance attributes, not HepMC3 ids. The writer must
preserve the input event record, charged-particle momenta, existing vertices,
and event weights exactly. It must not apply recoil.

The durable first HepMC3 topology is attributed stable root photons: each
LowMC-added photon is stored as an ordinary stable photon with `LowMC:`
provenance attributes directly on the photon, while all original event
vertices remain unchanged. This is accepted because it round-trips cleanly
through the local HepMC3 ASCII writer/reader and does not imply a physical
production vertex or recoil prescription. The originally planned dedicated
no-incoming afterburner vertex is deferred because it did not round-trip
cleanly in the local HepMC3 ASCII test.

The default HepMC3 writer mode is physical output with enhancement factor `1`
and ordinary unweighted photons. It should attach `LowMC:`-prefixed provenance
attributes to the appended photons, including LowMC origin, leading-Low model,
no-recoil policy, photon mode, sampled `k_T`, `eta`, and `phi`, and available
validation diagnostics. If adapter selection, charge conservation, yield
validation, shape-envelope validation, or accepted-candidate availability
fails, the writer must write no Low photons for that event or fail by explicit
user policy. It must not write partial, biased, weighted, or ambiguous photons
as ordinary final-state particles.

An explicit enhanced HepMC3 output mode is accepted and implemented in the
afterburner diagnostic. In that mode LowMC draws
`N_gamma ~ Poisson(enhancementFactor * lambda_event)`. The durable metadata
fallback is per added photon: `LowMC:photon_mode = "enhanced"`,
`LowMC:enhancement_factor = F`, and `LowMC:weight = 1 / F`. The first compact
duplicate is event-level metadata with the same mode, enhancement factor, and,
in enhanced mode, weight on events that receive LowMC photons. Skipped events
remain unchanged. Direct photon counting gives the enhanced sample yield, while
filling histograms with `LowMC:weight` gives the physical Low yield. Run-level
metadata is deferred until a stable file-level `GenRunInfo` convention is
needed.

The current `lowmc_hepmc3_afterburner_diagnostic` keeps streaming output as its
default behavior and can leave readable partial diagnostic output if
`--failure-policy fail` aborts after earlier events were written. This is
accepted only for the diagnostic default. The same diagnostic also supports
explicit atomic output through `--atomic-output`, writing to `<output>.tmp` and
publishing the requested output path only after all events succeed. A stable
production afterburner command or alias should use atomic output by default.
The repository quickstart helper follows the production-style default and runs
with atomic output unless `--no-atomic-output` is requested explicitly.

Weighted histogramming outside event-record output may still use an explicit
integer `oversamplingFactor`. Oversampling multiplies the number of returned
weighted photon samples and divides each photon weight by the same factor. It
is a variance-reduction and histogram-smoothing convenience. It should remain
conceptually distinct from enhanced HepMC3 event samples, even though both use
the same weight idea.

The first unweighting work should remain diagnostic-only. It may use empirical
per-event envelopes to test accept/reject efficiency and envelope failures, but
accepted diagnostic photons should not yet be written as physical HepMC3
particles. Production event output requires a settled multiplicity and
normalization convention in addition to a stable shape-unweighting envelope.

The first generator-independent shape-unweighting API should remain separate
from both multiplicity and event-record writing. It should return
shape-unweighted photon candidates plus envelope diagnostics, not physical
event-record photons. The first policy should use the probed-cell envelope with
`11 x 11` interior probes, safety factor `3`, and the suspect-cell hybrid
fallback. Event-record production should eventually use a strict failure
policy: if any trial exceeds its envelope, ordinary unweighted photons from
that event must not be silently written. Diagnostic modes may still return
accepted candidates together with failure diagnostics for algorithm studies.

The API diagnostic `lowmc_shape_unweighting_diagnostic` is the preferred path
for future production-like shape-unweighting stress studies. The older
standalone `lowmc_unweighting_diagnostic` may remain temporarily as an internal
envelope-development tool, but it should not be the user-facing validation
path. The API path uses the corrected cell-local envelope construction: sample
cells proportional to `q_c M_c`, sample uniformly within the selected cell, and
accept with probability `w / M_c`.

The first physical multiplicity convention for unweighted event-record photons
is Poisson sampling with mean equal to the event-level integrated Low yield in
the configured `PhotonPhaseSpace`,
`N_gamma ~ Poisson(lambda_event)`. Conditional on `N_gamma`, photon momenta
should be drawn from the normalized Low shape for that event. This is a first
independent-emission convention built from the leading-soft density; it is not
a recoil model or a full multi-photon exponentiation treatment.

`oversamplingFactor` remains diagnostic-only. For diagnostic histograms,
LowMC may sample `N_sample ~ Poisson(oversamplingFactor * lambda_event)` and
assign each generated photon analysis weight `1 / oversamplingFactor`, so the
expected weighted yield remains `lambda_event`. This must not change
`lambda_event`, the physical Poisson mean, or the number of ordinary unweighted
photons written to an event record.

For varying cell bounds `M_c`, sampling from the original proposal and accepting
with `w / M_c` is not enough by itself; it remains useful only as a bound and
efficiency diagnostic.

For the corrected transverse `p = 1000 GeV` stress case, a large fixed safety
factor can make the tested envelope hold, but the resulting acceptance can be
too low for production use. A fixed `21 x 21` probe lattice with safety factor
`20` is a conservative fallback baseline, not the preferred default. The
preferred next direction is adaptive or targeted probing of boosted/suspect
cells.

The first boosted-cell probe policy is accepted only as a diagnostic prototype.
It may promote cells near charged-particle directions or cells where calibration
exceeds the base probe maximum, then use denser probes and a larger local
safety factor for those cells. This preserves the corrected cell-local
unweighting construction. The first tested setting cuts probe evaluations
relative to full `21 x 21` probing and holds in the small `p = 1000 GeV` scan,
but its acceptance remains far below one percent, so it is not a production
default.

LowMC should describe boosted-particle limitations through validated operating
ranges, stress tests, and event-level diagnostics rather than through a hard
charged-particle momentum cutoff. The physics formula has no artificial
momentum limit, but the numerical sampling problem becomes harder as `m / E`
gets small. For unweighted event-record output, envelope failure or very poor
efficiency must be reported explicitly; weighted generation remains the
fallback for difficult events.

For the boosted-cell prototype, refining the event-specific grid can be a more
useful lever than simply promoting more cells or using denser promoted-cell
probes. The first refined-grid candidate, `48 x 48` with four refinement
levels and promoted `21 x 21` probes, lowers the `p = 1000 GeV` envelope
integral while keeping the five-seed scan failure-free. It is accepted as a
diagnostic candidate only: the probe cost is high and the accepted statistics
remain too low for production unweighted output.

Targeted charged-direction refinement is now preferred over global refined-grid
growth for the next diagnostic prototype. It splits only cells containing
selected charged-particle directions before the envelope is built, preserving
the corrected sampler that draws cells proportional to `q_c M_c` and accepts
with `w / M_c`. The current best tested setting,
`targeted_l4_p11_to_p21_s3_local5_n0`, gives much better `p = 1000 GeV`
acceptance and lower cost than the global `48 x 48`, level-4 candidate. It is
still a diagnostic candidate, not a production default, until higher-statistics
and broader-topology validations are complete.

After the first multi-seed candidate validation, `targeted_l4_p11_to_p21_s3_local5_n0`
is accepted as the default diagnostic shape-unweighting candidate. It remains a
diagnostic default only: production event-record output still requires a
settled yield integration path, strict envelope-failure behavior, broader event
topologies, and documented release thresholds.

Named validation-script presets are preferred over changing the low-level
diagnostic executable defaults. The candidate-validation script may default to
the accepted targeted candidate and expose a comparison preset for older
baselines, while `lowmc_shape_unweighting_diagnostic` should continue to expose
explicit tuning knobs in its command-line interface.

The first less-symmetric validation topology is an explicit two-pion case with
unequal momenta and non-back-to-back directions. It is accepted as a diagnostic
case because it is still charge-neutral and inspectable by eye. It is not a
replacement for later multi-particle or generator-event validation.

After the first topology-dependence scan, the preferred diagnostic
shape-unweighting candidate is
`targeted_l5_halo0.5_p11_to_p21_s3_local5`. The small targeted split halo is
kept because it improves the asymmetric two-pion case, and the fifth targeted
split level is kept because it improves the symmetric `p = 1000 GeV` stress
case. This supersedes `targeted_l4_p11_to_p21_s3_local5_n0` as the diagnostic
default, but it is still not a production majorant guarantee.

For the current toy validation suite, provisional envelope gates may be used:
zero envelope failures, zero fallback uses, and maximum trial/bound ratio below
one with visible margin. Shape-comparison L1 and projection RMS values should
continue to be reported, but they are not release-blocking yet because the
fixed-grid references can be resolution limited for narrow boosted structures.
These gates are diagnostic gates for the current toy suite, not production
event-record guarantees.

The executable gate interface belongs in the candidate-validation script rather
than in the lower-level C++ diagnostic. The script can aggregate over cases and
seeds, mark CSV/Markdown summaries, and optionally fail in strict mode, while
the C++ executable remains a focused source of per-run diagnostics.

The provisional toy-suite gate threshold is `max_trial_bound_ratio < 0.8`.
Current runs are well below this value, but the threshold still leaves visible
space below the hard limit of one. Do not wire this gate into routine CTest as a
release-blocking check yet; keep strict mode available for explicit validation
runs until more event-like charged-particle topologies have been tested.

Use a named fixed `event4` topology as the first multi-particle validation
case, rather than introducing a general particle-list parser immediately. The
named topology keeps the result reproducible, inspectable, and easy to run from
both the fixed-grid and shape-unweighting diagnostics. A general direct-particle
input format remains useful later, but it is not needed for the first
event-like gate check.

After the higher-statistics five-case check, the current targeted candidate is
accepted for the next validation stage without retuning. The provisional
`0.8` gate may be used in explicit strict validation commands, but it should
still not be part of routine CTest until at least one more event-like topology
or generator-derived charged-particle sample has been checked.

Use `forward4` as the next named validation topology. It stresses a different
feature from `event4` by placing one charged particle close to the forward edge
of the photon-eta acceptance, while preserving a small charge-conserving final
state that remains easy to inspect. The first gated run does not require
retuning the current targeted candidate.

After the higher-statistics six-case check, `forward4` is accepted as a stable
diagnostic topology for the current targeted candidate. The next validation
step should stress a different feature, preferably higher charged multiplicity
or a clustered jet-like pattern, rather than retuning for the forward case.

Use `multi8` as the first higher-multiplicity validation topology. It keeps the
named-topology approach but doubles the charged multiplicity relative to
`event4` and `forward4`, giving a stronger coherent-current and cell-refinement
stress without introducing a general particle-list input format yet. The first
gated run does not require retuning the current targeted candidate.

After the higher-statistics seven-case check, `multi8` is accepted as a stable
diagnostic topology for the current targeted candidate. The next validation
step should either add a clustered jet-like topology or introduce a minimal
direct-particle text input format to avoid continuing to hardcode every new
diagnostic topology.

Introduce the minimal direct-particle text input format before adding more
built-in validation topologies. The named-topology phase has already covered
symmetric two-pion, asymmetric two-pion, forward, four-particle, and
eight-particle inspectable cases. A line-oriented fixture format is now the
cleaner canonical path for further event-like sampler studies because it keeps
new validation inputs readable, reviewable, and shared between the fixed-grid
truth diagnostic and the shape-unweighting diagnostic.

The first diagnostic text format uses one charged particle per line with
columns `charge inOutSign mass p eta phi`, plus blank lines and `#` comments.
The format is diagnostic-only and does not define HepMC3 output, production
particle selection, recoil, or photon multiplicity. Future clustered or
generator-derived validation cases should prefer checked-in text fixtures over
new hardcoded topology names unless the topology needs special programmatic
construction.

Use a one-cluster-plus-balancing-particles fixture as the first clustered
jet-like validation case. This keeps the topology readable while adding a
different stress from `multi8`: several charged particles are close in
`eta-phi`, but the event remains a small outgoing, charge-conserving pion list.
Incoming-particle clustered fixtures should wait for a separate test of
`inOutSign` semantics rather than being mixed into the first clustered case.

Use a minimal abstract beam-to-final-state fixture as the first direct
incoming/outgoing sign-convention test. The fixture deliberately uses two
incoming `+1` particles and two outgoing `+1` particles, so the ordinary charge
sum is not zero while the signed Low-current sum is zero only if incoming
particles carry `inOutSign = incoming`. This tests the convention that all
four-momenta have positive energy and the incoming minus sign lives only in
`inOutSign`.

Use a tiny HepMC3 fixture with slightly off-axis charged beams as the first
generator-derived sampler-validation case. Exactly beam-axis particles cannot
be represented in the current diagnostic text format, which uses `p eta phi`;
the small transverse components keep the event finite in that format while
still exercising `event.beams()`, final stable selection, MEV-to-GeV momentum
conversion, and PDG-id charge lookup. The adapter-derived selected-particle
text fixture is checked in for reproducible sampler validation, and the
optional HepMC3 build regenerates and compares it against the checked-in file.
This remains diagnostic-only and does not write photons back into HepMC3.

The current eleven-case validation suite is broad enough to move to the first
event-level yield and multiplicity prototype. It covers boosted two-particle,
forward, asymmetric, four-particle, eight-particle, clustered, incoming/outgoing
sign-convention, and HepMC3-derived adapter input cases. The first-prototype
readiness gate for the current diagnostic suite is: all cases pass with
`status=ok`, zero envelope failures, zero fallback uses, and
`max_trial_bound_ratio < 0.8`. Shape L1 and eta/phi projection RMS remain
informative diagnostics only, because fixed-grid references can be resolution
limited for narrow boosted structures.

These readiness criteria are provisional. They justify implementing
`lambda_event` estimation and Poisson multiplicity diagnostics, but they are
not a production event-record guarantee. Before ordinary unweighted photons are
written to HepMC3, LowMC still needs a robust yield-integration diagnostic,
strict event-level behavior for invalid selections or failed envelopes, and a
settled output-writer policy.

Use a separate deterministic midpoint integration routine as the first
`lambda_event` estimator, rather than reusing the weighted proposal samples or
the shape-envelope integral. The integration variables are
`(log k_T, eta, phi)` so the physical measure is included as
`k_T^2 dlog(k_T) deta dphi`. This keeps the event yield independent from the
current shape proposal and gives a transparent baseline for controlled
fixtures. The accompanying half-resolution coarse estimate is a diagnostic
convergence indicator only, not a final uncertainty model.

The first multiplicity diagnostic draws only counts, not photons:
`N_gamma ~ Poisson(lambda_event)` for physical counts and
`N_sample ~ Poisson(oversamplingFactor * lambda_event)` with sample weight
`1 / oversamplingFactor` for diagnostic oversampling. This is still
event-record-output preparation; it does not write HepMC3 photons.

Keep the deterministic midpoint yield estimator as a baseline and regression
diagnostic only. The first yield-validation scan shows good agreement with the
weighted reference for the inspectable `p = 1 GeV` toy case, but large
under-resolution for boosted `p = 10`, `100`, and `1000 GeV` cases and visible
grid sensitivity for the text fixture. The production multiplicity path will
need an adaptive yield-integration strategy before LowMC can attach physical
unweighted photons to event records.

Use deterministic adaptive-cell quadrature as the first adaptive
`lambda_event` prototype. It mirrors the existing grid-adaptive eta-phi
machinery and integrates over the final refined cells rather than relying on a
uniform eta-phi grid. This method is natural for LowMC because the normalization
integral is dominated by the same narrow charged-particle neighborhoods as the
shape-sampling problem.

The current best diagnostic candidate is an adaptive `32 x 32` initial grid
with seven refinement levels, refine fraction `0.1`, and `2 x 2` midpoint
sub-bins per final cell. It agrees with the weighted reference for the current
`p = 1`, `10`, and `100 GeV` two-pion cases and for the first text fixture.
The `p = 1000 GeV` case improves strongly but remains a stress case because the
L5-to-L7 movement is still sizable and the weighted reference has low
effective sample fraction. Do not promote this to a production multiplicity
gate until the high-boost convergence criterion is tightened.

The first tuning gate separates ordinary validation from extreme high-boost
stress validation. For ordinary controlled cases, adaptive L7 is provisionally
acceptable when it differs from the L9 tuning reference by less than `1%` and
the weighted-reference comparison has `ESS/N >= 0.05`. For the `p = 1000 GeV`
stress case, L8 or L9-level checks are required; the current weighted reference
is treated as secondary because its ESS/N remains very low. This lets LowMC
move to a count-plus-shape diagnostic prototype for ordinary cases without
claiming production readiness for the most extreme boost.

The first count-plus-shape path is diagnostic-only. It may combine adaptive
`lambda_event`, Poisson multiplicity sampling, and corrected shape unweighting
for toy and text-fixture inputs, but it must only write diagnostic CSV outputs.
If `N_gamma` or oversampled `N_sample` is larger than the number of accepted
shape candidates, or if the shape envelope fails, the event is reported as
incomplete. This avoids confusing diagnostic photon candidates with production
event-record photons.

Use a small multi-seed count-plus-shape validation as the pre-writer gate. The
current provisional diagnostic gate is: all rows have `status=ok`,
`envelope_ok=1`, `enough_candidates=1`, and `max_trial_bound_ratio < 0.8`.
Use `oversamplingFactor = 100` for this validation so that count fluctuations
produce enough retained candidates for inspection while keeping the analysis
weight explicit. Passing this gate supports developing the HepMC3 writer
contract next; it still does not authorize writing production event-record
photons.

The reduced base-probe compensation study should be presented for independent
review as a diagnostic optimization candidate, not as a default-policy change.
For the thirty bounded Pythia one-event seeds `24681..24710`, the full
`11 x 11` boosted reference took `79.87 s` total, while `5 x 5` reduced probing
with compensation factors `1.30`, `1.35`, and `1.50` all took about `37 s` and
retained the same `174/174` requested photons with zero envelope failures.
Factor `1.30` is the less conservative performance candidate: it stayed below
the matching full-reference trial/bound ratio on every seed, but its narrowest
margin was only `0.000044`. Factor `1.35` should not be recommended because it
overshot the matching full reference at seed `24684`. Factor `1.50` is the
single conservative review candidate because it had no per-seed overshoot and
reduced the maximum trial/bound ratio to `0.321687`, but it also lowered
accepted trial count from `90602` to `79004` for the same `428000` trials.
This tradeoff should remain explicit: larger compensation buys envelope margin
by inflating the proposal envelope, not by changing the Low density. Before any
reduced-probe setting becomes configurable outside diagnostics, it should pass
an independent fresh-seed validation block and an external review of the
safety/efficiency criterion.

The independent fresh-seed block `24711..24740` confirms the conservative
recommendation. Factor `1.30` is no longer a clean safety candidate: it
overshot the matching full `11 x 11` reference on four seeds, although photon
accounting and envelope-failure diagnostics stayed clean. Factor `1.50`
remained below the full reference on every seed, with maximum trial/bound ratio
`0.351192`, while preserving the same requested photons and roughly the same
wall time as factor `1.30`. The accepted decision is therefore to present
factor `1.50` as the only current reduced-probe candidate for independent
review. Factor `1.30` should remain in the comparison table as the more
efficient but less safe performance reference. Reduced probing remains
diagnostic-only.

### Diagnostics and plotting

The core library, core tests, and first diagnostics should not depend on ROOT.

Diagnostic programs should write simple machine-readable outputs such as CSV or
JSON. Plotting should be handled by optional scripts, with Python/matplotlib as
the preferred first path.

ROOT support may be added later as an optional adapter or macro layer for HEP
workflows, but it should not be a dependency of the core.

## Deferred

### Particle selection policy

The default selection is settled, but alternative selections should remain
configurable.

Open options include:

- particles before weak decays
- particles before hadronic decays
- beam remnants
- generator-specific particle selections

### Recoil handling

No recoil prescription is chosen yet. The first model assumes an unchanged hard
scattering and unchanged hadronic final state.

### Multiple photons

Soft-photon exponentiation is a future extension, not part of the first
skeleton.

### Dielectron extension

No implementation path is chosen yet for virtual photons or dielectrons. The
current tendency is to make the first prototype pragmatic and factorized:
produce a soft virtual photon with the Low current, sample its invariant mass,
apply the Kroll-Wada conversion factor, and decay it to an electron-positron
pair. This should be treated as a first-prototype tendency, not as a final
physics decision. A later fully differential implementation should revisit the
full $J_\mu L^{\mu\nu}J_\nu$ tensor structure, especially for lepton-pair
angular distributions, polarization information, and longitudinal virtual
photon contributions.

### Existing QED radiation

The interface to shower photons, PHOTOS, or other QED-radiation sources is
unresolved. A future design must avoid double counting.

PHOTOS is related as a QED-radiation afterburner, but its primary domain is
radiation in decays. LowMC's distinguishing target is an event-wide coherent
current over a configurable set of charged particles.

### Weighted versus unweighted generation

Weighted generation is acceptable initially. The unweighting strategy remains
open.

### Photon phase-space API

The first photon phase-space API should expose physical acceptance limits:
`kTMin`, `kTMax`, `etaMin`, `etaMax`, `phiMin`, and `phiMax`.

This object should describe where photons may be generated, not how they are
sampled internally. Sampler-specific proposal maps and adaptive grids should be
separate implementation details.

### Eta-phi sampling strategy

The fixed-`k_T` `eta-phi` distribution can contain narrow structures around
charged-particle directions. The width depends on particle velocity and on the
mass-over-energy scale.

LowMC should not rely on a fixed grid combined with a fixed safety-factor
rejection envelope. That approach may work for typical charged-particle
momenta, but it is not robust for highly boosted particles.

The diagnostics with boosted and forward charged-pion toy events make this
visible: high momentum and large rapidity can produce structures much narrower
than a practical broad grid.

The preferred future strategy is weighted, proposal-based sampling with an
event-specific proposal map, adaptive refinement, and explicit sampler
diagnostics. Exact `S(k)` evaluation remains required for every sampled photon
point.

The fixed-`k_T` `eta-phi` proposal interface should return sampled `eta`,
sampled `phi`, and the normalized proposal density
`q_eta_phi(eta, phi | k_T, event)` with respect to `d eta d phi`. It should not
return a photon weight and should not be treated as a majorant. The weighted
photon layer combines this density with the one-dimensional `k_T` proposal and
the exact Low density.

The first adaptive implementation is accepted as a diagnostic baseline only:
a mixture of a uniform floor and local uniform patches around charged-particle
directions. It is easy to normalize and inspect, but it is not yet the final
proposal strategy.

A grid-adaptive proposal may also be prototyped as a diagnostic baseline. It
should sample from an event-specific cell map built from `S(k)` evaluations and
retain a uniform floor. It should remain weighted and diagnostic-only until its
reference agreement and performance are better understood.

A smooth local-kernel proposal may also be prototyped as a diagnostic baseline.
The first version should stay analytically normalizable and weighted, with a
uniform floor and exact `S(k)` evaluation at every sampled point. It should be
compared against the tuned grid-adaptive map before any production role is
considered.

### Analytic Low bounds remain audit-only

The exact global and rectangular cell-local Low-current bounds are accepted as
diagnostic infrastructure, including stable finite-mass denominators and
explicit invalid-leg fallback. They are not accepted as the production
rejection envelope.

The predeclared held-out test found no coverage violations, but the median
cell-envelope integral was `8.34` times the current envelope for hadronic Z
events and `40.5` times for SoftQCD. The global bound was billions of times too
large. LowMC therefore keeps the existing sampler and failure semantics. The
cosh-cos inverse-CDF proposal is deferred; it should be reconsidered only with
a materially tighter construction.

### Keep pseudorapidity for charged-direction targeting

True rapidity is the exact center of each massive-leg fixed-`kT` denominator,
but it is not a better default for the present shape-unweighting algorithm.
The grid is already constructed from exact Low intensities; the center
convention affects only optional targeted refinement and boosted-cell probing.

In a thirty-event paired Pythia audit, true rapidity changed the median
hadronic-Z envelope by `+0.03%` and the SoftQCD envelope by `+1.77%`, with no
failures in either convention. It also increased the median SoftQCD probe
count. Keep pseudorapidity as the default and retain true rapidity as an
explicit diagnostic option for future algorithm studies.
