# ALGORITHM

## Inputs

For each event, LowMC needs a list of charged particles:

- four-momentum `p_i`
- charge `Q_i`, in units of the positron charge `e`
- incoming/outgoing sign, denoted `eta_i` in formulas, with `eta_i = -1` for
  incoming particles and `eta_i = +1` for outgoing particles

In code and APIs, this sign should be called `inOutSign`, not `eta`, to avoid
confusion with pseudorapidity.

The default selection is incoming beam particles plus final stable charged
particles. Other selections should be supplied by adapter or user policy.

The photon candidate is represented by a lightlike four-vector `k`.

## Current Evaluation

For a photon momentum `k`, compute the coherent Low current

$$
J^\mu(k) =
\sum_i \eta_i Q_i \frac{p_i^\mu}{p_i \cdot k}
$$

The scalar product uses the metric `(+,-,-,-)`.

The implementation should compute the current first, then form the intensity:

$$
S(k) = -J_\mu(k)J^\mu(k)
$$

This keeps one photon-direction evaluation at `O(N)` in the number of charged
particles.

## Photon Weight

The leading-order density is

$$
\frac{dN_\gamma}{k_T\,dk_T\,d\eta\,d\phi}
=
\frac{\alpha}{(2\pi)^2}S(k)
$$

LowMC treats this as a density with respect to the physical photon measure

$$
d\mu = k_T\,dk_T\,d\eta\,d\phi.
$$

For weighted generation, a sampled photon coordinate

$$
x = (k_T,\eta,\phi)
$$

is drawn from a normalized proposal density

$$
q_\mu(x)
$$

defined with respect to the same measure `dmu`. The corresponding weighted
photon carries

$$
w(x)
=
\frac{
\frac{\alpha}{(2\pi)^2}S_\mathrm{event}(k)
}{
q_\mu(x)
}.
$$

Here `Photon::weight` is an event-level contribution to the expected photon
yield in the configured phase-space region. Histograms filled with this weight,
and then averaged over input generator events, estimate inclusive Low photon
spectra.

If an implementation samples `k_T` using an ordinary density `r(k_T)` with
respect to `dk_T`, rather than with respect to `k_T dk_T`, then the density with
respect to the physical measure is

$$
q_{k_T,\mu}(k_T)=\frac{r(k_T)}{k_T}.
$$

Equivalently, the event weight may be written with an explicit factor of `k_T`
in the numerator when proposal densities are expressed with respect to
`dk_T deta dphi`. The implementation should document which convention an API
uses and convert internally to the physical-measure convention above.

The first implementation may remain weighted.

Efficient proposal distributions are part of the algorithmic design. At fixed
`k_T`, the density can be sharply structured near charged-particle directions,
so the first implementation should keep the sampling strategy explicit and
testable.

The sampling strategy should be understood as an application of established
adaptive Monte Carlo ideas to the Low-photon problem, not as a new generic
sampling method. Weighted importance sampling, adaptive cell maps, variance
diagnostics, and later accept/reject unweighting are standard tools in Monte
Carlo event generation. LowMC's algorithmic task is to choose and validate such
tools for the event-specific `eta-phi-k_T` structures of the coherent Low
current.

A proposal density is not a majorant. It must have support wherever the target
density can contribute, but it does not need to bound the target from above. A
majorant or envelope is only required for rejection sampling or unweighted
generation, where one needs a validated bound on `w(x)` or on the target density
relative to the proposal.

The intended first proposal factorization is

$$
q_\mu(k_T,\eta,\phi)
=
q_{k_T,\mu}(k_T)
\,
q_{\eta\phi}(\eta,\phi \mid k_T,\mathrm{event}),
$$

where `q_kT,mu` is normalized with respect to `k_T dk_T`, and `q_eta_phi` is
normalized with respect to `deta dphi` at fixed `k_T`. This keeps the relatively
simple one-dimensional `k_T` problem separate from the event-dependent
two-dimensional `eta-phi` sampling problem.

## kT Proposal

The default one-dimensional `k_T` proposal should be logarithmic in `k_T`.

At fixed photon direction, the Low current scales as

$$
J^\mu(k) \propto \frac{1}{k_T},
$$

and therefore

$$
S(k) = -J^2 \propto \frac{1}{k_T^2}.
$$

Since the physical measure contains `k_T dk_T`, the generated yield behaves as

$$
dN_\gamma
\propto
k_T\,dk_T\,S(k)
\propto
\frac{dk_T}{k_T}.
$$

The natural default is therefore uniform sampling in

$$
u = \ln k_T
$$

over the configured range `[kTMin, kTMax]`.

With respect to the ordinary measure `dk_T`, the normalized proposal density is

$$
r(k_T)
=
\frac{1}{k_T \ln(kTMax/kTMin)}.
$$

With respect to the physical `k_T dk_T` measure used by the weighted-photon
convention, the corresponding one-dimensional proposal density is

$$
q_{k_T,\mu}(k_T)
=
\frac{r(k_T)}{k_T}
=
\frac{1}{k_T^2 \ln(kTMax/kTMin)}.
$$

This is the quantity that enters the factorized proposal density

$$
q_\mu(k_T,\eta,\phi)
=
q_{k_T,\mu}(k_T)
\,
q_{\eta\phi}(\eta,\phi \mid k_T,\mathrm{event}).
$$

The logarithmic `k_T` proposal is a default importance-sampling choice, not a
physics cutoff. The physical `k_T` range is still defined by `PhotonPhaseSpace`
and by the validity of the soft approximation.

## Photon Phase Space

The first public phase-space concept should be a simple physical acceptance
rather than a sampler-specific object.

A tentative `PhotonPhaseSpace` should contain:

- `kTMin`
- `kTMax`
- `etaMin`
- `etaMax`
- `phiMin`
- `phiMax`

The user-facing API should define the physical region in which photons may be
generated. The sampler is then free to choose an internal proposal distribution
inside that region.

The spectrum is differential in

$$
k_T\,dk_T\,d\eta\,d\phi,
$$

so photon weights must specify both the physical phase-space region and the
proposal density used to sample it.

The current diagnostics suggest that a broad fixed `eta-phi` grid is useful for
visualization but not sufficient as the main sampling strategy. At high momentum
or forward rapidity, structures near charged-particle directions can be much
narrower than a practical global grid cell.

The first sampler should therefore be weighted and proposal-based. Adaptive
proposal details remain open, but the phase-space API should not assume a
particular sampler.

The first adaptive fixed-`k_T` proposal is deliberately simple: a mixture of a
uniform floor over the full `eta-phi` acceptance and local uniform patches
around charged-particle directions. The local patch width is scaled from
`m / E` and clamped between configurable minimum and maximum half-widths. This
proposal is not a majorant and does not replace exact `S(k)` evaluation; it is
only an importance-sampling density with support everywhere through the uniform
floor.

A first grid-adaptive diagnostic proposal is also useful for learning. It
evaluates the actual event-level `S(k)` at cell centers on a coarse `eta-phi`
grid, refines high-mass cells, samples cells proportional to their estimated
mass, and samples uniformly inside the chosen cell. This is still a weighted
proposal, not a majorant: exact `S(k)` is evaluated again at the sampled point,
and a uniform floor keeps support nonzero over the full acceptance.

The grid-adaptive diagnostic can also choose refinement candidates with a
variance-aware score. The first such score probes a small lattice of points
inside each candidate cell, estimates the variance of `S(k)` inside that cell,
and uses `area^2 Var[S(k)]` as the refinement score. The proposal density after
refinement is still built from estimated cell mass, and exact Low weights are
still evaluated at every sampled point. This keeps the estimator weighted and
unbiased while testing a refinement rule closer to adaptive Monte Carlo
variance reduction.

For the first full phase-space grid diagnostic, the proposal factorizes into
the existing log-uniform `k_T` proposal and a grid-adaptive `eta-phi` proposal.
The grid is built once per event at the geometric mean of the configured `k_T`
range. This is justified as a diagnostic because the fixed-direction Low
current has the expected `1 / k_T` scaling, so the angular shape of `S(k)` is
independent of `k_T` up to an overall normalization. The sampled photon is
still assigned the full proposal density
`q_mu = q_kT,mu q_eta_phi`, and exact Low weights are evaluated at the sampled
`k_T`, `eta`, and `phi`.

The first event-level weighted API keeps generator policy outside the core. The
caller passes already selected `ChargedParticle` objects and a
`PhotonPhaseSpace`; the core builds a `LowCurrent`, constructs the tuned
grid-adaptive angular proposal, samples the requested number of weighted
photons with an explicit random seed, and returns the weighted photons plus
basic diagnostics. HepMC/Pythia record interpretation, particle-selection
policy, recoil, unweighting, and multiple-photon semantics remain outside this
minimal API.

The event-level API also supports an integer `oversamplingFactor` for weighted
histogramming. The returned sample count is
`samples * oversamplingFactor`, and each returned photon weight is divided by
`oversamplingFactor`. This increases histogram statistics without changing the
expected event contribution. It is a weighted-analysis convenience, not a
physical enhancement of the Low photon yield and not a production HepMC3
unweighting convention.

A first unweighting diagnostic now probes whether weighted samples can be
converted into unweighted photon candidates by accept/reject sampling. It uses
one calibration pass to estimate empirical maximum weights, multiplies those
maxima by a safety factor to define an envelope, and then tests the envelope on
a fresh trial pass. Every trial photon still receives an exact Low weight from
the weighted generator before the accept/reject decision.

Three diagnostic envelope modes are available:

- `global`: one empirical maximum weight for the whole event and phase-space
  region,
- `cell`: one empirical maximum weight for each refined grid-adaptive
  `eta-phi` cell, with a global fallback for cells not populated in the
  calibration pass,
- `probed-cell`: one maximum weight for each refined grid-adaptive `eta-phi`
  cell from an explicit lattice of exact `w = f/q` evaluations inside the
  cell. If random calibration finds a larger weight in a cell than the probe
  lattice, the diagnostic uses the larger calibration value for that cell and
  counts the cell as suspect.

The diagnostic reports the calibration maximum, global envelope weight, trial
maximum, global and local trial/calibration ratios, accepted photon count,
observed and expected acceptance efficiency, fallback count for uncalibrated
cells, suspect-cell count, and whether any trial sample exceeded its envelope.

This is only a diagnostic for the shape-unweighting step. An accepted photon in
this diagnostic does not yet define the physical number of photons to write to
an event record. Production HepMC3 output still requires a separate multiplicity
and normalization convention, and any envelope failures indicate that a
stronger event-specific bounding pass is needed.

The first small global-envelope scan uses `lowmc_unweighting_diagnostic` over
transverse two-pion toy events with momenta `p = 1`, `10`, and `100 GeV`, three
seeds, and safety factors `1.5`, `3`, and `10`. With `1000` calibration and
`1000` trial samples per point, the global envelope held for all `p = 1` and
`10 GeV` runs. At `p = 100 GeV`, safety factor `1.5` failed in one of three
runs, with two trial samples above the envelope; safety factors `3` and `10`
held in this small scan. Acceptance efficiency decreased strongly with safety
factor and with momentum: at `p = 100 GeV`, the mean observed efficiency was
about `1.8%`, `1.0%`, and `0.27%` for safety factors `1.5`, `3`, and `10`.

The first per-cell scan uses the same `p = 100 GeV` toy benchmark. With only
`1000` calibration samples, the per-cell envelope is much more efficient but
not reliable: many refined cells are empty, safety factor `10` still has
failures in two of three seeds, and the largest local trial/calibration ratio
is about `33`. With `10000` calibration samples and `1000` trial samples,
the per-cell envelope becomes much more promising. Safety factors `3` and `10`
hold in the three-seed scan; at safety factor `3`, the mean observed
acceptance efficiency is about `22%`, compared with about `0.5%` for the
global envelope in the same scan. The remaining warning sign is that roughly
`314` of `2257` cells are still empty on average, and trial samples falling in
those cells use the global fallback.

The first explicit probing pass also uses the `p = 100 GeV` benchmark. It
probes each final refined cell on a regular lattice and uses the largest exact
weight found in that cell as the local pre-safety bound. A `3 x 3` lattice is
not enough: the largest local trial/probed ratio is about `8.1`, so safety
factor `3` still fails. A `7 x 7` lattice is much better in this benchmark:
the largest local ratio drops to about `1.94`, and safety factor `3` holds in
the three-seed `1000 + 1000` trial scan when combined with the suspect-cell
hybrid fallback. The observed acceptance at safety factor `3` is about `21%`,
compared with about `1%` for the global envelope and about `14%` for the
sample-calibrated per-cell envelope in the same scan.

An edge-inclusive `7 x 7` probe rule has also been checked. It reduces the
number of suspect cells in the `p = 100 GeV` toy, but it does not improve the
worst trial/bound ratio or the acceptance efficiency in the current scan. With
`10000` trial samples per seed at safety factor `3`, both interior and
edge-inclusive probing hold for three seeds; the interior rule remains slightly
more efficient in this benchmark.

This suggests that local envelopes tied to the adaptive cell map are a
promising production direction, and explicit per-cell probing with a
suspect-cell hybrid fallback is already a clear improvement over sparse random
per-cell calibration. It is still not a certified production majorant: the
next step should test more seeds, forward/large-rapidity charged particles, and
more extreme momenta before treating it as the production shape-unweighting
path.

The first broader stress scan keeps the `7 x 7` interior `probed-cell` bound
with safety factor `3` and uses five seeds with `1000` calibration samples and
`2000` trial samples per case. It holds for the transverse `p = 100 GeV` toy
and for a forward `eta = 4`, `p = 100 GeV` toy. The forward case is not the
first failure mode in this scan; it has no envelope failures and an observed
acceptance near `30%`. The first failure is instead the deliberately extreme
transverse `p = 1000 GeV` toy, where three of five runs fail and the maximum
trial/bound ratio is about `7.1`. Increasing the safety factor to `10` holds
for this `p = 1000 GeV` scan but lowers the observed acceptance to about `6%`.

The current interpretation is therefore:

- `7 x 7` probed-cell bounds with safety factor `3` remain the leading
  candidate for typical and high-energy diagnostic cases tested so far.
- Very extreme boosted cases expose the need for either a larger safety factor,
  finer/adaptive probing, or recursive treatment of suspect cells.
- Production unweighting should not rely on one global safety factor until this
  boosted-tail behavior is better controlled.

The first targeted boosted-tail improvement is simply a denser per-cell probe
lattice for the same hardened `probed-cell` bound. For the transverse
`p = 1000 GeV` toy, increasing the interior lattice from `7 x 7` to `11 x 11`
removes the failures in the same five-seed, `1000 + 2000` sample scan at safety
factor `3`. The observed acceptance remains near `20%`, much better than the
about `6%` obtained by keeping `7 x 7` probes and raising the safety factor to
`10`. Running the same `11 x 11` setting on the transverse `p = 100 GeV` and
forward `eta = 4`, `p = 100 GeV` regression cases also gives no failures and
preserves the previous acceptance levels within the diagnostic precision.

This points to finer local probing, rather than a blunt safety-factor increase,
as the better next production candidate. The open design question is how to
choose the probe density event by event: a fixed `11 x 11` lattice is simple
but more expensive than `7 x 7`, while a future adaptive policy could increase
probe density only for highly boosted or otherwise suspect cells.

The recommended first production-like shape-unweighting policy is therefore:

- use the tuned grid-adaptive proposal map as the proposal,
- build a `probed-cell` envelope on the same final refined cells,
- use `11 x 11` interior probes per cell,
- use the suspect-cell hybrid fallback
  `cellBound = max(probedCellMax, calibrationCellMax)`,
- use safety factor `3` as the first diagnostic default,
- return envelope diagnostics with every event-level unweighting attempt.

For the current tuned grid in the two-pion diagnostics, the final map has
`2257` cells. The probe cost is therefore:

| probe lattice | exact probe evaluations |
| ---: | ---: |
| `7 x 7` | `110593` |
| `11 x 11` | `273097` |

These evaluations are only for building the local envelope; every proposed
photon still receives its own exact Low weight before the accept/reject
decision. A fixed `11 x 11` policy is deliberately simple and robust for the
tested benchmarks, but it is not necessarily the final performance optimum.
The likely later improvement is an adaptive policy that starts from a cheaper
lattice and refines cells using boost, cell size, local weight spread, or
suspect-cell diagnostics.

The first public unweighting configuration should expose only a small number of
controls: safety factor, probe bins, probe rule, and a strict/diagnostic
failure policy. If any accepted event-level path observes `w > M`, it must
report the failure count and maximum violation ratio. Production event-record
output should either reject that event or fall back to weighted diagnostic
output; it should not silently write failed-envelope photons as ordinary
unweighted final-state particles.

## Shape-Unweighting API

The first generator-independent shape-unweighting API should expose photon
candidates with unweighted shape, but it should not decide the physical photon
multiplicity and should not write event records.

A tentative configuration object should contain:

- `PhotonPhaseSpace phaseSpace`
- `int trials`
- `uint64_t seed`
- `GridAdaptiveEtaPhiProposalConfig gridConfig`
- `int calibrationTrials`, default `1000`
- `int probeBins`, default `11`
- `ProbeRule probeRule`, default `Interior`
- `double safetyFactor`, default `3`
- `EnvelopeFailurePolicy failurePolicy`

The first `EnvelopeFailurePolicy` values should be:

- `Strict`: no physical unweighted candidates are returned if any envelope
  violation is observed,
- `Diagnostic`: return diagnostics and accepted candidates, but mark the result
  as not suitable for event-record output.

A weighted fallback policy remains deferred. If added later, it should be
explicit and diagnostic until the downstream output convention is settled.

The result object should contain:

- accepted shape-unweighted photon candidates,
- the sampled `k_T`, `eta`, and `phi` for accepted candidates,
- charge-conservation diagnostics,
- envelope diagnostics,
- proposal/grid-build diagnostics,
- a status enum such as `Ok`, `EnvelopeFailed`, or `NoTrials`.

Accepted shape-unweighted candidates should have ordinary photon four-momenta
and no per-photon analysis weight. They are not yet physical event-record
photons: a later multiplicity layer decides how many such candidates to request
or keep for an event. Rejected trial samples should not be returned by default,
but the diagnostics must summarize them well enough for validation.

Envelope diagnostics should include at least:

- trial count,
- accepted count,
- observed and expected acceptance efficiency,
- failed-envelope count,
- maximum trial/bound ratio,
- probe evaluations,
- grid cell count,
- suspect cell count,
- fallback count,
- safety factor,
- probe bins and probe rule.

The API should compute exact Low weights for all trial photons before
accept/reject. For `probed-cell` bounds, the local bound is

$$
M_c =
\mathrm{safetyFactor}
\times
\max(w^\mathrm{probe}_c, w^\mathrm{calib}_c),
$$

where `c` is the refined proposal cell. A trial in cell `c` is accepted with
probability `min(1, w / M_c)`.

For event-record production, the default failure policy should eventually be
`Strict`: an envelope failure means ordinary unweighted photons from that event
are not trustworthy. During diagnostics and development, `Diagnostic` is more
convenient because it preserves the evidence needed to improve the envelope.

`lowmc_shape_unweighting_diagnostic` is now the preferred executable for
stress-testing this API. It runs the generator-independent
`generateShapeUnweightedPhotonCandidates(...)` path directly and reports the
same quantities needed for envelope validation: status, accepted candidates,
observed and expected acceptance efficiency, failed-envelope count,
maximum trial/bound ratio, envelope integral, suspect-cell count, fallback
count, and probe evaluations. The older standalone
`lowmc_unweighting_diagnostic` remains useful as an internal
envelope-development comparison tool, but future production-like stress studies
should exercise the public shape-unweighting API.

The API diagnostic reports `max_trial_bound_ratio = max(w / M_c)` after the
safety factor has already been applied. A value below one means the tested
trials stayed inside the current envelope. This differs from the older
standalone diagnostic's local trial/probe ratio, which was closer to a
pre-safety diagnostic for the probed cell maximum.

The current API sampler uses the corrected cell-local envelope construction:
after per-cell bounds `M_c` are built, it samples cells with probability
proportional to `M_c` times the original proposal probability in that cell,
samples uniformly inside the selected cell, and accepts with probability
`w / M_c`. This makes the accepted density proportional to `f`, while the
reported envelope integral is the normalizing constant of the proposal
envelope. The older diagnostic samples from the original proposal and is now
only a bound/efficiency development comparison, not the production shape
sampler.

The corrected sampler is now validated by a first eta-phi shape-comparison
script. It generates accepted candidates with the API sampler, bins their
normalized `eta-phi` distribution, and compares it with a fixed-grid Low
intensity reference through an L1 shape distance and one-dimensional projection
RMS values. These metrics are diagnostic only; release thresholds remain open.

The first boosted-tail tuning scan for the corrected sampler shows that the
transverse `p = 1000 GeV` toy is substantially harder than it looked in the
older original-proposal diagnostic. Fixed `21 x 21` interior probing with
safety factor `20` holds in a five-seed, `1000 + 2000` trial scan, but the
observed acceptance is only about `0.14%`. This is a useful conservative
fallback baseline, not an attractive production default.

A first adaptive boosted-cell probe policy now exists as a prototype. It first
probes all cells with the base lattice, runs the usual random calibration, and
then promotes selected cells to a denser lattice and a larger local safety
factor before constructing the corrected envelope. The current promotion rule
is deliberately simple and deterministic:

- promote cells whose calibration maximum exceeds the base probe maximum,
- promote cells near selected charged-particle directions.

The corrected sampler is otherwise unchanged: cells are still drawn
proportional to `q_c M_c`, sampled uniformly inside the chosen cell, and
accepted with `w / M_c`.

In the first transverse `p = 1000 GeV` scan, base `11 x 11` probing with
promoted `21 x 21` cells and local safety factor `15` holds across five seeds.
It uses about `291000` exact probe evaluations, compared with about `995000`
for probing every cell with `21 x 21`. The observed acceptance is still only
about `0.17%`, so the policy is a useful cost reduction and diagnostic
prototype, not yet a satisfactory production unweighting solution.

The first refinement of this prototype exposes the angular grid itself as a
control. Narrowing the charged-particle promotion neighborhood from one-cell
halo to only cells containing the charged-particle direction reduces the number
of promoted cells, but it does not significantly reduce the envelope integral:
the same few high-bound cells dominate. A more useful improvement is to build a
finer grid-adaptive proposal before constructing the envelope. In the transverse
`p = 1000 GeV` scan, a `48 x 48` initial grid with four refinement levels,
base `11 x 11` probes, promoted `21 x 21` probes, local safety factor `15`, and
zero charged-neighbor halo holds across five seeds. It lowers the mean envelope
integral from about `315` to about `248` and raises the expected acceptance from
about `0.12%` to about `0.13%`, while increasing the exact probe count from
about `281000` to about `806000`. This is a real but modest efficiency gain; it
turns grid refinement, not denser promoted probes alone, into the next
algorithmic lever.

The next refinement applies that lever locally. Before probing the envelope,
LowMC can split only the grid cells that contain selected charged-particle
directions. This targeted charged-direction refinement keeps the corrected
cell-local sampler unchanged: the final cells are still drawn proportional to
`q_c M_c` and accepted with `w / M_c`. It only changes the cell map on which
the local bounds are built.

In the transverse `p = 1000 GeV` stress scan, four targeted split levels on
the default `32 x 32`, level-3 grid are better than the global `48 x 48`,
level-4 refinement. With base `11 x 11` probes, promoted `21 x 21` probes,
zero charged-neighbor halo, and promoted local safety factor `5`, the five-seed
scan has no envelope failures, mean envelope integral about `8.0`, maximum
trial/bound ratio about `0.36`, observed acceptance about `4.6%`, about `2353`
grid cells, and about `292000` exact probe evaluations. The older global
`48 x 48`, level-4 candidate used about `6588` cells and `806000` probe
evaluations while still giving only about `0.16%` observed acceptance in the
same short scan. Targeted charged-direction refinement is therefore the leading
diagnostic candidate, but it still needs higher-statistics and broader-topology
validation before it becomes a production default.

A first multi-seed validation supports using this targeted setting as the
default diagnostic candidate. In a three-seed, `30000`-trial-per-seed
comparison over transverse `p = 1 GeV`, transverse `p = 1000 GeV`, and forward
`eta = 4`, `p = 100 GeV` two-pion cases, it has zero envelope failures in all
cases. For the `p = 1000 GeV` stress case it accepts `4394` candidates out of
`90000` trials, compared with `127` for the previous `32 x 32`, level-3
boosted-cell baseline and `164` for the global `48 x 48`, level-4 baseline.
A focused five-seed, `100000`-trial-per-seed `p = 1000 GeV` check also has
zero failures, `24270` accepted candidates, mean envelope integral about `8.0`,
and maximum trial/bound ratio about `0.365`. These are still diagnostic
statistics, not a mathematical proof of a global majorant.

The diagnostic workflow now promotes this candidate through named script
presets rather than by hiding the lower-level executable options.
`scripts/validate_shape_unweighting_candidates.py` defaults to the targeted
candidate set, while `--candidate-set comparison` keeps the previous
boosted-cell and global refined-grid baselines available. The
`lowmc_shape_unweighting_diagnostic` executable itself keeps explicit tuning
options so individual algorithm choices remain visible in command lines and
logs.

The same workflow includes a first less-symmetric two-pion topology: a
positive pion with `p = 100 GeV`, `eta = 0.35`, `phi = 1.1` and a negative
pion with `p = 35 GeV`, `eta = -0.75`, `phi = 4.0`. This is still
charge-neutral and easy to inspect, but it is no longer the symmetric
back-to-back benchmark. In the first targeted three-seed validation with
`30000` trials per seed, it has zero envelope failures, about `9.7%` observed
acceptance, mean envelope integral about `2.8`, and maximum trial/bound ratio
about `0.65`. The case should be treated as an early topology-dependence check;
it does not yet prove that the targeted candidate is optimal for all event
shapes.

The first topology-dependence scan shows that two independent refinements help
different cases. Adding a small targeted split halo improves the asymmetric
case, while adding one more targeted split level improves the symmetric
`p = 1000 GeV` stress case. Their combination,
`targeted_l5_halo0.5_p11_to_p21_s3_local5`, improves both without a large cost
increase. In a three-seed, `30000`-trial-per-seed scan, it gives zero failures,
about `11.4%` observed acceptance for the symmetric `p = 1000 GeV` case and
about `20.9%` for the asymmetric case, with about `2377` cells and about
`294000` probe evaluations. A focused five-seed, `100000`-trial-per-seed check
also has zero failures in both cases, with maximum trial/bound ratios below
`0.39`. This combined setting is therefore the preferred diagnostic
shape-unweighting candidate for the next validation round.

The first broader high-statistics diagnostic pass keeps this candidate
failure-free in all four current toy cases. With five seeds and `100000` trials
per seed, transverse `p = 1 GeV`, transverse `p = 1000 GeV`, forward
`eta = 4`, `p = 100 GeV`, and the asymmetric two-pion case all have zero
envelope failures, zero fallback uses, and maximum trial/bound ratio below
`0.39`. This supports provisional envelope-validation gates for toy
diagnostics: all configured validation cases should have no envelope failures,
no fallback uses, and a comfortable post-safety maximum trial/bound margin
below one. Shape L1 and projection RMS values should still be reported, but
they should not yet be release-blocking because the fixed-grid truth references
can be grid-resolution limited for narrow boosted structures.

The candidate-validation script now has an optional envelope gate layer for
these toy diagnostics. `--gate-policy diagnostic` records pass/fail columns in
the CSV rows and Markdown summary without changing the process exit status.
`--gate-policy strict` uses the same checks but returns nonzero if any
candidate-case summary fails. The maximum trial/bound threshold is explicit via
`--gate-max-trial-bound-ratio` and must be below one when the gate is enabled.
This keeps exploratory scans unconstrained by default while making the
provisional toy gates executable when desired.

A first gated toy-suite run with three seeds, `30000` trials per seed, and
`--gate-max-trial-bound-ratio 0.8` passes all four current toy cases. The worst
case is still below `0.38`, so a looser near-one threshold such as `0.95` would
not discriminate among the current diagnostics. For now, `0.8` is the
provisional diagnostic threshold: it is well above the observed toy-suite
maximum while still encoding a real margin below the mathematical limit of one.
This threshold should remain diagnostic until more event-like topologies are
covered.

The first event-like direct-input topology is the named `event4` case. It is a
fixed, inspectable, charge-conserving four-pion final state:

| particle | charge | p [GeV] | eta | phi |
| :--- | ---: | ---: | ---: | ---: |
| pi 1 | `+1` | `80` | `0.20` | `0.70` |
| pi 2 | `-1` | `55` | `-0.80` | `2.40` |
| pi 3 | `+1` | `25` | `1.60` | `4.30` |
| pi 4 | `-1` | `12` | `-2.20` | `5.40` |

Both the fixed-grid truth diagnostic and the shape-unweighting diagnostic
accept it as `--topology event4`. The candidate-validation workflow includes it
as `event4_inspectable`, alongside the earlier two-pion cases. This is still a
controlled topology, not a generator-event sample, but it exercises more than
two charged directions and avoids relying only on back-to-back structure.

The higher-statistics expanded-suite check uses all five current cases with
five seeds and `100000` trials per seed. The preferred targeted candidate
passes the provisional `0.8` diagnostic gate in every case. The worst
maximum trial/bound ratio is about `0.392` in `event4_inspectable`; all cases
have zero envelope failures and zero fallback uses. This supports keeping the
current candidate and threshold for the next validation step, while still
treating the gate as diagnostic rather than release-blocking.

The second event-like direct-input topology is the named `forward4` case. It
keeps the same fixed, inspectable four-pion structure, but stresses a very
forward charged direction:

| particle | charge | p [GeV] | eta | phi |
| :--- | ---: | ---: | ---: | ---: |
| pi 1 | `+1` | `350` | `4.50` | `0.40` |
| pi 2 | `-1` | `120` | `2.80` | `1.35` |
| pi 3 | `+1` | `45` | `-0.60` | `3.30` |
| pi 4 | `-1` | `30` | `-2.00` | `5.10` |

It is exposed as `--topology forward4` in both diagnostics and as
`forward4_inspectable` in the candidate-validation workflow. In the first
three-seed, `30000`-trial-per-seed gated run over all six current cases,
`forward4_inspectable` passes the provisional `0.8` gate with zero envelope
failures, zero fallback uses, about `22.9%` observed acceptance, and maximum
trial/bound ratio about `0.372`.

The higher-statistics six-case check uses five seeds and `100000` trials per
seed. All `30` candidate-case-seed rows pass the provisional `0.8` diagnostic
gate, again with zero envelope failures and zero fallback uses. The largest
maximum trial/bound ratio is about `0.396` in `forward4_inspectable`; the
earlier `event4_inspectable` case remains close, at about `0.392`. This keeps
the current targeted candidate viable for the next validation stage without
retuning.

The third event-like direct-input topology is the named `multi8` case. It
stresses higher charged multiplicity while keeping the event deterministic and
inspectable:

| particle | charge | p [GeV] | eta | phi |
| :--- | ---: | ---: | ---: | ---: |
| pi 1 | `+1` | `140` | `1.20` | `0.25` |
| pi 2 | `-1` | `95` | `-1.10` | `0.95` |
| pi 3 | `+1` | `70` | `2.40` | `1.80` |
| pi 4 | `-1` | `50` | `-2.70` | `2.55` |
| pi 5 | `+1` | `35` | `0.30` | `3.35` |
| pi 6 | `-1` | `25` | `-0.40` | `4.10` |
| pi 7 | `+1` | `18` | `3.10` | `5.00` |
| pi 8 | `-1` | `12` | `-3.20` | `5.75` |

It is exposed as `--topology multi8` in both diagnostics and as
`multi8_inspectable` in the candidate-validation workflow. In the first
three-seed, `30000`-trial-per-seed gated run over all seven current cases,
`multi8_inspectable` passes the provisional `0.8` gate with zero envelope
failures, zero fallback uses, about `24.0%` observed acceptance, and maximum
trial/bound ratio about `0.382`. The mean envelope integral is about `5.89`,
larger than the four-pion topologies as expected for a higher-multiplicity
current.

The higher-statistics seven-case check uses five seeds and `100000` trials per
seed. All `35` candidate-case-seed rows pass the provisional `0.8` diagnostic
gate, with zero envelope failures and zero fallback uses. The `multi8` case
remains stable: maximum trial/bound ratio about `0.384`, observed acceptance
about `24.1%`, and mean envelope integral about `5.89`. The suite-wide maximum
is still the `forward4_inspectable` value of about `0.396`.

## Road To Unweighted Event Output

The long-term event-record goal is to add ordinary unweighted Low photons to an
input generator event:

```text
input HepMC3 event
-> selected charged particles
-> weighted Low photon proposal samples
-> shape unweighting
-> photon multiplicity/normalization
-> output HepMC3 event with added stable photons
```

The current weighted sampler provides unbiased weighted samples from the Low
distribution. These samples are suitable for histograms, especially with
`oversamplingFactor`, but they are not yet ordinary event-record photons.

To turn the weighted proposal samples into unweighted photon directions and
momenta, LowMC needs an envelope for the importance weight

$$
w(x) = \frac{f(x)}{q(x)},
$$

where `x = (k_T, eta, phi)`, `f` is the Low photon density in the configured
phase space, and `q` is the proposal density. If `M >= max_x w(x)`, a proposal
sample can be accepted with probability `w(x) / M`. The accepted samples then
have the desired Low shape without per-photon analysis weights. This statement
is exact for one global envelope constant `M`.

For cell-local envelopes `M_c`, the production sampler must be formulated with
the correct envelope density. Sampling from the original proposal `q` and
accepting with `w / M_c` would produce a density proportional to `f / M_c`, not
to `f`, whenever `M_c` varies between cells. LowMC's corrected cell-local
construction therefore samples cells from the normalized envelope density
proportional to `q_c M_c`, where `q_c` is the original proposal probability in
cell `c`, and then accepts with probability `w / M_c`. This preserves the target
shape `f / lambda_event` when the envelope is valid.

For speed, the current implementation keeps this corrected construction but
avoids repeated full scans over the refined cell list. Envelope-cell sampling
uses the same cumulative-mass rule as the linear reference, implemented with a
prefix-sum table and binary search. Proposal-density lookup for sampled
`eta-phi` points uses an eta-bucket cell index that returns the same first
matching cell as the linear reference, including boundary-adjacent test points.
This is intended as a bookkeeping optimization only: it must not change the
proposal density, envelope mass, acceptance probability, or accepted Low shape.
Focused tests compare the indexed path against the linear reference before the
optimization is trusted.

A first four-seed Pythia-produced timing scan after this change keeps the same
photon counts and envelope diagnostics as the pre-index checks while reducing
total event time by roughly a factor of four for the tested seeds. The
bookkeeping index is therefore part of the default shape-unweighting path. It
should still be regarded as a semantics-preserving implementation detail: any
future changes to cell lookup or envelope-cell sampling must continue to be
tested against the linear reference.

Shape unweighting is separate from event multiplicity. After the photon shape
can be sampled without weights, LowMC still needs a normalization rule for how
many photons to add to a given input event in the configured phase space. The
event yield is

$$
\lambda_\mathrm{event}
= \int_{\mathrm{phase\ space}} k_T\,dk_T\,d\eta\,d\phi\,
\frac{\alpha}{(2\pi)^2}S_\mathrm{event}(k).
$$

The first physical multiplicity convention is

$$
N_\gamma \sim \mathrm{Poisson}(\lambda_\mathrm{event}).
$$

Conditional on `N_gamma`, photon momenta should be sampled independently from
the normalized shape density

$$
p_\mathrm{shape}(k) =
\frac{
\frac{\alpha}{(2\pi)^2}S_\mathrm{event}(k)
}{
\lambda_\mathrm{event}
}.
$$

This is the natural first event-record convention, but it should be documented
as an independent-emission model used with the leading-soft density, not as a
full multi-photon exponentiation or recoil treatment.

The integrated yield estimate used for multiplicity must be recorded with
diagnostics before ordinary photons are written. At minimum, the event-record
path should report the configured phase space, the estimated
`lambda_event`, its numerical uncertainty or deterministic integration status,
the shape-envelope status, the maximum trial/bound ratio, and the random seeds
or streams used for multiplicity and shape sampling.

The first implemented `lambda_event` estimator is a deterministic midpoint
integral in variables

$$
(\log k_T,\eta,\phi).
$$

For each midpoint it evaluates

$$
\frac{\alpha}{(2\pi)^2}S_\mathrm{event}(k)
$$

and multiplies by the measure Jacobian

$$
k_T\,dk_T = k_T^2\,d\log k_T.
$$

The diagnostic reports both the requested-grid estimate and a simple
half-resolution coarse estimate. The difference is only an early convergence
indicator; it is not yet a production uncertainty model.

The first validation scan compares this deterministic estimate with the mean
weight from the existing full weighted sampler. It shows that a uniform
midpoint grid can work for the broad `p = 1 GeV` toy case, but under-resolves
the narrow boosted structures already by `p = 10 GeV` and remains inadequate
for `p = 100` and `1000 GeV` at the tested grids. This is consistent with the
eta-phi sampling problem: the yield integral is dominated by the same narrow
charged-particle neighborhoods that made naive uniform sampling inefficient.
The next normalization step should therefore use an adaptive yield integration
strategy rather than simply increasing the uniform grid until it becomes
prohibitively expensive.

The first adaptive yield prototype reuses the grid-adaptive eta-phi cell map.
For each midpoint in `log k_T`, LowMC builds the same kind of refined eta-phi
cell partition used by the grid-adaptive sampler, then performs deterministic
midpoint quadrature inside each final cell. This keeps the estimator
deterministic and transparent while concentrating angular resolution near the
charged-particle directions.

The prototype method is reported as
`adaptive-log-kt-grid-eta-phi-midpoint`. In the first validation scan, an
adaptive `32 x 32` initial grid with seven refinement levels and `2 x 2`
sub-bins per final cell agrees with the weighted reference for `p = 1`, `10`,
and `100 GeV` two-pion cases and for the first text fixture. The
`p = 1000 GeV` stress case is much improved but still shows sizable movement
between refinement levels, so it remains a convergence stress test rather than
a settled production gate.

For the leading-Low real-photon model, the angular shape is independent of
`k_T` after multiplying by the physical-measure factor `k_T^2`. Therefore the
same adaptive angular yield can be evaluated once at a representative
`k_T`, chosen as the geometric mean of `kTMin` and `kTMax`, and then multiplied
by `ln(kTMax / kTMin)`. LowMC implements this as
`factorized-adaptive-eta-phi-midpoint`. It is mathematically equivalent to the
multi-slice adaptive estimator for the present leading-Low density, while
avoiding repeated eta-phi grid construction at each log-`k_T` midpoint. The
multi-slice adaptive estimator remains available as an explicit reference and
as a natural fallback if future non-leading effects introduce nontrivial
`k_T` dependence.

The first tuning scan uses the deeper adaptive L9 result as a deterministic
reference for the current controlled fixtures. This is not an analytic truth;
it is a practical convergence reference for choosing prototype settings. With
the current `32 x 32`, refine-fraction-`0.1`, `2 x 2` subcell family,
adaptive L7 is stable at the few-times-`10^{-4}` level relative to L9 for
`p = 1`, `10`, and `100 GeV`, and for the text fixture. The `p = 1000 GeV`
case needs deeper refinement: L7 is about `6%` below L9, while L8 is about
`1.3%` below L9.

The provisional diagnostic gate is therefore split into two categories:

* ordinary controlled cases: adaptive L7 or better, deterministic difference
  below `1%` relative to the L9 tuning reference, and weighted-reference
  `ESS/N >= 0.05`;
* extreme high-boost stress cases: require L8/L9-level convergence checks, and
  do not rely on the current weighted reference alone when its `ESS/N` is low.

`oversamplingFactor` remains a diagnostic weighted-histogram tool. In a
diagnostic oversampled mode with factor `F`, LowMC may draw

$$
N_\mathrm{sample} \sim
\mathrm{Poisson}(F\,\lambda_\mathrm{event})
$$

and assign each generated photon an analysis weight

$$
w_\mathrm{sample} = \frac{1}{F}.
$$

Then the expected weighted count is still physical,
`E[sum w_sample] = lambda_event`, while histograms receive more entries. These
oversampled photons are weighted diagnostic samples. They must not be written
as ordinary unweighted event-record photons, and `F` must not change the
physical Poisson mean used for event-record output.

The first count-plus-shape diagnostic keeps these two pieces explicit. It
computes `lambda_event` with the adaptive-cell yield estimator, draws
`N_sample` from the Poisson rule, then runs the corrected shape-unweighting API
to obtain at least that many unit-shape candidates. The diagnostic retains the
first requested candidates and assigns the diagnostic analysis weight
`1 / oversamplingFactor`. If the shape envelope fails, or if too few accepted
shape candidates are available, the diagnostic reports an incomplete event
rather than silently writing photons.

This remains a diagnostic composition of the normalization and shape samplers.
It does not append photons to HepMC3, does not apply recoil, and does not yet
define the production event-record failure policy.

The first validation matrix for this diagnostic uses `oversamplingFactor = 100`
and three seed pairs for each controlled fixture. The provisional pass
condition is: every row has `status=ok`, `envelope_ok=1`,
`enough_candidates=1`, and `max_trial_bound_ratio < 0.8`. The stress case
`p = 1000 GeV` is still labelled separately because it relies on the L9
adaptive yield setting, but it is included in the same diagnostic matrix.

Only after both pieces are settled should the HepMC3 writer append stable
PDG-`22` photons to the event record. The first writer should preserve the
input generator event, should not apply recoil, and should mark the added
photons or their afterburner vertex with LowMC provenance and `recoil = none`.

A first smooth local-kernel diagnostic proposal is also available for
comparison with local boxes and grid-adaptive maps. It uses the same charged
particle directions and `m / E` width scale as the local-box proposal, but
replaces each flat patch by a separable compact triangular kernel in `eta` and
`phi`, truncated to the configured acceptance and normalized analytically. A
uniform floor again keeps support nonzero over the full acceptance. This
proposal is diagnostic-only and weighted; it is not a majorant.

## Eta-Phi Sampling

At fixed `k_T`, the main sampling problem is the two-dimensional distribution in
`eta` and `phi`.

For a relativistic charged particle,

$$
p_i \cdot k
\sim
E_i k_T
\left[
\cosh(\eta-\eta_i)-\cos(\phi-\phi_i)
\right].
$$

Near the charged-particle direction this becomes approximately

$$
p_i \cdot k
\sim
E_i k_T
\frac{(\Delta\eta)^2+(\Delta\phi)^2+\theta_{0,i}^2}{2},
$$

with a mass scale

$$
\theta_{0,i} \sim \frac{m_i}{E_i}.
$$

The structures in `eta` and `phi` therefore become narrower for highly boosted
charged particles. A fixed grid with a fixed rejection-envelope safety factor can
work for typical momenta but is not robust across a wide range of particle
velocities.

The preferred direction for LowMC is:

- build an event-specific `eta-phi` proposal map,
- use the selected charged-particle directions and mass-over-energy scales to
  guide refinement,
- sample cells according to estimated cell integrals,
- evaluate the exact coherent current at every proposed photon point,
- return weighted photons first, with optional unweighting only after the
  proposal and envelope have been validated.

The full coherent factor can contain interference suppressions and enhancements,
so charged-particle-centered proposals are useful for efficiency but cannot
replace exact evaluation of `S(k)`.

Sampler diagnostics should include weight spread, acceptance rate, refined-cell
counts, and checks that any local envelope is conservative.

## Fixed-kT Eta-Phi Proposal Interface

The fixed-`k_T` `eta-phi` proposal is responsible for the conditional part of
the factorized proposal density,

$$
q_{\eta\phi}(\eta,\phi \mid k_T,\mathrm{event}),
$$

with respect to the measure

$$
d\eta\,d\phi.
$$

For a given event, `k_T`, and `PhotonPhaseSpace`, it should sample only inside
the configured `eta` and `phi` acceptance and return:

- sampled `eta`,
- sampled `phi`,
- proposal density `q_eta_phi(eta, phi | k_T, event)`.

The proposal density must be normalized over the accepted `eta-phi` region at
fixed `k_T`:

$$
\int_{\eta_\min}^{\eta_\max}
d\eta
\int_{\phi_\min}^{\phi_\max}
d\phi\,
q_{\eta\phi}(\eta,\phi \mid k_T,\mathrm{event})
=
1.
$$

The full proposal density used in `Photon::weight` is then

$$
q_\mu(k_T,\eta,\phi)
=
q_{k_T,\mu}(k_T)
\,
q_{\eta\phi}(\eta,\phi \mid k_T,\mathrm{event}).
$$

The fixed-`k_T` proposal should know the selected charged particles, because the
useful proposal shape depends on their directions and mass-over-energy scales.
It must still return a density, not a weight and not a majorant. Exact
evaluation of `S(k)` remains the responsibility of the weighted photon layer.

The `phi` coordinate is periodic. Implementations should handle wrap-around at
the chosen `phi` acceptance boundary consistently, for example when a charged
particle direction lies close to `phiMin` or `phiMax`.

The minimal baseline proposal may be uniform in `eta-phi`:

$$
q_{\eta\phi}
=
\frac{1}{(\eta_\max-\eta_\min)(\phi_\max-\phi_\min)}.
$$

This is useful as a correctness and diagnostics baseline, but it is not expected
to be efficient for highly boosted particles.

Expected diagnostics for a fixed-`k_T` proposal include:

- support coverage inside the configured acceptance,
- sampled-weight spread after combining with the exact density,
- effective sample size or an equivalent weight-quality measure,
- optional refined-cell counts and cell-integral summaries for adaptive
  implementations,
- optional envelope checks only for later unweighting studies.

## Diagnostic Particle Text Input

The fixed-grid and shape-unweighting diagnostics can read a small deterministic
charged-particle text format. It is intended for hand-written toy fixtures and
sampler studies, not as a replacement for event-record adapters.

Each non-comment line contains one charged particle:

```text
charge inOutSign mass p eta phi
```

where `charge` is in positron-charge units, `inOutSign` is `incoming` or
`outgoing`, `mass` and momentum magnitude `p` are in GeV, and `eta, phi` define
the momentum direction. Blank lines and `#` comments are ignored. The
diagnostic converts this representation into the core `(E, px, py, pz)`,
charge, and `inOutSign` input.

The format is deliberately explicit about mass and sign convention. It does not
define HepMC3 output, particle-selection policy, or physical photon
multiplicity. Its purpose is to let future validation topologies be checked in
as readable fixtures instead of adding a new hardcoded C++ topology for every
sampler study.

## Numerical Requirements

- use physical particle masses through the supplied four-momenta
- do not add an artificial angular regulator
- preserve event-wide interference by summing all selected charged particles in
  one current
- treat near-zero `p_i . k` denominators carefully and document any numerical
  protection that is introduced

## Expected Tests

- charge-neutral two-particle toy event
- finite-mass collinear regularization
- `O(N)` current evaluation rather than explicit `O(N^2)` pair summation
- sharply boosted charged-particle direction that tests sampler stability

## Open Questions

- Adaptive `eta-phi` proposal strategy for efficient weighted sampling
- Exact data structures for adaptive cell maps and refinement criteria
- Which proposal diagnostics should be release-blocking
- Numerical threshold policy for very small `p_i . k`
- Whether charge conservation should be checked, warned, or left to the caller
