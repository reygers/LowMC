# PHYSICS

## Purpose

This document describes the physical model implemented in LowMC.

The goal of LowMC is to generate coherent soft electromagnetic radiation from
arbitrary scattering events using the leading-order Low theorem.

The emphasis is on a faithful implementation of the underlying physics rather
than on a particular event generator.

---

# Scope

Current scope:

* leading-order Low theorem
* coherent soft-photon emission
* leading-current soft virtual photons and weighted dielectron conversion
* arbitrary charged external particles
* event-by-event calculation

Not yet included:

* multiple-photon exponentiation
* recoil prescriptions
* higher-order corrections
* virtual-photon recoil or finite-$q$ corrections
* unweighted dielectron event-record output

---

# Soft-photon theorem

For a process containing charged external particles

$$
\{p_i, Q_i\},
$$

the leading soft-photon amplitude factorizes as

$$
\mathcal M_\gamma
=
e J^\mu(k)
\varepsilon_\mu
\mathcal M_0,
$$

where

$$
J^\mu(k)
=
\sum_i
\eta_i Q_i
\frac{p_i^\mu}{p_i\cdot k}.
$$

Here

* $Q_i$ is the particle charge in units of the positron charge $e$,
* $\eta_i=-1$ for incoming particles,
* $\eta_i=+1$ for outgoing particles.

All four-momenta are physical four-momenta with positive energy. The incoming
minus sign is carried by $\eta_i$; incoming momenta should not be negated.

The current satisfies

$$
k_\mu J^\mu = 0
$$

provided charge is conserved.

For a complete coherent current, this corresponds to

$$
\sum_i \eta_i Q_i = 0.
$$

If a configured particle selection violates this condition, LowMC should report
that the selected current is not charge conserving.

---

# Photon intensity

Summing over photon polarizations gives

$$
S(k)
=
-J_\mu(k)J^\mu(k).
$$

Equivalently,

$$
S(k)
=
\sum_{i,j}
\eta_iQ_i
\eta_jQ_j
\frac{-p_i\cdot p_j}
{(p_i\cdot k)(p_j\cdot k)}.
$$

This is the quantity implemented in LowMC.

For computational efficiency the current

$$
J^\mu
$$

should always be computed first before forming

$$
S=-J^2.
$$

This reduces the computational cost of evaluating one photon direction from
$O(N^2)$ to $O(N)$.

---

# Event-generator ensemble

The inclusive Low formula used in experimental comparisons contains an integral
over the hadronic final-state distribution,

$$
\int d^3\vec p_1 \ldots d^3\vec p_{N_f}
\frac{dN_\mathrm{hadrons}}
{d^3\vec p_1 \ldots d^3\vec p_{N_f}}.
$$

In LowMC this distribution is supplied by the input event generator or event
sample.

LowMC does not model the underlying hadronic distribution itself. Instead, each
input event is treated as one Monte Carlo draw from that distribution. For that
event, LowMC evaluates the coherent factor

$$
S_\mathrm{event}(k) = -J_\mathrm{event}^2(k).
$$

Inclusive photon spectra are obtained by averaging the generated photon weights
or yields over the input event ensemble.

This is why the physics core can remain generator-independent: Pythia, HepMC, or
other event sources provide the hadronic event sample, while LowMC provides the
event-level soft-photon factor.

---

# Differential spectrum

The leading-order photon density is

$$
\frac{dN_\gamma}
{k_T\,dk_T\,d\eta\,d\phi}
=
\frac{\alpha}{(2\pi)^2}
S(k).
$$

The implementation samples this distribution within user-defined phase-space
limits.

At fixed photon direction the leading soft current scales as $J \propto 1/k_T$,
so $S(k) \propto 1/k_T^2$. Combined with the physical measure
$k_T\,dk_T$, this gives the familiar soft behavior

$$
dN_\gamma \propto \frac{dk_T}{k_T}.
$$

This motivates logarithmic `k_T` sampling inside the configured physical range.

For weighted generation, photon weights should be interpreted with respect to
the same physical measure. If photons are sampled from a normalized proposal
density $q_\mu(k_T,\eta,\phi)$ defined with respect to
$k_T\,dk_T\,d\eta\,d\phi$, then the event-level photon weight is

$$
w =
\frac{\alpha}{(2\pi)^2}
\frac{S_\mathrm{event}(k)}
{q_\mu(k_T,\eta,\phi)}.
$$

These weights are averaged over the input event-generator ensemble to obtain
inclusive spectra.

---

# Interference

An essential feature of LowMC is the inclusion of the full interference pattern.

Unlike decay-based photon generators, the current is constructed from all charged
particles contributing to the event.

Consequently,

$$
S
=
\sum_i S_i
+
\sum_{i\neq j}S_{ij},
$$

contains all interference terms automatically.

This event-wide coherence is one of the defining characteristics of LowMC.

---

# Finite masses

Charged particles are treated with their physical four-momenta.

No massless approximation is used.

Because

$$
p_i^2=m_i^2,
$$

the apparent collinear divergence is naturally regulated.

Near the particle direction,

$$
p_i\cdot k
\simeq
E_i\omega
\frac{\theta^2+\theta_0^2}{2},
$$

where

$$
\theta_0
=
\frac{m_i}{E_i}.
$$

Thus the familiar dead-cone suppression appears automatically.

No additional angular regulator should be introduced.

---

# Validity

The leading Low theorem is valid for photons whose energies are small compared
with the characteristic energy scales of the underlying process.

The implementation therefore assumes

* soft photons,
* unchanged hard scattering,
* unchanged hadronic final state.

The generated photons should be regarded as an infrared approximation rather
than exact matrix-element calculations.

---

# Event definition

The default LowMC current should include:

* incoming beam particles,
* final stable charged particles.

Incoming particles carry $\eta_i=-1$, while outgoing particles carry
$\eta_i=+1$.

This default represents the coherent current change from the incoming state to
the selected final charged state.

Other configurable choices may include

* particles before weak decays,
* particles before hadronic decays,
* beam remnants,
* generator-specific particle selections.

Different choices correspond to different physical approximations and may depend
on event-generator record conventions.

The implementation should therefore make the particle selection configurable.

---

# Extensions

## Multiple photons

Soft-photon exponentiation.

## Virtual photons and dielectrons

The first weighted dielectron implementation uses the same coherent-current
framework with a timelike virtual photon. It deliberately remains a leading
soft, no-recoil approximation.

Replace

$$
k^2=0
$$

by

$$
q^2=M_{ee}^2.
$$

For each event, the same Low current is evaluated at the virtual-photon
momentum,

$$
J^\mu(q)
=
\sum_i
\eta_i Q_i
\frac{p_i^\mu}{p_i\cdot q}.
$$

The corresponding coherent factor is

$$
S_\mathrm{event}(q) = -J_\mathrm{event}^2(q).
$$

The virtual photon converts through

$$
\gamma^\ast
\rightarrow
e^+e^-.
$$

In the pointlike approximation, the conversion from a soft virtual photon to a
dielectron pair is controlled by the Kroll-Wada factor,

$$
dN_{ee}
\simeq
dN_{\gamma^\ast}
\frac{\alpha}{3\pi}
\frac{dM_{ee}^2}{M_{ee}^2}
\sqrt{1-\frac{4m_e^2}{M_{ee}^2}}
\left(1+\frac{2m_e^2}{M_{ee}^2}\right).
$$

Balek, Pisutova, and Pisut discuss the same soft-dilepton problem as soft
virtual-photon bremsstrahlung from charged hadrons followed by conversion to an
electron-positron pair. Their mass-dependent conversion factor is equivalent to
the Kroll-Wada factor written in terms of $dM_{ee}$ rather than
$dM_{ee}^2$. Their treatment also separates transverse and longitudinal
virtual-photon contributions; the longitudinal term vanishes as
$M_{ee} \rightarrow 0$ but should be revisited before a fully differential
dielectron implementation is designed.

`generateWeightedDielectrons` samples this factorized internally converting
virtual-photon picture. Its density and proposal are defined with respect to

$$
dM_{ee}^2\,q_Tdq_T\,dy\,d\phi.
$$

Mass and transverse momentum are proposed logarithmically, while pair rapidity
and azimuth are proposed uniformly. The returned Monte Carlo weight includes
the exact proposal density in that measure.

The two-body decay is proposed isotropically in the virtual-photon rest frame.
For a conserved leading current, its angular weight is

$$
w_\Omega =
\frac{1-\beta_e^2(\hat{J}\cdot\hat{n}_e)^2}
     {1-\beta_e^2/3},
\qquad
\beta_e^2=1-\frac{4m_e^2}{M_{ee}^2},
$$

where $\vec J$ is evaluated in the pair rest frame. This is the normalized
$J_\mu L^{\mu\nu}J_\nu$ angular dependence: it averages to one, so integrating
the lepton angles reproduces the Kroll--Wada result. A charge-conserving input
selection is required for the Ward identity and this polarization
interpretation; the API reports the signed charge sum for every generation
call.

As for real photons, the hadronic event ensemble is supplied by the input event
generator or event sample. LowMC evaluates the event-level virtual-photon
factor, and inclusive dielectron spectra would be obtained by averaging over
events.

This extends LowMC's weighted core to soft electromagnetic pair production.

Limitations and open questions for this extension include:

* validity range in pair mass and pair transverse momentum,
* possible source-dependent form factors or corrections.

---

# Open physics questions

* Which particles define the coherent current?
* Best treatment of beam remnants?
* How should recoil be handled?
* Can weighted generation be used throughout?
* How should exponentiation be implemented?
* Which recoil prescription and finite-$q$ current should define a later
  beyond-leading-soft dielectron mode?
* How should LowMC interface to existing QED radiation (e.g. shower photons or
  PHOTOS) without double counting?

This document should evolve together with the physics understanding of the
project.
