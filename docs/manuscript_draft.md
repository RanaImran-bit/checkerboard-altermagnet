# Manuscript draft: Introduction, Model and Method, Results

Status 2026-08-04. Prose is LaTeX-ready. Numbers come from `data/chi_master_all.csv`
(L = 6, 8, 10, 12 x U = 0, 2, 4, 6, 8 x 6 fillings x 5 anisotropies x 6 seeds).
**All quoted results are L = 12.** L = 6 is used only where a finite-size failure is
the point. Sections marked PENDING wait on runs in progress.

---

## I. INTRODUCTION

```latex
Altermagnets are collinear magnets with compensated sublattice moments whose bands are
spin split at generic momenta, with the splitting following a form factor of definite
symmetry fixed by the crystal \cite{Smejkal2022a,Smejkal2022b}. The two magnetic
sublattices are related by a rotation rather than by translation or inversion, which
allows a momentum-dependent spin splitting without a net moment. Candidate materials
include RuO$_2$, MnTe and KRu$_4$O$_8$, and the resulting spin-split Fermi surfaces have
been proposed as a route to spin-current generation without stray fields
\cite{Smejkal2022b,Bai2024}.

Most microscopic treatments of altermagnetism build the spin splitting into the
Hamiltonian from the start, either through spin-dependent hopping or through a
mean-field decoupling that fixes the ordered moment. Both approaches answer what an
altermagnet does once it exists. Neither answers whether electronic correlations alone
produce one, or what pairing follows when they do. That distinction matters for
superconductivity, because the pairing symmetry selected by an altermagnet depends on
the symmetry of the magnetic form factor, and a form factor imposed by hand carries no
information about the interaction that would also drive pairing.

We study a model in which the altermagnetism is not imposed. The Hamiltonian is a
single-band Hubbard model on the anisotropic checkerboard (planar pyrochlore) lattice
with strictly spin-independent hopping. The lattice has two sublattices related by a
$90^\circ$ rotation, distinguished only by which diagonal of a crossed plaquette carries
the larger hopping. The anisotropy $\delta$ that distinguishes them is a geometric
parameter of the one-body problem and produces no spin splitting on its own. Spin
splitting appears only when the interaction $U$ generates magnetic order on a lattice
where $\delta \neq 0$, so both ingredients are necessary and the altermagnetism is
emergent.

The geometry also fixes the symmetry channel. On the checkerboard lattice the two
sublattices are exchanged by the same rotation that changes the sign of
$\sin k_x \sin k_y$, so the altermagnetic form factor has $d_{xy}$ symmetry. This
differs from the $d_{x^2-y^2}$ form factor of square-lattice altermagnet models, and it
raises a question that the imposed-splitting studies cannot pose: if the magnetic order
carries $d_{xy}$ symmetry, which pairing channel survives?

We answer this with constrained-path quantum Monte Carlo (CPQMC), benchmarked against
exact diagonalization, on lattices up to $12 \times 12$. We measure the connected vertex
pairing susceptibility in four channels, on-site $s$, extended $s$, $d_{x^2-y^2}$ and
$d_{xy}$, across filling, anisotropy and interaction strength. Three results follow.

First, the pairing vertex vanishes identically at $U = 0$, to a numerical precision of
$10^{-13}$, at every lattice size and every $\delta$. Whatever the anisotropy does to
the band structure, it produces no pairing vertex without interactions.

Second, the anisotropy suppresses both $d$-wave channels, but it suppresses $d_{xy}$
several times more strongly than $d_{x^2-y^2}$. At $n \approx 0.78$ and $U = 4$ on the
$12 \times 12$ lattice, $\chi_{d_{xy}}$ falls from $-1.68$ to $-5.99$ as $\delta$ goes
from $0$ to $0.4$, while $\chi_{d_{x^2-y^2}}$ stays positive throughout. The channel
that shares the symmetry of the magnetic order is the channel that is destroyed, so the
selection is symmetry resolved rather than a uniform weakening of pairing.

Third, this selection is stable against finite-size effects. It holds at every
combination of the four lattice sizes and four nonzero interaction strengths we studied,
and the margin grows with $L$.

We also report a finite-size failure that bears on how such results should be read. On a
$6 \times 6$ lattice, $\chi_{d_{x^2-y^2}}$ appears to be enhanced by the anisotropy, at
a statistical significance of $5.3\sigma$. Every larger lattice reverses the sign of that
trend. We give the numbers in Sec.~\ref{sec:finite_size} because the smallest lattice
sizes remain common in this literature.
```

**Notes for revision.** The `\cite` keys are placeholders and need the real bibliography.
The existing abstract in `title_abstract_captions.md` states that the splitting "scales as
the product $M\delta$". That claim needs revising before submission. The measurement at
$U = 0$ gives $M = 0.24$ to $0.43$ where there is no order at all, so $M$ as defined is
the structure-factor amplitude and not an order parameter. See Sec.~V.

---

## II. MODEL AND METHOD

### A. Model

```latex
We study the single-band Hubbard model
\begin{equation}
  H = -\sum_{\langle ij \rangle \sigma} t_{ij}
        \left( c^{\dagger}_{i\sigma} c^{\phantom{\dagger}}_{j\sigma} + \mathrm{H.c.} \right)
      + U \sum_{i} n_{i\uparrow} n_{i\downarrow} ,
  \label{eq:H}
\end{equation}
on the checkerboard lattice. Here $c^{\dagger}_{i\sigma}$ creates an electron of spin
$\sigma$ on site $i$, $n_{i\sigma} = c^{\dagger}_{i\sigma} c^{\phantom{\dagger}}_{i\sigma}$,
and $U > 0$ is the on-site repulsion. All energies are in units of the nearest-neighbor
hopping $t = 1$.

The hopping $t_{ij}$ equals $t$ on every nearest-neighbor bond. Diagonal bonds exist only
on alternating plaquettes, and their amplitude is
\begin{equation}
  t_{\pm} = t' \pm \delta ,
  \label{eq:tpm}
\end{equation}
where $t_{+}$ acts along one diagonal and $t_{-}$ along the other. Which diagonal carries
$t_{+}$ alternates between the two sublattices, defined by the parity of $x + y$. We take
$t' = -0.3$ throughout.

Two properties of Eq.~\eqref{eq:tpm} matter for what follows. First, $t_{ij}$ is
independent of spin, so the one-body problem has no spin splitting for any $\delta$.
Second, at $\delta \neq 0$ the two sublattices are inequivalent and are exchanged by a
$90^\circ$ rotation, not by a translation. The Hamiltonian is then invariant under the
diagonal translation $(1,1)$ but not under the single-site translation $(1,0)$, so the
unit cell contains two sites and the spectrum has two bands. We have verified this
directly from the hopping matrix. At $\delta = 0$ the two diagonals become equal, the
sublattices merge, and the single-site translation is restored.

Because the sublattice exchange is the same $90^\circ$ rotation under which
$\sin k_x \sin k_y$ changes sign, any sublattice-staggered spin polarization carries
$d_{xy}$ symmetry. The altermagnetic form factor is therefore fixed by the lattice and is
not a free choice.
```

### B. Constrained-path quantum Monte Carlo

```latex
We solve Eq.~\eqref{eq:H} with constrained-path quantum Monte Carlo (CPQMC)
\cite{Zhang1997}, a ground-state auxiliary-field method in which the fermion sign problem
is controlled by constraining the random walk to the half-space of Slater determinants
with positive overlap against a trial wave function. The constraint removes the
exponential growth of the variance and introduces a bias that depends on the trial state.
Results reported here are therefore obtained within the constrained-path approximation
rather than being numerically exact, and we benchmark against exact diagonalization on
lattices where that is feasible.

The trial wave function is the free-electron ground state of the hopping matrix in
Eq.~\eqref{eq:H}. We use a Trotter step $\Delta\tau = 0.05$, $160$ walkers, $60$
equilibration steps and $40$ measurement blocks, with population control and
reorthogonalization every ten steps. Each reported quantity is averaged over six
independent random seeds, and error bars are the standard error over those seeds.
Lattices are $L \times L$ with periodic boundary conditions and $L = 6$, $8$, $10$ and
$12$.

Pairing is measured through the imaginary-time-displaced pair correlation function
\begin{equation}
  C_{a}(\tau) = \langle \Delta_{a}(\tau) \Delta^{\dagger}_{a}(0) \rangle ,
  \label{eq:Ctau}
\end{equation}
with the singlet pair operator
\begin{equation}
  \Delta^{\dagger}_{a}
    = \sum_{ij} f_{a}(i,j)\, c^{\dagger}_{i\uparrow} c^{\dagger}_{j\downarrow} .
  \label{eq:Delta}
\end{equation}
The pair carries zero center-of-mass momentum, so all susceptibilities reported here are
evaluated at the $\Gamma$ point, and the internal symmetry is carried entirely by the form
factor $f_{a}$. We use four form factors: on-site $s$, with $f$ diagonal; extended $s$,
with $f = 1$ on all nearest-neighbor bonds; $d_{x^2-y^2}$, with $f = +1$ on $x$ bonds and
$-1$ on $y$ bonds; and $d_{xy}$, with the corresponding sign structure on the two
diagonals.

The susceptibility is the integral
$\chi_{a} = \int_{0}^{\tau_{\mathrm{max}}} C_{a}(\tau)\, d\tau$, evaluated by the
trapezoidal rule over $\tau_{\mathrm{max}} = 16 \Delta\tau$. We report the connected
vertex contribution
\begin{equation}
  \chi^{\mathrm{vertex}}_{a} = \chi_{a} - \chi^{\mathrm{bubble}}_{a} ,
  \label{eq:vertex}
\end{equation}
where the bubble is built from the ensemble-averaged time-displaced Green functions. The
vertex isolates the interaction-induced part of the pair correlation, and a positive value
indicates an attractive contribution in that channel. Earlier benchmarking against exact
diagonalization showed that the equal-time vertex is sensitive to the constrained-path
bias while the integrated vertex of Eq.~\eqref{eq:vertex} is not, so we use the integrated
form throughout.

Equation~\eqref{eq:vertex} provides a stringent internal check. At $U = 0$ the full
correlation function is its own bubble, so $\chi^{\mathrm{vertex}}_{a}$ must vanish
identically. We recover this to between $10^{-13}$ and $10^{-14}$ in all four channels at
every lattice size, filling and anisotropy.
```

---

## III. RESULTS: CHANNEL SELECTION

```latex
Figure~\ref{fig:maps} shows the four vertex susceptibilities over the filling-anisotropy
plane at $L = 12$ for $U = 0$ to $8$. The $U = 0$ row is zero within numerical precision
and serves as a validation of Eq.~\eqref{eq:vertex} rather than as data.

The on-site $s$ channel is repulsive at every point with $U > 0$, as expected for an
on-site interaction. The extended $s$ channel is attractive over most of the plane. The
two $d$-wave channels behave differently from each other, and that difference is the
central result.

At the doped filling $n \approx 0.78$ and $U = 4$, the anisotropy suppresses both $d$-wave
channels. Table~\ref{tab:channels} gives the values at $L = 12$. Over the range
$\delta = 0$ to $0.4$, $\chi_{d_{x^2-y^2}}$ falls from $+3.79$ to $+0.87$ but remains
positive, while $\chi_{d_{xy}}$ starts negative at $-1.68$ and falls to $-5.99$. The two
channels are not in competition in the sense that one grows while the other shrinks.
Both weaken, and $d_{x^2-y^2}$ prevails because it survives.

This is what a $d_{xy}$ magnetic form factor should do. The altermagnetic order parameter
transforms as $\sin k_x \sin k_y$, so it couples directly to the pairing channel of the
same symmetry and leaves the orthogonal channel comparatively untouched. The selection is
therefore symmetry resolved.

Comparing channels directly, $d_{x^2-y^2}$ exceeds $d_{xy}$ at $n \approx 0.78$ in all
sixteen combinations of the four lattice sizes and four nonzero interaction strengths,
and the margin grows with $L$.

The interaction dependence is not monotonic. Both $d$-wave channels are strongest at
intermediate coupling and weaken by $U = 8$. Locating the maximum for each channel gives
$U = 4$ for extended $s$ and $U = 2$ for $d_{x^2-y^2}$, and these values are identical at
all four lattice sizes. Weak coupling produces no vertex, and strong coupling localizes
the moments, so pairing is largest in between.
```

### PENDING for Sec. III
- Half-filling behavior. At $L = 12$ the sign of $\chi_{d_{xy}} - \chi_{d_{x^2-y^2}}$ at
  $n = 1$ changes near $\delta \approx 0.2$ at every $U$. Periodic boundary conditions make
  $n = 1$ an open-shell case with a degenerate trial state, so this needs the twist-averaged
  run before it can be stated.
- Dense-filling scan (run in progress) will let us quote where in $n$ the selection changes,
  currently resolved only to $\pm 0.1$.

---

## IV. RESULTS: FINITE-SIZE BEHAVIOR
`\label{sec:finite_size}`

```latex
The suppression of $d_{xy}$ by the anisotropy strengthens monotonically with system size.
At $n \approx 0.78$ and $U = 4$, the change in $\chi_{d_{xy}}$ between $\delta = 0$ and
$\delta = 0.4$ is $-1.68$ at $L = 6$, $-2.45$ at $L = 8$, $-3.40$ at $L = 10$ and $-4.31$
at $L = 12$. The effect is not a small-cluster artifact.

The $d_{x^2-y^2}$ channel behaves differently, and the difference is instructive. On the
$6 \times 6$ lattice, $\chi_{d_{x^2-y^2}}$ increases with $\delta$, by
$+0.48 \pm 0.09$ between $\delta = 0$ and $\delta = 0.4$. Taken alone, this supports a
picture in which the anisotropy transfers weight from one $d$-wave channel to the other.
The larger lattices do not support it. The same quantity is $-2.75 \pm 0.37$ at $L = 8$,
$-1.85 \pm 0.42$ at $L = 10$ and $-2.93 \pm 0.66$ at $L = 12$. All three are suppressions,
each more than four standard errors from zero, and the $L = 6$ enhancement is more than
five standard errors in the opposite direction. The trend at $L = 6$ is significant and it
is wrong.

The magnitudes are also strongly size dependent even where signs agree. At $\delta = 0$,
$U = 4$ and $n \approx 0.78$, $\chi_{d_{x^2-y^2}}$ is $+0.43$ at $L = 6$ and between
$+2.61$ and $+3.79$ on the larger lattices. Quoting a magnitude from the smallest lattice
would understate the attraction by a factor of six.

A second discrepancy appears in the leading-channel assignment. At $U = 4$, ranking all
four channels point by point over the $(n,\delta)$ plane puts extended $s$ first at $13$
of $25$ points on the $6 \times 6$ lattice, whereas $L = 8$, $10$ and $12$ put
$d_{x^2-y^2}$ first at $13$, $15$ and $15$ points respectively. At $U = 8$ the four sizes
agree.

We report these comparisons because $6 \times 6$ clusters remain in use for this class of
model, and because the failure is not visible from within a single lattice size. The
$L = 6$ trend has a small error bar and a consistent sign across seeds. Only the
comparison across sizes reveals it.
```

---

## V. RESULTS: PAIRING AND THE MAGNETIC ORDER PARAMETER

**PENDING.** The $L = 12$ magnetic run is in progress. Current status of the argument:

- $M = \sqrt{S(\pi,\pi)}$ is **not** an order parameter. At $U = 0$ we measure
  $M = 0.24$ to $0.43$ across fillings, where by construction there is no order. $M$ is the
  amplitude of the structure factor.
- $S(\pi,\pi)$ is extensive once long-range order develops, so the size-comparable moment is
  $m = \sqrt{S(\pi,\pi)/N}$. The new runs record both.
- The candidate order parameter is $m_{\mathrm{AM}} = \Delta m\, \delta$ with
  $\Delta m = m - m_{U=0}$, which vanishes at $U = 0$ and at $\delta = 0$ as required.
- At $L = 6$, $\chi_{d_{xy}}$ tracks $m_{\mathrm{AM}}$ at fixed filling with Spearman
  $\rho = -0.91$, $-0.83$ and $-0.68$ at $n = 0.556$, $0.667$ and $0.778$, while
  $\chi_{d_{x^2-y^2}}$ shows no such dependence ($|\rho| \le 0.29$). **This is $L = 6$ only
  and is not reportable until reproduced at $L = 12$**, for the reasons in Sec.~IV.
- A projection $\Psi = N^{-1} \sum_{\mathbf q} \sin q_x \sin q_y\, S(\mathbf q)$ was tested as
  a general order parameter and **failed** its null test. At $L = 14$ its $\delta = 0$ value
  is the same size as its $\delta \neq 0$ values and shows no monotonic trend. The reason is
  structural: $(\pi,\pi)$ is invariant under the $90^\circ$ rotation, so any rotation-odd
  weight vanishes exactly where the magnetic weight sits. We do not use it.

---

## VI. NOT YET WRITTEN

- **Sec. VI, Spectral properties.** CPQMC cannot reach this. The usable imaginary-time
  window is $\tau \approx 1$, beyond which $G(\mathbf{k},\tau)$ becomes negative. Reaching
  $\tau \approx 3$ would require of order $10^{5}$ times more samples. This is the
  exponential signal-to-noise problem for time-displaced fermion Green functions, and the
  route is a determinant QMC calculation with analytic continuation rather than CPQMC.
- **Sec. VII, Discussion and conclusions.**
- **Appendix, ED benchmark.** Numbers exist in the earlier validation scripts and need
  collecting into a table.
```
