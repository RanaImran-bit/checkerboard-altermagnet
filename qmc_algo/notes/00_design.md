# Design note 00 — Evolutionary / Neural trial-state AFQMC

Status: DESIGN (no code yet). Author seed: owner's preliminary ideas —
(1) "walker as parents" (evolutionary algorithm), (2) "neural-network quantum
state as walker". This note turns those into a concrete, ED-gateable method
aimed at the repo's open blocker: the **d-wave constrained-path (CP) bias**.

---

## 1. The unifying observation

Constrained-path AFQMC population control is *already* an evolutionary algorithm:

| Evolutionary algorithm | AFQMC population control |
|---|---|
| individual / genotype     | walker = Slater determinant \|φ_w⟩ |
| mutation                  | auxiliary-field propagation e^{-Δτ ĥ(x)} |
| fitness / selection       | weight w_w from importance sampling + branching |
| reproduction (asexual)    | high-weight walker spawns copies; low-weight dies |
| **crossover (missing)**   | — walkers never interbreed |
| **rich genotype (missing)** | — genotype is a *single* determinant |

So the two preliminary ideas are exactly the two missing generalizations:

- **Idea 1 — "walkers as parents":** add *crossover*. Recombine two parent
  walkers into offspring. Combining two determinants is *not* a determinant — it
  is a 2-determinant state. So crossover **natively produces multideterminant
  structure**, which is precisely what kills the CP bias (this repo,
  `validate_casci_ed.py` / `validate_casscf.py`: full CAS → ED exactly; reduced
  CAS(2,2) cut the bias ~85× with NO ED).
- **Idea 2 — "NQS as walker":** enrich the genotype from a single Slater
  determinant to a **neural-network quantum state** (or neural-backflow
  determinant), which can carry pairing/backflow correlation a single
  determinant cannot.

The open problem the repo named — *"generate good multideterminant trials
without ED"* — is an **optimization/search** problem. Evolution and neural
generative models are two ways to solve it. That is the thesis.

---

## 2. Why this targets the d-wave bias specifically

The d-wave failure (HANDOFF: half-filling + anisotropy, equal-time vertex
ED +23 vs CPMC −2, sign-flipped) is a **trial-node** problem: the fixed
free-fermion single determinant — and even the adaptive natural-orbital trial —
have the wrong nodal structure for the pairing channel. The proven fix is
multideterminant, but the natural-orbital (normal-state) construction does not
reach the *pairing* manifold.

The pairing nodal structure lives in determinants connected by **pair
excitations** (k↑,−k↓ → k'↑,−k'↓). An evolutionary search whose gene pool is
**seeded with a projected-BCS / d-wave reference** and whose crossover *combines
pair-excited determinants* builds exactly that manifold — and the variational
fitness will retain it iff pairing is genuinely favored in the anisotropic
regime. That is the mechanism by which this can succeed where adaptive failed.

---

## 3. The unbiasedness firewall (critical design choice)

Selection toward low energy is *non-variational* and biased — this is the same
caveat that dogged the adaptive trial. We neutralize it with a **two-population
firewall**:

- **P1 — sampling walkers.** The ordinary AFQMC ensemble. Estimators
  (mixed / back-propagated) are computed from P1 exactly as today. **Untouched,
  unbiased given the constraint.**
- **P2 — trial archive.** The evolved/neural population. Its ONLY job is to
  define the trial ψ_T used for (a) the constraint and (b) the importance
  function. It is freely optimized toward low energy.

Key fact that makes this safe: a *better* trial can only ever **reduce** the CP
bias, and the bias → 0 as ψ_T → exact (monotone, demonstrated in the CASCI
sweep). P2's optimization bias never leaks into the estimator — it only tightens
the constraint. This is the same logic that justifies any improved trial; we are
just generating a better one automatically.

---

## 4. Algorithm A — Evolutionary Multideterminant Trial (EMT-AFQMC)

The gate-able, first-to-build method. Reuses `casci.py` (Slater–Condon CI,
integral transform, basis-invariant E0 = ED to 1e-14) and `set_multidet`.

**Genotype (two options, build option (a) first):**
- (a) *Occupation string* in a fixed orbital basis (e.g. natural orbitals) →
  **selected CI**. Discrete, classic-GA friendly. Crossover = exchange occupied
  orbitals between two strings; mutation = single particle–hole excitation.
  Fitness via Slater–Condon (reuse `casci.py` directly).
- (b) *Nonorthogonal determinant* (continuous orbital matrix) → nonorthogonal
  CI. Matches AFQMC walkers (which ARE nonorthogonal determinants), so current
  P1 walkers can be injected as genes. Crossover = column/orbital swap or
  orbital averaging; mutation = small orbital rotation OR one AFQMC propagation
  step (a physically meaningful, energy-lowering mutation). Fitness via
  generalized (Löwdin) eigenproblem.

**Fitness = marginal energy lowering** ΔE when a candidate is added to the
archive (the CIPSI / Monte-Carlo-CI selection criterion). Natural, cheap
(perturbative estimate), no ED.

**Loop:**
1. **Seed** archive with: free-fermion GS determinant + symmetry-broken
   references (AFM, and a **projected d-wave BCS** determinant) to plant pairing
   genes.
2. **Reproduce:** generate offspring by crossover + mutation from archive
   members *and from current P1 AFQMC walkers* ("walkers as parents").
3. **Score** offspring by ΔE.
4. **Select:** add top offspring; cull to archive size K enforcing **diversity**
   (drop determinants with overlap > threshold — niching, prevents genetic
   collapse).
5. **Install** ψ_T = Σ_m c_m |D_m⟩ (c from diagonalizing the K×K CI) as the
   AFQMC constraint for the next window; iterate.

**ED gates (reuse existing harness):**
- K → full archive must reproduce ED exactly (hard gate, via `casci.py`).
- Reduced K must drive CP bias monotonically down (replicate the 85× CAS(2,2)
  result) with NO ED.
- **NEW claim to prove:** evolutionary selection reaches that bias reduction
  with fewer determinants, no ED, AND captures the **d-wave nodes at
  half-filling** where the natural-orbital adaptive trial failed. This is the
  publishable result.

---

## 5. Algorithm B — Neural trial / neural walker (NQS)

**B1 — Self-learning generative trial (recommended NN entry point).**
A pure NQS in the configuration basis cannot be cheaply overlapped with an AFQMC
determinant walker (⟨ψ_T|φ_w⟩ = Σ_R ψ_T*(R)·det-amp(R), exponential). The clean
coupling: use the NQS as a **generative model for the multideterminant
expansion** — sample/enumerate its large-amplitude determinants and feed them as
the gene pool for Algorithm A. The network *learns the distribution of important
determinants*; evolution refines and the CI makes it exact-in-the-limit. This
**folds B into A** and keeps every overlap tractable. (Conceptually: a
self-learning Monte Carlo proposal for the trial.)

**B2 — Neural-backflow Slater trial.**
ψ_T = det[Φ(R) + NN-backflow(R)] or det[Φ]·exp(J_NN(R)). More expressive single
object; train by VMC, then use as constraint. Overlap with AFQMC walkers is
tractable only for specific determinant-expandable forms — flagged as a research
sub-track, not the first build.

**B3 — NQS as the walker (the literal "idea 2", research-grade).**
Each walker is an NQS; projector imaginary-time evolution in NQS space
(t-VMC / stochastic-reconfiguration style), auxiliary-field one-body operators
applied then re-fit. Heavy, with its own re-fitting bias. Keep as a later
exploratory track once A is validated.

---

## 6. Staging

1. **A(a) — discrete-GA selected-CI trial.** Smallest step on top of `casci.py`.
   Gate: full = ED; reduced beats adaptive on a small cluster.
2. **A(b) — nonorthogonal-determinant GA with walker injection.** Wires in
   "walkers as parents" literally; targets the d-wave half-filling case.
3. **B1 — neural generative gene pool.** Replace/augment the GA proposal with a
   learned generative model (start with something light: an autoregressive or
   RBM amplitude over occupation strings).
4. **B2 / B3 — neural-backflow trial / neural walker.** Exploratory.

Each stage is independently ED-gated and independently useful.

---

## 7. Pitfalls (name them up front)

- **Genetic collapse / loss of diversity** → enforce niching (overlap-threshold
  culling, fitness sharing). A QMC trial needs spread, not a single genotype.
- **Non-variational selection bias** → handled by the §3 two-population firewall.
- **Crossover of determinants ≠ determinant** → genotype must be multidet (CI)
  or crossover restricted to orbital swaps; spell this out per representation.
- **NQS↔determinant overlap cost** → why B1 (generative → CI) is the tractable
  coupling, not a pure configuration-basis NQS trial.
- **None of this removes the sign problem.** It reduces the *constrained-path
  bias* by improving the trial. State that honestly in any writeup.
- **Cost.** Fitness eval (CI diagonalization / nonorthogonal overlaps) must stay
  cheap (perturbative ΔE, capped archive K) or the outer loop dominates.

---

## 8. Literature anchors (to position the work)

- Carleo & Troyer, *Science* 2017 — neural-network quantum states (RBM).
- Luo & Clark, *PRL* 2019 — neural backflow for fermions; Pfau et al.
  (FermiNet) 2020; Hermann et al. (PauliNet) 2020 — neural VMC.
- Liu, Qi, Meng, Fu, *PRB* 2017 — self-learning Monte Carlo (learned proposals).
- Selected-CI / CIPSI / Monte-Carlo-CI — determinant-selection by ΔE fitness
  (the classical analogue of our GA fitness).
- AFQMC trial-state improvement / phaseless-bias reduction via better/
  multideterminant trials (this repo's CASCI result is the in-house anchor).

---

## 9. Recommended first move

Build **A(a)** — a discrete genetic-algorithm selected-CI trial on top of
`casci.py` — because it (i) reuses a validated exact gate, (ii) directly tests
the core thesis (evolution generates good multidets without ED), and (iii) is
small. Then escalate to A(b) to attack the d-wave half-filling case with walker
injection, which is the result that would actually matter for the manuscript.
