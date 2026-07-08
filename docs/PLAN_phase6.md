# Phase 6 Plan — Observable parity, susceptibility, adaptive trial WF

Status: PLANNED (2026-06-24). Scope = three workstreams requested by the owner.
Driver: the PRL LE20050 altermagnet appeal needs the Python port (`pyqmc/cpqmc.py`)
to measure the *same* order parameters the Fortran code does (d-wave pairing,
magnetization), plus dynamical susceptibilities, and to test whether an adaptive
(self-consistent) trial wavefunction reduces the constrained-path bias in the
sign-problematic doped regime.

Current state going in:
- `pyqmc/cpqmc.py` (461 lines, single `CPMC` class) validates energy on all three
  interaction channels (uxx/uxy/v) vs ED, and equal-time correlations (G^s_ij,
  charge <n_i n_j>, spin <S^z_i S^z_j>) via back-propagation (`run_bp_obs`).
- Trial WF is FIXED: free-electron ground state `psiT_up/dn = eigh(K)[:, :nup/ndn]`,
  built once in `__init__`; the constraint/importance sampling uses it every step.
- Fortran `code/src/meas.f90` (`correl` + `pair`) measures CDW, SDW-z, SDW-x and a
  large pairing-symmetry set (s, d, p, extended-s `so`, and BCS-like `*b` channels,
  plus orbital-resolved d1/d2/d12) via the 4-bond form factors in `Initpair`
  (sf/df/pf/ppf/ddf). ED references: `ed/altermagnet_ed.py` (`ed_energy`,
  `ed_correlations`).

Validation philosophy (unchanged): every new observable gets an exact ED target on a
small cluster (2x2/4+4, anisotropic to avoid the degeneracy trap), compared via
`tools/compare.py` z-scores, with back-propagation to remove mixed-estimator bias.

---

## Workstream 1 — Magnetization + pairing parity (pyqmc <- Fortran), OOP tidy-up

Goal: the Python port measures the same magnetization and pairing order parameters
as `meas.f90`, structured cleanly with OOP.

### 1a. Refactor (OOP) — DONE (2026-06-24)
Split `CPMC` into `LatticeModel`, `TrialWF`, `WalkerEnsemble`, `Propagator`,
`Estimators` + a thin `CPMC` orchestrator (public surface run/run_bp/run_bp_obs
and the constructor unchanged). Regression-gated bit-for-bit by
`pyqmc/regression.py` (5 fixed-seed configs: hubbard/alt mixed, alt uxx/uxy bp,
correlation checksum) — identical before/after. The `TrialWF.update()` hook is now
the seam for WS3 (adaptive trial). CLI + benchmark import paths smoke-tested.

Original design notes (kept for 1b/1c):
Split the monolithic `CPMC` into composable pieces (single file or a small
`pyqmc/` package — decide at implementation; keep the CLI import path stable):
- `Lattice` / `Model` — holds `K`, dims, interaction term list (`_build_terms`),
  HS constants. Pure problem definition.
- `Walker` ensemble state — `phi_up/dn`, weights, overlaps; `reorthogonalize`,
  `pop_control`.
- `TrialWF` — encapsulates `psiT_up/dn` and the overlap/Green's-function machinery
  (`_ov_spin`, `_overlap`, `_green`). This seam is REQUIRED for Workstream 3
  (swap fixed -> adaptive trial behind one interface).
- `Propagator` — `step`, `_step_record`, `_term_ratio` (force bias / HS sampling).
- `Estimators` — energy (mixed + BP), correlations, and the new pairing/mag
  observables. Each estimator takes a Green's function and returns a dict.
- `CPMC` becomes the orchestrator (`run`, `run_bp`, `run_bp_obs`) wiring the above.
Keep public methods/CLI behaviour byte-identical; re-run the existing 3-channel
energy + correlation validation as a regression gate after the refactor (must still
PASS at the same z-scores) BEFORE adding anything new.

### 1b/1c. STATUS — DONE (2026-06-24)
Implemented `Estimators.pairmag_block` + `CPMC.run_bp_pairmag` (single-band square
lattice): s-wave + d_{x^2-y^2} singlet pairing structure factors, spin S(pi,pi),
local moment^2, all back-propagated with per-block error bars. ED reference =
number-conserving product operator O^dagger O ("+-|+-", QuSpin) in
`pyqmc/validate_pairmag.py`. Verified EXACT at U=0 (non-degenerate nup=ndn=1, all
observables 0.0% vs ED). Interacting (U=3): spin S(pi,pi) 0.2%, d-wave + moment ~0
captured; s-wave SF +4% residual constrained-path bias (estimator exact, bias is
CP-node from the free trial -> motivates WS3). Gotcha: avoid open-shell degenerate
fillings (e.g. nup=ndn=2 at U=0) -- use closed-shell non-degenerate. See
VALIDATION.md. Altermagnet pairing geometry / Fortran-convention match = deferred.

### 1b. Magnetization (original design notes)
- Local moment m_i^z = <n_{i,up} - n_{i,dn}> (diagonal of BP G^up - G^dn).
- Staggered magnetization M_s = (1/N) sum_i e^{iQ·r_i} m_i^z, Q=(pi,pi) (single
  band) and the altermagnet's orbital-resolved staggered moment (sublattice/orbital
  sign pattern from `sublatt`).
- Spin structure factor S(q) = (1/N) sum_ij e^{iq(r_i-r_j)} <S^z_i S^z_j>
  (reuse the already-validated `szsz` from `run_bp_obs`; just add the FT + Q-peak).
- ED target: extend `ed/altermagnet_ed.py` with m_i^z and S(Q) operators (QuSpin
  number/spin ops). Validate on 2x2/4+4 (anisotropic), uxx=2, BP.

### 1c. Pairing correlations
- Port the singlet pairing operator Delta_alpha(i) = sum_{b in bonds} f_alpha(b)
  c_{i,up} c_{i+b,dn}, pair-pair correlation P_alpha(r) = <Delta^+_alpha(i+r)
  Delta_alpha(i)>, for alpha in {s, extended-s, d_{x^2-y^2}}. The d-wave channel
  (df = +1,+1,-1,-1 on the 4 bonds) is THE order parameter for the appeal.
- Wick form from the walker BP Green's functions (mirror `pair` subroutine:
  P ~ G_up(m,n) G_dn(m1,n1) + G_up(m1,n1) G_dn(m,n) with the bond form factors).
  Start with s + d singlet (the `swave`/`dwave` paths); defer the triplet/`*b`
  BCS-vertex channels unless needed.
- Report the long-range / vertex-subtracted pairing (P_d(r) minus the
  uncorrelated G·G piece) to expose enhancement — that is the paper's claim.
- ED target: pair-pair correlator for s and d on 2x2 (exact). Validate z-scores.

Deliverables: new estimators + CLI flags (`--mag`, `--pair`), ED references,
a 3-way (ED/Fortran/Python) parity table appended to VALIDATION.md, platform
hook optional (a "Pairing/Magnetization" tab) — defer UI to a follow-up.

---

## Workstream 2 — Dynamical susceptibility (unequal-time)

### STATUS — DONE (2026-06-24)
Implemented Estimators.chi_block + CPMC.run_bp_chi (single-band): imaginary-time-
displaced staggered spin correlation C(tau_l)=<O(tau_l)O(0)>, O=sum_i (-1)^{x+y}
S^z_i, l=0..bp, and static chi_s=integral C dtau. Per-walker forward propagators
B_(l) from recorded fields; displaced GFs P=B(I-g^T), H=B^{-T}g, g(tau)=B^{-T}g B^T.
ED ref = Lehmann (pyqmc/validate_chi.py). EXACT at U=0 (caught a transpose bug that
U=0 symmetry hid -- always validate dynamics on an interacting case). U=3: C(tau)
tracks ED across the window (z<~4), windowed chi_stag ED 0.0510 vs 0.0511 PASS.
C(tau_0)=N*S(pi,pi). Pairing chi_d (4-point) = follow-on. See VALIDATION.md.

Original design notes:

Goal: measure imaginary-time-displaced correlators and the corresponding
susceptibilities chi (the unequal-time observables the Fortran does not emit).

Physics: chi_O = integral_0^beta dτ <O(τ) O(0)>, with the imaginary-time
displaced Green's function G(τ) = <T c(τ) c^+(0)> built by inserting extra
propagation steps between the back-propagated bra and the ket. In projector
(ground-state) AFQMC this is a τ-displaced measurement on the central slice:
  G_s(τ; i,j) = < psi_T | c_i e^{-τ H} c^+_j e^{-(β-τ) H} | psi_T > / <...>
approximated per walker by propagating the ket forward τ extra slices while
holding the BP bra fixed, then Wick-contracting.

Targets:
- Spin susceptibility chi_s(q,τ) and uniform chi_s = sum_τ; staggered chi_s(Q).
- Pairing susceptibility chi_d(τ) for the d-wave vertex (the dynamical analogue
  of Workstream 1c — directly supports "enhancement of d-wave pairing").
- Single-particle G(k,τ) -> spectral hints (optional; defer A(k,ω) analytic
  continuation).

Implementation:
- Extend the BP machinery: after forming `ket` and the BP `bra`, generate a chain
  of τ-displaced Green's functions G(τ_l) by reusing the recorded HS fields
  (`_step_record` already stores them) to propagate the ket forward slice by slice.
- New method `run_bp_chi(taus=...)` returning {tau, chi_s(tau), chi_d(tau)} and
  the τ-summed susceptibilities, weighted-averaged over walkers.
- Numerical care: τ-displaced estimators are noisier and sign-sensitive; use the
  same anisotropic small cluster + larger NWLKRS-equivalent (`nwalkers`).

ED target: exact chi via Lehmann representation in `altermagnet_ed.py`
(eigen-decompose H, chi_O(τ) = sum_{mn} |<m|O|n>|^2 e^{-τ(E_n-E_0)} ... ),
on 2x2. This gives an exact unequal-time curve to validate against — the cleanest
possible check for a brand-new dynamical estimator.

Deliverables: `run_bp_chi`, ED Lehmann reference, a chi_s(τ)/chi_d(τ) ED-vs-pyqmc
overlay (numbers in VALIDATION.md; optional UI chart later).

---

## Workstream 3 — Adaptive (self-consistent) trial wavefunction

### STATUS — DONE (2026-06-24)
Implemented `trial="adaptive"` in `TrialWF.update` (natural-orbital self-
consistency: weighted ensemble 1-RDM -> symmetrize -> damp toward current trial
projector -> top-occupation natural orbitals as new trial), wired into CPMC with
`trial_every` / `trial_mix` + CLI flags. Harness `pyqmc/validate_adaptive.py`.
RESULT (4x4, nup=ndn=2, U=4 degenerate open shell, ED -11.530): adaptive reduces
the bias vs ED by ~17% (mixed: +0.206->+0.170) / ~22% (bp: +0.241->+0.189);
hypothesis CONFIRMED. Moderate damping (mix=0.5) beats aggressive (mix=0.25, only
~10%) -- noisy mixed RDM. Residual bias remains (single-determinant trial). Fixed
stays the validated default; regression gate bit-for-bit identical. Next refine:
back-propagated RDM trial + multi-determinant. See VALIDATION.md. Original design:



Goal: make the trial WF switchable between FIXED (current free-electron GS) and
ADAPTIVE, where the constraint boundary is redefined each step from the current
walker generation, and test whether ADAPTIVE moves the energy closer to ED in the
sign-problematic regime (away from half-filling / finite doping), where the
constrained-path bias from a poor fixed trial is largest.

Background / rationale: constrained-path AFQMC removes the sign problem by
projecting out walkers that cross the node of <psi_T|phi>=0. The systematic bias
scales with how far psi_T is from the true ground state. A self-consistently
updated trial (the AFQMC analogue of self-consistent / released-node ideas) can
reduce that bias. This is exactly where doped clusters (where the fixed
free-fermion trial is worst) should show improvement.

Design (built on the Workstream-1a `TrialWF` seam):
- `trial="fixed"` (default, current behaviour) vs `trial="adaptive"`.
- Adaptive update: each step (or every k steps), rebuild psi_T from the current
  walker ensemble. Candidate definitions (make it a strategy, test both):
  (i) weighted average of walker Slater determinants projected to an orthonormal
      n-particle frame (population's dominant single-determinant — e.g. leading
      left singular vectors of the weighted walker stack);
  (ii) density-matrix natural orbitals: build the ensemble 1-RDM
      rho_s = sum_i w_i G^s_i / sum_i w_i, diagonalize, take the nup/ndn highest-
      occupation natural orbitals as the new psi_T^s.
  Option (ii) is cleaner/cheaper and is the recommended first cut.
- Recompute overlaps/Green's functions against the new psi_T after each update;
  guard against discontinuities (the importance function changes -> reweight or
  damp the update, e.g. psi_T <- mix(old, new); start with full replace + small
  step, monitor weight stability).
- Stability: adaptive trials can destabilize population control; reuse the
  reorthogonalization + pop-control cadence and watch for walker collapse
  (the same failure mode seen historically in the v channel).

Test protocol (the actual question being asked):
- Pick clusters WITH a sign problem: away from half filling. For the single-band
  Hubbard 4x4 at, e.g., nup=ndn=5..7 (doped), U=4..8, and the altermagnet doped
  away from 16+16. Compute exact ED energy (small enough cluster) as ground truth.
- Run FIXED vs ADAPTIVE at matched statistics (same nwalkers, dt, BP). Compare
  |E_QMC - E_ED| and the z-score for each. SUCCESS = adaptive reduces the bias
  (central value closer to ED, ideally within error) in the doped/sign-problem
  case WITHOUT regressing the already-PASS half-filled cases.
- Report a table: (filling, U) x {fixed delta/z, adaptive delta/z} in VALIDATION.md.

Caveat to document honestly: a self-consistent trial can REINTRODUCE a bias of its
own (it is no longer a fixed external constraint; the result can depend on the
update scheme / become non-variational). The test is empirical — we report whether
it helps against ED, and we keep FIXED as the validated default.

Deliverables: `TrialWF` with `fixed`/`adaptive` strategies, `--trial` CLI flag,
the doped-regime fixed-vs-adaptive-vs-ED comparison table, honest write-up of when
it helps / hurts.

---

## Sequencing & acceptance gates

1. WS1a refactor -> regression gate (3-channel energy + correlations still PASS).
2. WS1b/1c magnetization + pairing -> ED parity on 2x2 (z<~2, BP).
3. WS2 susceptibility -> ED Lehmann overlay on 2x2.
4. WS3 adaptive trial -> doped-regime fixed-vs-adaptive-vs-ED table.
Each workstream: one (or few) focused commit(s), VALIDATION.md updated, HANDOFF.md
"Phase 6" section kept current. Commit per logical step; push to origin.

## Risks / open questions
- Pairing/susceptibility sign-to-noise on small clusters may need more walkers
  (now easy on the Fortran side via QMC_NWLKRS; pyqmc takes `nwalkers`).
- τ-displaced estimator correctness is the highest-risk new physics -> the ED
  Lehmann reference is the gate; build the 2-site toy check first (mirror
  pyqmc/vee_check.py methodology).
- Adaptive trial may not beat fixed at half-filling (expected: no sign problem,
  fixed trial already good) — the interesting signal is strictly in the doped runs.
- OOP refactor must not change numerics — lock it behind the regression gate.
