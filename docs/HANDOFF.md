# QMC Validation Platform — Handoff

Last updated: 2026-06-29 (three-method d-wave campaign + tau-integrated pairing eigenvalue +
Fortran-referenced stabilization; READ THE 2026-06-29 BLOCK FIRST. Older blocks below for history.)

================================================================================
# 2026-06-29 — d-WAVE PAIRING IN THE ALTERMAGNET HUBBARD MODEL (READ FIRST)
================================================================================

Goal: assess the PRL-appeal claim "the altermagnet (tam) enhances d-wave pairing", with
numerically controlled QMC, across methods/lattice/temperature/filling. Verdict + how to
continue below. Everything committed; cluster compute on profhokin@thkclusters:250-258.

--------------------------------------------------------------------------------
## (1) METHOD DEVELOPMENT STATUS
--------------------------------------------------------------------------------
- **ED (finite-T + ground state, python/QuSpin)** — DONE/validated. `code/dqmc_py/dqmc.py:
  ed_finite_T` (Tr e^{-bH}O/Z: energy, density, S_d/S_s structure factors, chi_d/chi_s Kubo).
  Also `pyqmc/ed_aniso_scan.py`, `pyqmc/validate_*.py`. Used to gate all QMC. 2x2/4x4 only
  (full Fock space). The arbiter for correctness.
- **DQMC-python (finite-T BSS)** — DONE/validated, the workhorse. `code/dqmc_py/dqmc.py`.
  Discrete-HS, ASvQRD-stabilized EQUAL-TIME Green's (brute-validated 1e-13), Sherman-Morrison
  sweep, sign tracking. Phase-1 added STABLE time-displaced `_green_tau` (B(l,0)(I+B)^-1) ->
  suscV converges at large beta (was 13095 garbage). Measures equal-time corr + tau-integrated
  susc, d-wave AND ext-s, FULL + connected-VERTEX, k-resolved (maxk/k0/r0/rgt). Plus
  `pair_eig.py`: equal-time AND tau-integrated k-space PAIRING EIGENVALUE (the SC indicator).
  Gates: U=0 exact, energy/S_d/chi_d vs ED, beta->inf. SIGN-FREE only at tam=0 half-filling.
- **CP-DQMC (finite-T constrained-path AFQMC, python)** — DONE/validated, the KEY method
  (sign=1 at half-filling all tam). `code/ftcpqmc_py/ftcpmc.py`. Force-biased per-site
  heat-bath (`run_fb_stab`), constraint via trial density matrix, Phase-2 incremental-UDV /
  reuse of DQMC `_green_tau` for stability. Measures everything DQMC does (corrV_d/s,
  suscV_d/s) PLUS the tau-integrated pairing eigenvalue (`paireig_tau` flag). Gates: U=0
  exact; vs ED at doping (density/S_d ~1%); vs DQMC across t1-scan. Expensive (L<=10 practical;
  L=12 not feasible).
- **CPQMC (T=0 constrained-path / BP-AFQMC, python)** — DONE. `pyqmc/unified_scan.py`
  (`one_pop`). Free-particle trial; canonical (fixed N). Measures corr/susc FULL+VERTEX,
  d-wave + ext-s, k-resolved; Phase-3 BP stabilization (`modgs`-style) allows longer bp.
  Added tau-integrated pairing EIGENVALUE (`paireig` flag). CAVEAT: tau-window = bp*dt
  (default 0.8 = TOO SHORT for true T=0 susc); half-filling is OPEN-SHELL (use closed shell
  near half, or tam>0 which gaps it). The WEAKEST method for the leading-channel question.
- **DQMC-fortran (BSS)** — REFERENCE ONLY (frozen, `code/record/benchmark/DQMC/BSS.f90`).
  Source of the UDV/ASvQRD stabilization we ported (Phases 0-1: `makeg`,`udvb`,`udvbt`,
  `udvBg`,`matinv` big/small split; time-displaced `FTGt`). Documented in `docs/STABILIZATION.md`.
  Not built/run this round (gfortran INTEGER-kind shims needed; the python port is gated
  against ED instead).
- **CPQMC-fortran** — REFERENCE ONLY (frozen, `code/record/benchmark/CPQMC/cp.f90` etc.;
  `modgs` reorthonormalization is the BP-stabilization reference for Phase 3). This is the
  ORIGINAL manuscript code; appeal data lives on the cluster at amax@node-253:~/run/CPQMC/
  dopingAM/ (READ-ONLY). Not rebuilt; the python unified_scan reproduces its observables.

--------------------------------------------------------------------------------
## (2) PHYSICS FOUND + DATA LOCATIONS
--------------------------------------------------------------------------------
Model: single-band ALTERMAGNET Hubbard, U=4. tam = NN spin-dependent (d_x2-y2) hopping;
t1 = NNN spin-dependent (d_xy, anisotropic t'); tp = NNN spin-independent (isotropic t').
`am_hopping(lx,ly,t0,tam,t1,tp)` in pyqmc/cpqmc.py + code/dqmc_py/dqmc.py.

THE VERDICT (robust, multi-method, sign-controlled where it matters):
1. **d-wave is the dominant pairing channel NEAR HALF-FILLING.** Projected vertex
   susceptibility suscV_d >> suscV_s (DQMC half b10 ~11x; CP-DQMC half ~6x; CPQMC
   closed-shell 1.6-2.4x; grows with L,beta). The CORRECT leading-channel test
   (TAU-INTEGRATED pairing eigenvalue) at tam=0 half-filling gives the leading eigenvector
   = D-WAVE (DQMC d-overlap 0.70 / CP-DQMC 0.80; ext-s ~0). NB the EQUAL-TIME eigenvalue
   instead leads ext-s -> wrong observable; use the tau-integrated one.
2. **tam does NOT enhance d-wave SC -- it SUPPRESSES it.** With sign=1 (CP-DQMC) the
   tau-integrated pairing eigenvalue's d-wave character collapses immediately with tam
   (d-overlap 0.58->0.15 by tam=0.2) and lambda falls to ~0 by tam=0.5. suscV_d falls with
   tam in every clean run. (The DQMC "tam=0.1-0.2 enhancement" was a SIGN ARTIFACT, sign 0.3-0.6.)
3. **What tam DOES do = enhance the EQUAL-TIME d-wave CORRELATION** (amplitude), while
   suppressing the susceptibility -> "preformed pairs without coherence" (amplitude up,
   coherence down). tam = momentum-dependent PAIR-BREAKING (spin-split FS -> pair energy
   mismatch hardens the d-wave mode). Standard pseudogap-like physics, not a contradiction.
4. **These preformed pairs are NOT useful for an SC/PDW claim:** (a) SHORT-RANGE -- the
   equal-time corrV decays to ~0 by |R|>2 (|R|>2 avg ~0.0002 vs R=0 ~0.07; NO ODLRO);
   (b) LOW/falling susceptibility (no instability); (c) chi_d(q) DOES peak at finite q under
   tam (maxk>k0) but that peak is SMALL and does NOT grow with beta (peaks ~b6 then falls) ->
   NO PDW. All three independent tests agree: no finite-q d-wave instability.
5. **Filling dependence:** d-wave leads near half; at OVER-DOPING (n=0.72 closed shell) ext-s
   leads even at long tau-window. The SC-dome doping (n~0.87, L=14/16) was launched but may
   still be running -- check `~/qmc/results/cpqmc_nearhalf_L14.csv` on node-250.
6. **t1 (anisotropic t') SUPPRESSES d-wave** (dropped earlier); **tp (isotropic t')** weak/
   non-monotonic. So the appeal's "t'-enhancement" does not hold; the tam story is the
   amplitude-vs-coherence one above.

DATA (results/dqmc_scan/): campA_dqmc_L{68,12}.csv (DQMC backbone, tam=0 half, d-vs-s vs L,b);
campB_cpdqmc_L{6,8,10,8b10}.csv (CP-DQMC full tam, sign=1); campC_cpqmc_L{6,8,10}.csv (CPQMC
half) + campC_cpqmc_cs_L{6,8,10}.csv (CPQMC closed-shell near half); dqmc_vertex2d_{6x6,8x8_b
2/4/6}.csv + cpqmc_vertex2d_{6x6,8x8}.csv (tam x t1 vertex maps); paireig_8x8_doped_b2.csv +
paireig_tau_L6_b6.csv + petau_cpdqmc_L6_b{6,8}.csv (pairing eigenvalues); cpqmc_bp{conv,scan}
_*.csv (tau-window scans). Cluster mirror: ~/qmc/results on each of nodes 250-258 (home NOT
shared; per-node python: 250-252 ~/miniconda3, 253/258 system, 256/257 /opt/anaconda3).

--------------------------------------------------------------------------------
## (3) FIGURES + APPEAL EXPLANATION (docs/*.png)
--------------------------------------------------------------------------------
APPEAL-RELEVANT (the honest story):
- `paireig_8x8_doped.png` — pairing eigenvalue: the leading-channel test.
- `cpqmc_vertex_map_8x8.png` / `dqmc_vertex_map_6x6_b2.png` — d-wave VERTEX (corr+susc, 4
  reductions) over (t',tam): tam enhances peak-q equal-time corr, suppresses q=0.
- `cpqmc_full_map_8x8.png` — FULL (bubble-contaminated) for contrast.
- `dqmc_4x4_contradiction.png` — the susc-up/corr-down "contradiction" IS the bubble.
- `iso_vs_aniso_tprime.png` — only the (subleading) effects of t1 vs tp; t' does not enhance.
- `xval_three_method_t1.png`, `xval_dqmc_cpqmc_{t1,tam}.png` — cross-method agreement on trends.
- `dqmc_qstar_tam_t1.png` — peak-momentum q* of P(q): tam pushes the FULL peak off q=0 (bubble).
- `dqmc_chi_beta.png`, `dqmc_btrend_half_8x8.png` — beta-convergence (susc stabilization payoff).
- `cpqmc_landscape_heatmap.png`, `dqmc_landscape_heatmap.png` — (t',tam) susc/corr landscapes.
APPEAL FRAMING: the defensible claim is NEGATIVE/corrective -- "d-wave is the dominant pairing
channel near half-filling, but the altermagnet (tam) enhances only short-range, incoherent,
non-growing local d-wave AMPLITUDE and SUPPRESSES the coherent (susceptibility) channel; net
it works against d-wave SC." This CORRECTS a mean-field 'enhancement' claim (QMC adds the
beyond-MF coherence suppression). PRB/Comment-level, NOT a positive PRL discovery. Full
narrative in `docs/APPEAL_FINDINGS.md`.

--------------------------------------------------------------------------------
## (4) FURTHER DEVELOPMENT SUGGESTIONS
--------------------------------------------------------------------------------
Physics:
- FINISH the SC-dome doping test (n~0.87, L=14/16, CPQMC tau-integrated eigenvalue): does
  d-wave lead in BOTH the projected susc AND the eigenvalue where SC actually lives?
  (Launched on nodes 250/251; pull cpqmc_nearhalf_L{14,16}.csv.) Make-or-break for any claim.
- If chasing PDW: track chi_d(q) PEAK vs beta at fixed tam>0 -- need it to GROW, not just be
  finite-q. So far it does NOT grow (PDW route closed at L<=8). Larger L + lower T to confirm.
- CP-DQMC pairing eigenvalue at the SC-dome doping with sign=1 (the doped regime DQMC can't
  reach) -- the cleanest positive-or-negative SC statement.
Method:
- Build the Fortran (gfortran INTEGER-kind shims) for an independent cross-check of the python
  ports at manuscript scale; or port CP-DQMC pieces to Fortran for L=12+.
- Proper BSE (irreducible vertex Gamma = chi0^-1 - chi^-1) for the pairing eigenvalue, instead
  of the connected-pair-matrix proxy -- removes the full-bubble cancellation entirely.
- Distributed driver / use clusterq broker (server_code/clusterq) instead of manual per-node
  ssh (home is not shared; deploy code per node, watch per-node python paths).
- CPQMC longer-bp convergence is variance-limited; an importance-sampled/longer-projection BP
  or the AGP open-shell trial would let CPQMC reach the true T=0 susceptibility at half-filling.

================================================================================
# (history below) MANUSCRIPT REVISION WORK — d-wave vertex, k-resolution, run DB (2026-06-27)
================================================================================

## MANUSCRIPT REVISION WORK — d-wave vertex, k-resolution, run DB (2026-06-27, read first)

**Context.** The PRL appeal manuscript (Fortran CPQMC, data on cluster
`amax@node-253:~/run/CPQMC/dopingAM/`, plot in `plot_appeal.ipynb`) reports the d-wave
pairing **vertex** (`Vertex_Nppair.dat`, k-resolved → real-space `N_pp(R)`) INCREASING
with the altermagnet anisotropy `tam`, scanning DOPING (n=0.3..1.0) × tam (0..0.9) × U
(0..12) × L (6..22), at U=4, closed-shell fillings near half (e.g. n=0.984).

**THE APPARENT CONTRADICTION — RESOLVED (it was the OBSERVABLE/projection, no bug):**
our first pyqmc tam-scan (q=0 SUM of the connected vertex, U=6, exact half-filling =
open shell) showed tam SUPPRESSING d-wave. The manuscript plots the **local R=0** vertex
`N_pp(R=0)` at U=4 closed-shell, which RISES with tam. Different projections of the SAME
vertex: `P(q=0)=Σ_R N_pp(R)` (q=0 sum, falls) vs `N_pp(R=0)=Σ_q P(q)/N` (k-average,
rises). Verified bond-by-bond that pyqmc `am_hopping` == the Fortran `tk` (so no model
bug); the hopping SIGN is a sublattice gauge (irrelevant for NN-only / t1=0); my exact
half-filling = open-shell = the CP-bias-prone regime (their closed-shell n=0.984 is the
clean one). **pyqmc reproduces the manuscript's increase** once observable + U=4 +
closed-shell are matched (6x6 N17: N_pp(R=0) 0.10→0.15 with tam; U=0 gate = 0).

**NEW TOOLS (this session, all on master, cpqmc.py regression bit-for-bit):**
- `pyqmc/pair_rspace.py` — R-resolved equal-time d-wave vertex N_pp(R): local R=0 +
  long-range <N_pp>_{|R|>2}. (superseded by unified_scan.)
- `pyqmc/unified_scan.py` — **ONE CPMC walk → ALL observables**: energy, d-wave (and
  s-wave) CORRELATION (equal-time) and SUSCEPTIBILITY (τ-integral), each FULL and
  connected-VERTEX, fully **k-resolved** P(q). Outputs 4 reductions per quantity:
  **maxk** (peak in k), **k0** (q=0), **r0** (R=0 local), **rgt** (<|R|>2). CSV cols
  (0-idx): 14-17 corrV maxk/k0/r0/rgt, 22-25 suscV. KEY FINDING: the d-wave vertex
  PEAKS AT FINITE q (maxk ≫ k0). Gates: U=0 all vtx reductions=0; k0/r0/rgt reproduce
  the separate estimators. Has `--t1` (spin-dep/anisotropic NNN) and `--tp`.
- `am_hopping(... t1, tp)` — `tp` = spin-INDEPENDENT (isotropic) NNN, to separate the
  ANISOTROPIC altermagnet t' (t1) from a general t' (tp). Regression-safe (tp=0 default).
- `pyqmc/run_db.py` + `results/runs.db` (force-tracked) — **provenance registry**:
  campaigns (driver + git commit + cluster node + path + observables + status) and
  per-parameter points. `run_db.py find --lx 8 --U 4 --tam 0.3 --t1 0` answers "is this
  point done & where". 8 campaigns / ~1000 points logged. Export: docs/run_registry.md.
- Plotters: plot_rspace_tam, plot_pd_heatmap, plot_tt_heatmap, plot_cs_heatmap,
  plot_un_reductions. Figures in docs/: rspace_vertex_vs_tam, pd_heatmap_Cd0 (n×tam
  phase diagram), tt_heatmap_Cd0 ((t1,tam) 2 dopings), cs_heatmaps (corr+susc), and the
  un_t1_{maxk,k0,r0,rgt} four-figure set (pending plot).

**CLUSTER COMPUTE** (see memory [[cluster-compute-access]]): `ssh -p <node#>
profhokin@thkclusters.duckdns.org` (port=node#; my env on node-253/258; nodes 254/256
have a broken py3.8 = no distutils, skip). Manuscript data: `amax@...:~/run/CPQMC/
dopingAM` (READ-ONLY; different acct). amaz 48-core workstation still available too.

**KEY PHYSICS RESULTS (free trial, ED-checked estimator):**
1. Manuscript's local `C_d(0)=N_pp(R=0)` rises with tam — REPRODUCED.
2. d-wave vertex strongest NEAR HALF-FILLING (n×tam phase diagram, pd_heatmap), fades
   with doping; peaks at finite q.
3. q=0-sum vertex (correlation & susceptibility, same structure) FALLS with tam — the
   complementary projection.

**PLAN / IN PROGRESS (the 3-part manuscript-revision request):**
- (1) k-resolved q-grid output (corr/susc, full/vtx) — **DONE** (unified_scan).
- (2) **isotropic t' (tp) vs anisotropic t1**: claim only the ANISOTROPIC t' enhances
  d-wave. Code done; RUN the (tam,tp) grid (8x8 U4, two dopings) and compare to (tam,t1).
  NOT yet run (nodes were busy with the t1 grid).
- (3) **FOUR figures** (maxk/k0/r0/rgt), each 2x2 (corr/susc × n=0.97,0.81) over
  (tam,t1): (tam,t1) grid DONE on cluster (un_t1, node 258 N31 / 253 N26, U4, 36 pts
  each) → results/un_t1; NEXT: collect + `plot_un_reductions.py` → docs/un_t1_*.png.
NEXT STEPS: collect un_t1 → 4 figures; launch (tam,tp) grid → isotropic comparison;
ingest both into run_db; optional L16 n=0.984 direct overlay on their Vertex_Nppair.dat.
Honest scope: free-trial CP estimate (CP-bias caveat); report trends; the connected
vertex is the SC-relevant quantity; AGP/BCS trial helps only on small clusters (unusable
≥6x6) — large-lattice numbers use the validated free trial.

## CURRENT STATE — read first (2026-06-24)
Platform validates a Constrained-Path QMC against QuSpin ED. Fortran code/src
(two-orbital altermagnet, the PRL LE20050 model) validated vs ED on all three
interaction channels (uxx/uxy/v). Python port `pyqmc/cpqmc.py` is the active R&D
surface; every change is gated bit-for-bit by `pyqmc/regression.py` (5 fixed-seed
configs) — run it after any edit. Repo is on `master`, remote `origin` live, all
pushed. Tag `v1.0` = state before the multidet work. Theory write-up:
`docs/theory/cpqmc_theory.pdf`. Plan: `docs/PLAN_phase6.md`. Quantitative log:
`docs/VALIDATION.md` (the source of truth for results).

**pyqmc OOP structure (Phase 6 WS1a):** LatticeModel / TrialWF / WalkerEnsemble /
Propagator / Estimators / CPMC. Estimators implemented + ED-validated: energy
(mixed + back-propagated), equal-time correlations, s/d-wave pairing + magnetization
(`run_bp_pairmag`), unequal-time staggered spin susceptibility (`run_bp_chi`), and
the unequal-time singlet PAIRING susceptibility chi_s/chi_d (`run_bp_chid`,
`pyqmc/validate_chid.py` — ED Lehmann on the full Fock space). All exact at U=0;
interacting agrees within the expected constrained-path bias. The dynamic d-wave
pairing susceptibility (the PRL appeal's order parameter) is now covered.

**Trial wavefunctions (Phase 6 WS3):**
- `fixed` (free-electron GS) and `adaptive` (self-consistent natural-orbital from
  the walker 1-RDM) — both VALIDATED, the production options. Adaptive reduces the
  CP bias vs ED (~17-22%) and its bias DECREASES with walker count (fixed's does
  not). Moderate damping (trial_mix=0.5) beats aggressive.
- `sample1` (single weight-sampled walker) — tested, NOT recommended (non-
  variational, huge variance).
- `multidet` (weighted walker superposition) — tested, UNSTABLE as a method
  (walker superposition is ill-conditioned: W->0). BUT the multidet ENGINE is now
  correct + validated: feeding it the exact GS as a truncated CI
  (`pyqmc/validate_casci_ed.py`) drives the CP bias monotonically to 0 and full-N
  gives E0 with zero variance. Fixing it required a real two-body-energy bug fix
  (`Estimators._local_energy`: per-determinant Wick average, not Wick of averaged G).

**Scalable multidet trial = self-consistent natural-orbital CASCI (DONE 2026-06-24).**
`pyqmc/casci.py` (general Slater-Condon CI in an arbitrary orbital basis +
integral transform + active-space builder; basis-invariant E0=ED to 1e-14) +
`pyqmc/validate_casscf.py` (the SCF loop: walker 1-RDM -> natural orbitals ->
freeze core/CAS/drop virtual -> CASCI -> install frozen multidet -> iterate).
2x2/4+4/U=4: full active = ED exactly (gate); reduced CAS(2,2) cuts the CP bias
~85x (1.548->0.018) with 4 configs and NO ED — the trial NODES, not its energy,
drive the bias to ~0. Mini-CASSCF complete.

**SINGLE-BAND spin-dependent altermagnet (the MANUSCRIPT'S primary model) — DONE
2026-06-25.** The manuscript (PRL LE20050) is NOT the two-orbital model; it is a
single-band square Hubbard with SPIN-DEPENDENT anisotropic hopping (production code
= `code/record/benchmark/CPQMC`, not `code/src`). Built it in pyqmc:
- `am_hopping(lx,ly,t0,tam,t1)`: NN `tam` (K_up: weak-x/strong-y; K_dn rotated 90)
  AND NNN `t1` (spin-dependent d_xy diagonal). Spin-dependent hopping K_up!=K_dn
  added to LatticeModel/CPMC (`K_dn=`; regression bit-for-bit when K_dn=K).
- VERTEX (connected) pairing = full - disconnected, equal-time (`run_bp_pairmag`)
  AND susceptibility (`run_bp_chid` -> chi_a + chi_a_vertex), in pyqmc AND ED
  (`validate_am_single.py`, `validate_chid.py`). ED-validated: U=0 vertex=0 exactly;
  dilute interacting agrees within error.
KEY FINDINGS (`docs/single_band_altermagnet_report.md`):
1. The anisotropy that reproduces the manuscript trend (AFM down, d-wave vertex up)
   is the NNN `t1` (it BREAKS (pi,pi) nesting); NN `tam` alone PRESERVES nesting and
   gives the OPPOSITE trend. State the anisotropy type precisely in the reply.
2. ED and QMC MATCH for energy/s-wave/spin/susceptibility/U=0/dilute. The ONE
   exception is the **d-wave VERTEX at half-filling + anisotropy**: severe
   constrained-path bias (4x2: ED +23 vs CPMC -2, sign-flipped). ED supports the
   enhancement; CPMC d-wave-vertex magnitude/sign is unreliable in that regime.
   adaptive trial improves AFM but NOT the d-wave vertex. -> CP-bias reduction is the
   key open problem for a quantitative d-wave result (see options below).

**d-wave PAIRING SUSCEPTIBILITY SCAN — DONE 2026-06-25** (`docs/pairing_susceptibility_report.md`,
commit d0cde44). Ran the dynamic/tau-integrated singlet pairing susceptibility
chi_a=integral C_a(tau)dtau (full + connected VERTEX, d and s wave) vs the
nesting-breaking NNN `t1`, half filling, U=6, tam=0, across three lattices on all
available compute: 4x4 (laptop, 4 seeds), 6x6 (workstation amaz 48-core, 6 seeds),
8x8 (cluster node-250, 4 seeds). Driver `pyqmc/chid_tam_scan.py` (--t1 added; emits
chi_d/chi_d_vtx/chi_s/chi_s_vtx); one-process-per-point, xargs -P. CSVs under
`results/chid_scan/` (gitignored; all numbers are in the report).
FINDING: the NNN `t1` **selectively enhances the d-wave susceptibility** (chi_d full
AND vertex) over s-wave (flat). On the larger lattices it's a **DOME** — chi_d^vtx
peaks at intermediate `t1≈0.2` (6x6: 1.67->1.99->1.01; 8x8: 2.60->3.69(+42%)->2.09)
then is suppressed by large t1; 4x4 pushes the peak to t1≈0.4-0.5 (finite-size). So
there's an OPTIMAL anisotropy. METHOD POINT for the reply: the tau-integrated VERTEX
stays POSITIVE and is ED-validated within error bars on clean dilute clusters (d-wave
ED 0.0168 / QMC 0.0246, ~1.6sigma) — unlike the sign-flipped EQUAL-TIME vertex
(ED +23 / CPMC -2). Report the d-wave result as a TREND of the susceptibility, with
the CP-bias caveat on absolute magnitude. NOTE: this scan used the **fixed
free-fermion trial** (CPMC default, the open-shell single determinant); adaptive/CASCI
NOT used here — so the d-wave-vertex bias above is exactly this trial's node bias.

**EVOLUTIONARY / NEURAL TRIAL-STATE TRACK — A(a) DONE & ED-VALIDATED, A(b) WORKING
2026-06-25** (new workspace `qmc_algo/`, design `qmc_algo/notes/00_design.md`,
findings `qmc_algo/notes/01_findings_A_a.md`). A separate R&D line on the owner's
idea: use EVOLUTIONARY ALGORITHMS / NEURAL QUANTUM STATES to generate a better CPMC
trial (attacks the same d-wave CP-bias blocker, as a general method). Key framing:
CPMC population control IS already an evolutionary algorithm (branch=select,
auxiliary-field=mutate); the missing pieces are CROSSOVER ("walkers as parents")
and a richer GENOTYPE (neural state). Crossover of two determinants is natively a
MULTIDETERMINANT state — exactly what the repo proved kills the CP bias (the CASCI
result); the open problem "good multidets without ED" is a SEARCH problem.
- **A(a) — evolutionary selected-CI trial** (`qmc_algo/src/evo_trial.py`): a genetic
  algorithm over occupation-string configs, scored by Epstein-Nesbet marginal energy
  (CIPSI criterion) + memetic local single-excitation enumeration, NO ED. Builds a
  multideterminant trial for `TrialWF.set_multidet`. Reuses `pyqmc/casci.py`
  Slater-Condon helpers. TWO-POPULATION FIREWALL: estimators from the CPMC walkers
  (unbiased), the evolved archive only defines the constraint.
- **ED gates** (`qmc_algo/tests/validate_evo_trial.py`): Gate 1 exactness = full-space
  evolution reproduces ED to 1.6e-14; Gate 2 = CP bias falls monotonically to 0
  (2x2 half-filled U=4, the maximally-multireference worst case: free-det +1.49 ->
  +0.011 at k=2 -> 0 by k=16). The GA now matches the ED-optimal (oracle) config
  ranking at EVERY determinant count.
- **HEAD-TO-HEAD vs adaptive & CASCI** (`qmc_algo/tests/compare_trials.py`,
  `qmc_algo/results/compare_2x2_U4.txt`): at MATCHED determinant count the evolutionary
  trial WINS — ~5x lower CP bias than adaptive at 1 det (0.115 vs 0.548), ~15x vs
  CASCI at 4 det, ~900x vs CASCI at 9 det. Mechanism: selected-CI spends each
  determinant on the individually most-important config; CASCI wastes them filling a
  whole active window; adaptive is capped at one determinant.
- **A(b) — self-consistent WALKER INJECTION** (`qmc_algo/src/evo_trial_scf.py`): the
  literal "walkers as parents" — mine the live walker population for its dominant
  configs, inject as genes, re-evolve; + a PAIRING SEED hook (paired double
  excitations) for the d-wave target. Runs end-to-end (2x2: near-exact, injects
  walker configs each round). Validated demonstration on a LARGER/anisotropic case
  (where injection matters and full enumeration is infeasible) is the next step.
- THREE bugs found+fixed (notes/01): candidate-gen livelock; sign-flipped EN
  selection (rejects the key configs in a multireference state -> select by MAGNITUDE);
  basis-orthogonality nan (re-init walkers on the trial's leading determinant).
- pyqmc regression untouched (bit-for-bit green). Relates to [[qmc-platform-plan]];
  complementary to the AGP pairing trial below (that hard-codes pairing nodes; this
  discovers multireference/pairing configs by search).
**Open / next (evolutionary track):** larger-cluster scaling gate (GA explores <<
full space); A(b) on the anisotropic half-filled d-wave case with pairing seed +
walker injection vs the AGP trial; then idea-2 neural generative gene pool (B1).

**PAIRING-AWARE (AGP / number-projected-BCS) TRIAL — core DONE & ED-VALIDATED
2026-06-25** (`docs/agp_pairing_trial_report.md`, commits 5bd5d07 + fe76194). The
diagnosed fix for the d-wave vertex CP bias is a pairing trial; built it in pyqmc,
referencing the sister repo `thkaitools-cell/pph-qmc` (Nambu/BdG effort; see
[[pph-qmc-reference]]). The hard, rigorous part is complete:
- `pyqmc/agp.py`: AGP trial `|Ψ⟩=(Σ F_ij c†_i↑c†_j↓)^N|0⟩`, overlap with a
  number-conserving Slater walker = `det(Aᵀ F B)`. KEY PHYSICS: the projected AGP is
  NON-Gaussian, so normal-only Wick is WRONG; but since H conserves N and the walker
  is fixed-N, `⟨Ψ_AGP|H|Φ⟩/⟨Ψ_AGP|Φ⟩ = ⟨BCS|H|Φ⟩/⟨BCS|Φ⟩` on the UNPROJECTED Gaussian
  BCS, where generalized Wick is exact INCLUDING the anomalous pairing tensors
  `Kd=⟨c†c†⟩=F(I−Gdᵀ)`, `Ka=⟨cc⟩=−A M⁻ᵀ Bᵀ` (M=AᵀFB). Cross-spin term gets
  `−Kd·Ka`; this is exactly what pph-qmc's normal-only Nambu prototype (their §6
  bottleneck) was missing.
- `pyqmc/validate_agp.py`: brute-force QuSpin ED gate — overlap/Gu/Gd/Kd/Ka and the
  normal+anomalous local energy match to ~1e-14 (N=1,2 / half-filling / t1 / tam).
- d-wave geminal from BdG `F=−(U†)⁻¹V†` is SINGULAR (cond~1e16): d-wave NODES make
  deep orbitals unpaired (v/u→∞). USABLE construction = `augmented_geminal`
  `F = A Bᵀ + η Fd` (η=0 → free trial exactly, η>0 → bounded d-wave pairing).
- `pyqmc/agp_dwave.py`: CP-AFQMC with the AGP trial. VALIDATED: U=0 → energy exact AND
  d-wave VERTEX = 0 exactly; η=0 → reproduces free-trial CP-AFQMC energy (=ED on clean
  clusters, e.g. 2×2/1+1/U=4/t1=0.3: −7.256(4) vs ED −7.254). cpqmc.py untouched
  (regression bit-for-bit).
**BP-AGP BRA → 7x d-wave VERTEX BIAS REDUCTION — DONE 2026-06-25** (plan
`docs/agp_vertex_plan.md`; code `pyqmc/agp_bp_vertex.py` + full-Fock ED gate in
`agp_dwave.ed_reference`; commit d9af219). Decouple the trial's two roles: STABLE
free-electron determinant as the CP constraint/importance (AGP-as-constraint is
variance-unstable at half filling), + a BACK-PROPAGATED AGP (BCS) bra for the
measurement only. BP of the AGP bra = geminal transform F_bp = Bu_tot^T F Bd_tot from
the recorded forward propagators (since <AGP(F)|(prod B)|phi0>=det(phi0^T F_bp phi0)).
Connected vertex subtracts BOTH disconnected pieces (normal bubble + anomalous
<Delta><Delta^dag>) -> U=0 vertex=0 exactly at every eta (gate). RESULT (2x2/2+2/U=4/
t1=0.3; full-Fock ED vertex +22.41): free bra 2.3 (10% of ED) -> BP-AGP bra
17.0+/-0.65 at eta~0.5 (76% of ED), a 7x cut of the d-wave-vertex CP discrepancy,
stable plateau eta in [0.3,0.6]. The residual ~24% tracks the free-trial KET-constraint
bias (eta=0.5 BP energy -5.21 vs ED -5.657). NOTE qmc_algo/ (evolutionary multidet
trial, separate session) left untouched -- this is the standalone BCS thread.
RESIDUAL ATTACKED (2026-06-25, commit 175cacd): bp-length scan shows the vertex
SATURATES ~18 = 80% of ED by bp~32 (14.8->18.1->17.9 at bp 16/32/48), i.e. longer
back-propagation recovers 76%->80%; the saturation ceiling is the free-trial KET
constraint bias. Constraint RELEASE (agp_release_vertex.py: CP equil -> release ket +
fixed AGP bra; U=0 vertex=0 gate) makes the released vertex transiently REACH ED
(~22-30 at tau~0.14-0.42) but is VARIANCE-LIMITED (the <AGP|phi> denominator blows up
for released walkers; BP regularizes, direct mixed does not) -> inconclusive on the
pathological 2x2-degenerate cluster. NET: vertex cut 10%->80% of ED robustly; last
~20% is ket bias.
OPEN to close it: a STABLE improved constraint (variance-controlled AGP/BCS importance
fn) or a less-degenerate cluster (free constraint already good); push to the 4x2
target with a SPARSE full-Fock ED reference (dense 4^8 is 34GB; scipy sparse matvecs).

**Open / next:** finish the pairing-trial application above (BP AGP estimator + η
control + full-Fock ED pairing ref) for a quantitative d-wave-vertex bias reduction;
find production t0/tam/t1/U/filling (no in.dat in submission records); 400-core
large-cluster susceptibility scan for the TREND; Fortran StepMeas bug. The detailed
dated entries below are the full history; this block is the summary.

## Goal
Build a platform to validate changes to a Constrained-Path QMC (CPQMC) code by
comparing its results against (a) exact diagonalization (ED) ground truth and
(b) other QMC versions. Then port the QMC to Python. Owner is a condensed-matter
/ QMC physicist + scientific-software engineer.

## What the code is
- **CPQMC**: ground-state projector auxiliary-field QMC for 2D Hubbard-type
  lattice models (+ optional Holstein e-phonon), MPI-parallel over walkers.
- **DQMC (BSS)**: finite-T determinant QMC, used as a secondary cross-check.

## Repository layout
```
code/record/                 FROZEN baselines — never edited
  benchmark/CPQMC            validated single-band CPQMC (reference)
  benchmark/DQMC             BSS finite-T DQMC (cross-check)
  revise/250525, 251111      historical revision snapshots
code/src/                    ACTIVE dev tree (seeded from 251111); all edits here
ed/hubbard_ed.py             QuSpin ED reference (common schema)
tools/                       build.sh, run.sh, parse_out.py, compare.py, env.sh
docs/                        ENVIRONMENT.md, VALIDATION.md, HANDOFF.md
environment.yml              conda spec (qmc-platform)
build/  results/             gitignored (out-of-source builds / run outputs)
.tools/                      gitignored per-machine micromamba binary
```
Milestones are git **tags**, not new dated folders. `record/` stays byte-frozen.

## Environment (Phase 0 — DONE)
- conda env `qmc-platform` (Python 3.11) via the bundled `.tools/bin/micromamba`
  (classic conda solver too slow). Recreate: see `docs/ENVIRONMENT.md`.
- Supplies the missing system MPI: `openmpi` + `gfortran` 14.3 + `openblas`.
- QuSpin installed via **pip** (no conda-forge arm64 build) — native arm64.
- **Always** `source tools/env.sh` before building (sets SDKROOT for gfortran).

## Validation loop (Phases 1 & 2 — DONE)
```bash
source tools/env.sh
tools/build.sh benchmark-cpqmc
tools/run.sh   benchmark-cpqmc code/record/benchmark/CPQMC/in.dat 4
python tools/parse_out.py results/benchmark-cpqmc/<runid> -o qmc.json
python ed/hubbard_ed.py --lx 4 --ly 4 --nup 1 --ndn 1 --t 1 --U 3 -o ed.json
python tools/compare.py ed.json qmc.json
```
- **Common schema** (the contract ED + QMC + UI share): `{code, run_id,
  model{...}, observables{energy_total/kinetic/potential/per_site: {value,error}},
  kspace, correlations}`.
- `compare.py` reports delta / combined-sigma / z-distance / PASS-FAIL.

### Established result
4x4 Hubbard, Nup=Ndn=1, U=3, t=1: **total energy ED -7.8672 vs CPQMC -7.8701,
z=1.01 -> PASS**. Single-band CPQMC validated against ED. PE component fails
(z=15.7) = expected projector-QMC mixed-estimator bias for operators not
commuting with H; total energy is the primary observable. See VALIDATION.md.

## Operational gotchas (don't re-derive)
- gfortran on arm64: **no `-mcmodel=large`** (drop it).
- Run needs `ulimit -s 65520` + `OMP_STACKSIZE` (large arrays; cluster used
  -mcmodel=large) — handled inside `tools/run.sh`.
- ED needs single-thread + `KMP_DUPLICATE_LIB_OK` (handled in hubbard_ed.py).
- Frozen Intel sources need gfortran shims (char-array constructors etc.),
  applied to staged COPIES in `tools/build.sh` (record/ untouched).

## UI (Phase 3 — DONE)
FastAPI backend (`platform/backend/app.py`) + React/Vite frontend
(`platform/frontend`). Launch: `source tools/env.sh && tools/serve.sh`,
open http://localhost:5180. Backend: /api/versions, /api/runs,
/api/runs/{code}/{run_id}, /api/ed (isolated subprocess), /api/compare.
Frontend: ED form, run picker, energy chart (error bars) + delta/sigma/z/
PASS-FAIL table. See platform/README.md. (Vite on 5180 strictPort to avoid
an IPv4/IPv6 collision with another local dev server.)

## code/src build + physics (Phase 4 — in progress)
`code/src` is the PRL LE20050 altermagnet model: two-orbital (d_xz,d_yz)
square lattice, eps_x=-2t1 cos kx -2t2 cos ky -4t3 cos kx cos ky,
eps_y=(x<->y), eps_xy=-4t4 sin kx sin ky; on-site uxx, inter-orbital uxy,
neighbour v (signs: +v intra-orbital, -v inter-orbital).

Builds + runs via `tools/build.sh src` / `tools/run.sh src`.

VALIDATED (U=0): ed/altermagnet_ed.py free-fermion GS (analytic-verified) =
-96.0 for 4x4/16+16/t=-1. code/src now reproduces this exactly after fixing
a back-propagation Green's-function sign bug in calgf method 2 (commit
ef171d9; it built gx=+G but calenergy wants gx=detp*(I-G) like method 1).
The U=0 measured energy went +96 -> -96.

INTERACTING VALIDATED (on-site uxx channel): `tools/build.sh src-small`
(patches cpParameter to lx=ly=2, NUP=NDN=4 so ED is tractable). Compared
out.dat `totalEn` (the production energy) to interacting ED
(ed/altermagnet_ed.py ed_energy) on 2x2/4+4 with equilibration:
  uxx=0: -32.000 vs -32.000 (z=0)
  uxx=2: -28.275+/-0.007 vs -28.289 (z=1.85 PASS)
  uxx=4: -25.134+/-0.010 vs -25.128 (z=0.58 PASS)
The interacting two-orbital altermagnet CPQMC reproduces ED. USE out.dat
`totalEn` as the production energy (parse_out.py). The separately PRINTED
"energy is the:" number (StepMeas mixed estimator) is still wrong (e.g.
-19.98 at uxx=2) -- a secondary estimator-path bug, distinct from the
validated observable.

NEXT: validate uxy and v channels the same way; fix the StepMeas mixed
estimator; extend parse_out.py to read out.dat for src (multi-orbital);
then scale toward publication clusters.

## Known gaps / next
- **Phase 4 (continue)**: interacting validation via small-lattice ED (above).
- **Phase 5 (STARTED)**: `pyqmc/cpqmc.py` = from-scratch Python constrained-
  path AFQMC. Single-band Hubbard validated vs ED (U=3 -> z=1.66 PASS).
  EXTENDED to the two-orbital altermagnet (`--model altermagnet`, reuses
  ed/altermagnet_ed.build_hopping). THREE-WAY validation (2x2, 4+4) ED vs
  Fortran code/src vs Python:
    uxx=0: -32.000 / -32.000 / -32.000
    uxx=2: -28.289 / -28.275(7) / -28.289(1)
    uxx=4: -25.128 / -25.134(10) / -25.126(5)
  All three agree. NEXT: uxy and v interaction channels in pyqmc + ED + code.

## UI is parameter-driven (runs QMC live)
The UI Model&parameters panel drives both Compute ED and Run QMC for either
model. Backend: POST /api/qmc runs the Python CPMC port (subprocess) for the
given params; POST /api/ed dispatches hubbard vs altermagnet. Compare overlays
ED vs the live QMC run. NOTE: the Python port models on-site U (uxx) only; if
uxy/v != 0 the QMC response carries a red note (those channels are ED-only).

## ALL THREE CHANNELS VALIDATED (2026-06-24)
uxx, uxy, AND v all reproduce ED on small clusters:
  uxx=2: z~0.25 ; uxy=2: z~1.33 ; v (anisotropic 2x2): v=0.1 z=0.95,
  v=0.2 z=0.39, v=0.3 z=1.08 -- all PASS.
Key fixes: (a) backphi replayed each interaction term's 2nd operator at the
WRONG site (loop index instead of bond partner i+lxy/neighbour) -> rewrote to
mirror Vxy/Vznn bond structure; (b) rec_fields v-gate abs(uxy)->abs(v);
(c) InitV per-channel V sign; (d) meas.f90 energy estimator includes uxy/-v.
The v "collapse" on the SYMMETRIC cluster (t1=t2=t3=t4) was degeneracy +
level-crossing + poor free-fermion trial (a hard AFQMC regime), NOT a bug --
use anisotropic/non-degenerate clusters (build.sh src-small honours
QMC_LX/LY/NUP/NDN). 2-site reference: pyqmc/vee_check.py.
pyqmc now models uxx+uxy+v (general density-density HS); mixed estimator so v biases at larger coupling.

## uxy NOW VALIDATED (2026-06-24)
The 2-site reference pyqmc/vee_check.py proved Vee's forward determinant/GF
update is exact (1e-16), isolating the uxy bias to BACKPHI: it replayed each
term's 2nd operator at site j(=loop index) instead of the bond partner
(i+lxy for uxy; neighbour for v) -- correct only for same-site uxx. Fixed
backphi to mirror the Vxy/Vznn bond structure; fixed rec_fields v-gate
(abs(uxy)->abs(v)). Now: uxx=2 z=0.25 PASS (no regr), uxy=2 -24.614 vs ED
-24.601 z=1.33 PASS. ED<->Fortran validated for uxx AND uxy.
REMAINING: v collapses walkers (phi->0 ~step 41) -- a PRE-EXISTING
forward-propagation instability in the attractive -v channels (MkExpV HS
constants verified correct; suspect Vznn importance sampling / constraint /
stabilization). Separate from the fixed back-prop. pyqmc still uxx-only.

## uxy / v channel status (UPDATED 2026-06-24)
Two real bugs fixed (no uxx regression): (1) meas.f90 production energy
estimator now includes uxy + inter-orbital -v terms (correct Wick form);
(2) InitV per-channel V sign (was overwriting all 32 v-channels with one sign).
REMAINING bug is isolated to the Vee two-site HS PROPAGATION (i/=j path):
uxy=2 samples -25.41 vs ED -24.60 (below ground state); v=1 collapses walkers.
Next: hand-derive a single Vee update on a 2-site toy and compare; inspect the
same-spin rdet/GF update and the attractive-channel (-v) coeffv/constraint.
See VALIDATION.md for the full trace.

## uxy / v channel status (earlier)
Empirically confirmed the inter-orbital/neighbour interactions are applied in
Fortran PROPAGATION (Step calls Vxy/Vznn) but OMITTED from the energy
ESTIMATOR: calenergy's two-body term sums only uxx (Vlist(i,1)). Test
(2x2,4+4,uxy=2): Fortran totalEn=-32.88 vs ED=-24.60 (off by the missing
inter-orbital potential ~+8.3). ED supports uxy/v (ed/altermagnet_ed.ed_energy,
validated U=0). To validate uxy/v end to end:
  (1) Fortran: extend calenergy two-body term to inter-orbital (uxy) and
      neighbour (v) density-density with the Wick exchange piece
      <n_a n_b> = <n_a><n_b> - sum_sigma G_sigma(a,b) G_sigma(b,a),
      in the gx = detp*(I-G) normalization (mirror the existing uxx term).
  (2) Python pyqmc: add uxy/v via charge-channel HS (heavier; sign problem) --
      currently uxx only.

## Report & launching the UI
- Summary PDF: `python docs/report/generate_report.py` -> docs/report/
  validation_report.pdf (2 pages: platform, single-band check, bug fix,
  two-orbital three-way table+chart).
- Launch platform: `source tools/env.sh && tools/serve.sh` ->
  http://localhost:5180 (frontend; backend on :8000). vite uses port 5180
  (another local project holds 5173).
- `benchmark-dqmc`: gfortran rejects INTEGER(8)/(4) kind mismatches in
  BSS.f90 — deferred (secondary cross-check).

## Git
Local repo on `master`. One commit per logical change, messages end with the
Co-Authored-By trailer. **No remote yet** — run `gh auth login`, then a remote
will be added and history pushed.

## Phase 5 + platform additions (2026-06-24)
- pyqmc back-propagation: `run_bp()` / CLI `--bp <len>`. Removes the mixed-
  estimator bias (v=0.2: mixed -17.034 z=5.3 -> back-prop -17.050 z=0.15 vs ED).
- Platform more functional: /api/qmc gains `bp`; new /api/sweep returns the ED
  curve (+ optional QMC points) over U/uxy/v; UI has a bp input and a Parameter
  sweep panel with an E-vs-interaction line chart (ED line + QMC points).
  Launch: source tools/env.sh && tools/serve.sh -> localhost:5180.

## Speed + correlation benchmarks (2026-06-24)
Two new platform features (tabs in the UI; endpoints + CLI):

(1) SPEED — Fortran code/src vs Python pyqmc on the same small problem.
  - tools/benchmark.py ; backend POST /api/benchmark ; UI "Speed" tab.
  - Reports wall time + throughput (walker-updates/s) + speedup. NOTE small
    runs include MPI startup (~seconds) which UNDERSTATES the Fortran
    throughput; use larger --fort-nblk/--fort-nblkstps for a fair ratio.
    Fortran is compiled+MPI; Python is interpreted numpy with brute-force
    determinant ratios (so Python scales worse with system size).

(2) CORRELATIONS — equal-time GF + spin/charge, ED vs back-propagated CPQMC.
  - ED: ed/altermagnet_ed.py ed_correlations / `--corr` -> G^s_ij=<c+_is c_js>,
    charge <n_i n_j>, spin <S^z_i S^z_j> (densities from diagonal number ops;
    off-diagonal G from hopping ops; ~17s on 2x2).
  - CPQMC: pyqmc run_bp_obs / `--corr` -> same observables from each walker's
    back-propagated Green function via Wick, weighted-averaged.
  - backend POST /api/correlations compares element-by-element (max|dev|,
    mean, PASS if max<0.05); UI "Correlations" tab shows the summary table +
    GF heatmaps (ED vs CPQMC).
  - RESULT (2x2, 4+4, uxx=2): max|ED-CPQMC| green 0.004, charge 0.009,
    spin 0.002 -> the CPQMC GF and spin/charge correlations reproduce ED.
  - The Fortran out.dat also has correlations (cdw/sdwz/...) but matching its
    conventions to ED is future work; the ED<->Python correlation check is
    the validated one.

## NWLKRS convergence knob (2026-06-24)
`tools/build.sh src-small` now patches NWLKRS via env `QMC_NWLKRS` (default 1000),
alongside QMC_LX/LY/NUP/NDN. Lets ED-comparison runs raise the walker count for
tighter error bars without editing the source. Convergence demo (2x2/4+4, uxx=2,
ED=-28.289): NWLKRS 250/1000/4000 -> error 0.0187/0.0073/0.0033, central value
pinned on ED (~1/sqrt(N), no population bias). Same sweep on all three channels:
  uxx=2 (ED -28.289): -28.279(0.019)/-28.290(0.007)/-28.291(0.003)
  uxy=2 (ED -24.601): -24.580(0.015)/-24.597(0.015)/-24.602(0.004)
  v=0.2 anisotropic (ED -17.053): -17.060(0.013)/-17.053(0.006)/-17.051(0.004)
All channels: bias-free convergence onto ED. (commit d12f801)

## PHASE 6 (PLANNED 2026-06-24) — see docs/PLAN_phase6.md
Three workstreams requested by the owner, plan written down before implementation:
1. **Magnetization + pairing parity** (pyqmc <- Fortran meas.f90) with an OOP
   tidy-up FIRST (split CPMC into Lattice/Model, Walker, TrialWF, Propagator,
   Estimators; regression-gate the refactor on the existing 3-channel energy +
   correlation PASS). Port local moment m_i^z, staggered M_s, S(q), and singlet
   pairing P_alpha(r) for s / extended-s / d_{x^2-y^2} (the appeal's order
   parameter). ED references extended in altermagnet_ed.py.
2. **Dynamical susceptibility** (UNEQUAL-time): imaginary-time-displaced G(tau)
   from the BP machinery -> chi_s(q,tau), chi_d(tau) (dynamical d-wave pairing).
   ED Lehmann-representation reference on 2x2 as the exact gate. New run_bp_chi.
3. **Adaptive trial WF**: switch trial="fixed"(current free-electron GS) vs
   "adaptive" (constraint boundary redefined each step from the walker ensemble,
   e.g. 1-RDM natural orbitals). TEST: does adaptive reduce |E-E_ED| in the
   DOPED / sign-problem regime (away from half-filling) without regressing the
   half-filled PASS cases? Honest caveat: self-consistent trial is non-variational
   and can introduce its own bias; keep FIXED as the validated default.
Built on the TrialWF seam from WS1a. Sequencing + acceptance gates in PLAN_phase6.md.

## PHASE 6 IN PROGRESS (2026-06-24) — pyqmc observables + trial + theory
Plan: docs/PLAN_phase6.md. Quantitative log: docs/VALIDATION.md. All pushed.

**WS1a — OOP refactor (DONE).** pyqmc/cpqmc.py split into LatticeModel / TrialWF /
WalkerEnsemble / Propagator / Estimators + a thin CPMC orchestrator (run/run_bp/
run_bp_obs + ctor unchanged). Gated bit-for-bit by pyqmc/regression.py (5 fixed-
seed configs); identical before/after and after every later change.

**WS1b/1c — pairing + magnetization (DONE).** Estimators.pairmag_block +
CPMC.run_bp_pairmag (single-band square lattice): s-wave + d_{x^2-y^2} singlet
pairing structure factors, spin S(pi,pi), local moment^2, back-propagated with
per-block error bars. ED ref = number-conserving product operator "+-|+-"
(pyqmc/validate_pairmag.py) -- the pair op alone changes particle number and
vanishes in a fixed-Nf basis. EXACT at U=0; interacting shows few-% constrained-
path bias (s-wave +4%). GOTCHA: open-shell degenerate fillings are NOT clean ED
gates (ED vs trial pick different states in the degenerate manifold) -- use
closed-shell non-degenerate (nup=ndn=1) for estimator checks.

**WS3 — adaptive self-consistent trial (DONE).** trial="adaptive" in
TrialWF.update: ensemble 1-RDM -> symmetrize -> damp toward current trial (mix)
-> top-occupation natural orbitals as new trial. CLI --trial/--trial-every/
--trial-mix; frozen across each BP window; fixed mode is a no-op (regression
bit-for-bit). Test (pyqmc/validate_adaptive.py) on 4x4 nup=ndn=2 U=4 (degenerate
open shell, hard CP regime, ED -11.5303): adaptive cuts the bias ~17% (mixed
+0.206->+0.170) / ~22% (bp +0.241->+0.189). Moderate damping (mix=0.5) beats
aggressive (mix=0.25). WALKER SCALING (pyqmc/scan_walkers.py): the FIXED bias is
flat in Nw (~0.21, systematic); the ADAPTIVE bias DECREASES with Nw (+0.228 ->
+0.159 at Nw=500, ~25% better) -- it converts walker count into reduced systematic
bias (coherent effect fixed cannot give). Crossover at low Nw (~100): noisy RDM,
adaptive ~ fixed. adaptive is NON-VARIATIONAL -> empirical gain; FIXED stays the
validated default.

**Theory doc (DONE).** docs/theory/cpqmc_theory.tex (+ compiled cpqmc_theory.pdf):
full CP-AFQMC derivation (projection, Trotter, discrete Hirsch HS matching the
code's _hs, importance sampling + constrained path, mixed + back-propagated
estimators, pairing/magnetization, adaptive trial Eq.), a TikZ algorithm
flowchart, and a theory<->code map table. Build: cd docs/theory &&
pdflatex cpqmc_theory.tex (x3 for refs/TOC). TeX Live 2019 at /Library/TeX/texbin.

**REMAINING Phase 6:** WS2 unequal-time susceptibility chi_s/chi_d (tau-displaced
BP Green functions vs an ED Lehmann reference); altermagnet-geometry pairing /
Fortran-convention match; adaptive refinements (build trial from the back-
propagated RDM instead of the mixed RDM; multi-determinant trial).

## PHASE 6 WS2 DONE — dynamic susceptibility (2026-06-24)
Estimators.chi_block + CPMC.run_bp_chi (single-band): imaginary-time-displaced
staggered spin correlation C(tau_l)=<O(tau_l)O(0)>, O=sum_i (-1)^{x+y} S^z_i, and
static chi_s=integral C dtau. Per-walker forward propagators B_(l) reconstructed
from the recorded HS fields; displaced single-particle GFs P=B(I-g^T) [particle],
H=B^{-T}g [hole], g(tau)=B^{-T}g B^T [equal-time at tau]. ED ref = Lehmann
(pyqmc/validate_chi.py). EXACT at U=0; U=3 tracks ED across the tau window
(windowed chi_stag ED 0.0510 vs pyqmc 0.0511 PASS). KEY LESSON: the first version
used H=gB^{-1} / g(tau)=BgB^{-1}, correct only for symmetric B,g -> U=0 passed but
interacting C(tau) went negative; the transposes matter. ALWAYS validate dynamic
estimators on an interacting (asymmetric) case. Theory doc updated with the
time-displaced GF section. All four Phase 6 workstreams (WS1a/1b/1c/2/3) now DONE;
remaining: pairing susceptibility chi_d (4-point), altermagnet-geometry pairing,
multi-determinant / BP-RDM adaptive trial. Regression gate bit-for-bit throughout.

## MULTIDET TRIAL — resolved (2026-06-24)
The multi-determinant trial saga concluded. Two issues, now separated:
(1) REAL BUG (fixed): the multidet energy estimator computed the two-body term from
   the weight-averaged Green's function; <n_a n_b> is quadratic in G so a multidet
   trial needs the t_m-weighted average of PER-DETERMINANT Wick energies,
   E_L = sum_m t_m E^(m)/sum_m t_m (Estimators._local_energy). Masked by U=0 / k=1.
(2) Walker-superposition ill-conditioning (real, stands): a coherent superposition
   of sampled CP walkers is a bad multidet basis (W->0; near-duplicates; arbitrary
   phases) regardless of coefficient. So DON'T build multidet trials from walkers.

ENGINE VALIDATED (step a): pyqmc/validate_casci_ed.py builds the exact GS as a
truncated CI in our determinant convention (ci_ground_state, Slater-Condon, E0 =
QuSpin to 1e-14), installs it frozen (TrialWF.set_multidet), sweeps #determinants.
2x2 U=4: CP bias collapses monotonically +1.05 -> 0 and full-N gives E0 with ZERO
variance. So a GOOD multidet trial removes the CP bias; the engine is correct.
NOTE: this is a VALIDATION ONLY (needs ED -> small clusters), not a method.

NEXT (step b, NOT yet done): the SCALABLE multidet trial = self-consistent natural-
orbital CASCI (no ED). Plan: (i) natural orbitals from the QMC 1-RDM (already
computed in adaptive); (ii) bounded active space (fractionally-occupied orbitals
near E_F); (iii) integral transformation U_pqrs = U sum_i W_ip W_iq W_ir W_is
restricted to active orbitals; (iv) build + diagonalize the small CAS Hamiltonian
(Slater-Condon in OUR orbital determinant convention, to avoid QuSpin sign issues
as in step a) -> CI coeffs; (v) self-consistent iteration. This is a mini-CASSCF:
a focused implementation with its own Slater-Condon + validation (full active space
on a small cluster must reproduce ED). Reuse the verified set_multidet + _local_energy.

## REMOTE COMPUTE SERVER (2026-06-24)
A 48-core box is available for the heavy Python (pyqmc/CASSCF) runs:
  ssh amaz@192.168.101.54   (Ubuntu 24.04, 48 cores, 125 GB RAM, 1.7 TB free)
  -- has python3.12 + git; NO conda/compilers/MPI (Fortran code/src stays local).
Setup done: venv ~/qmcenv with numpy 2.5 + scipy 1.18 (pip via PyPI, reachable).
QuSpin has no cp312 wheel -> not installed; ED references run on the laptop (or
validate_casscf falls back to its own full-CI / takes --ed <E0>). Code lives in
~/qmc (rsync pyqmc ed tools docs from the laptop). 
GOTCHAS: (1) run anything long DETACHED on the server (nohup ... &) + poll a
logfile -- a foreground SSH command gets cut off when the local wrapper returns,
killing the remote process mid-run. (2) Python stdout to a file is block-buffered
-> output appears only at process exit. (3) Cross-machine reproducibility: the
altermagnet regression checksums match the laptop BIT-FOR-BIT, but the single-band
hubbard_mixed differs -- LAPACK degenerate-eigenvector phase conventions differ
across OpenBLAS builds (square K is degenerate; the anisotropic altermagnet K is
not). Physics unaffected. Use the server for embarrassingly-parallel scans
(independent runs via background jobs / xargs -P) -- pyqmc is single-threaded per
run, so parallelism = many concurrent runs, not BLAS threads.
