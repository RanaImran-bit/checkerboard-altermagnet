# qmc-platfrom — Student Handoff (START HERE)

**Last updated:** 2026-06-30 · **Repo:** `thkaitools-cell/qmc-platform` · **Branch:** `master`

This is a one-page orientation so a new person can pick the project up. It summarizes
(1) engine status, (2) physics + data locations, (3) figures + manuscript-reproduction
status, and (4) what to do next. For the **full dated history and every quantitative
result**, read `docs/HANDOFF.md` (the deep log) and the per-topic reports in `docs/`.

> **Project in one sentence.** A multi-method, ED-gated QMC platform built to test the
> claim *"the altermagnet enhances d-wave superconductivity in the doped single-band
> Hubbard model"* (the PRL/PRB appeal). Every QMC method is validated against exact
> diagonalization before it is trusted.

**Model.** Single-band **altermagnet Hubbard**, U = 4. Spin-dependent hopping:
`tam` = NN d_{x²-y²} anisotropic hopping (K_up weak-x/strong-y, K_dn rotated 90°);
`t1` = NNN spin-dependent (d_{xy}, anisotropic t′); `tp` = NNN spin-independent (isotropic t′).
Built by `am_hopping(lx,ly,t0,tam,t1,tp)` in both `pyqmc/cpqmc.py` and `code/dqmc_py/dqmc.py`.

**How to run / validate (read `docs/ENVIRONMENT.md` + `docs/VALIDATION.md` first):**
```bash
source tools/env.sh                 # sets up env / gfortran SDKROOT
python pyqmc/regression.py          # bit-for-bit gate (5 fixed-seed configs) — run after ANY pyqmc edit
python pyqmc/validate_chid.py       # ED vs CPMC d-wave susceptibility gate
# Fortran build/run + ED compare:
tools/build.sh src-small && tools/run.sh src-small ...   # then tools/compare.py ed.json qmc.json
```

---

## (1) ENGINE / METHOD DEVELOPMENT STATUS

| Engine | Location | Status | Role |
|---|---|---|---|
| **ED** (ground-state + finite-T, Python/QuSpin) | `ed/hubbard_ed.py`, `ed/altermagnet_ed.py`, `code/dqmc_py/dqmc.py::ed_finite_T`, `pyqmc/validate_*.py` | **DONE / the arbiter** | Exact ground truth that gates every QMC method. 2×2 / 4×4 only (full Fock space). |
| **CPQMC-python** (T=0 constrained-path / BP-AFQMC) | `pyqmc/cpqmc.py` (OOP: LatticeModel/TrialWF/WalkerEnsemble/Propagator/Estimators/CPMC), `pyqmc/unified_scan.py` | **DONE / active R&D surface** | The from-scratch Python port. All estimators ED-validated (energy mixed+BP, equal-time corr, s/d-wave pairing + magnetization, dynamic spin & d-wave pairing susceptibility). Every change gated bit-for-bit by `pyqmc/regression.py`. |
| **CPQMC-fortran (active dev)** | `code/src/*.f90` (cpCore, cpMain, meas, z_*, Vee, Vph …) | **VALIDATED on energy (uxx/uxy/v all 3 channels) vs ED**; one known bug | This is the two-orbital altermagnet (PRL LE20050) tree. Use `out.dat::totalEn` as the production energy. **Open bug:** the separately printed StepMeas *mixed* estimator is still wrong (distinct from the validated `totalEn`). |
| **CPQMC-fortran (frozen reference)** | `code/record/benchmark/CPQMC/cp.f90` etc. | **REFERENCE ONLY (byte-frozen)** | The ORIGINAL single-band manuscript code. Appeal data lives on the cluster (`amax@node-253:~/run/CPQMC/dopingAM/`, READ-ONLY). The Python `unified_scan` reproduces its observables. `modgs` reorthonormalization here is the BP-stabilization reference. |
| **CP-DQMC** (finite-T constrained-path AFQMC, Python) | `code/ftcpqmc_py/ftcpmc.py` | **DONE / validated — the KEY method** | sign = 1 at half-filling for **all** tam (the only sign-clean route into the interesting regime). Force-biased heat-bath, trial-density-matrix constraint, UDV stabilization. Measures everything DQMC does + the τ-integrated pairing eigenvalue. Practical to L ≤ 10. |
| **DQMC-python** (finite-T BSS) | `code/dqmc_py/dqmc.py` | **DONE / validated — workhorse** | ASvQRD-stabilized equal-time + time-displaced Green's; equal-time corr & τ-integrated susc, d-wave + ext-s, full + connected-vertex, k-resolved; `pair_eig.py` pairing eigenvalue. **Sign-free only at tam=0 half-filling** (away from there the sign problem bites — that's why CP-DQMC exists). |
| **DQMC-fortran** (BSS) | `code/record/benchmark/DQMC/BSS.f90` | **REFERENCE ONLY (frozen)** | Source of the ported UDV/ASvQRD stabilization (see `docs/STABILIZATION.md`). Not built this round (gfortran INTEGER-kind shims needed). |
| **Trial-WF R&D** (CP-bias reduction) | `pyqmc/agp*.py` (AGP/number-projected-BCS), `qmc_algo/` (evolutionary/CASCI multidet), `pyqmc/casci.py`, `pyqmc/validate_casscf.py` | **core DONE / ED-validated** | Attacks the one hard blocker — the d-wave-**vertex** CP bias at half-filling. BP-AGP bra recovers ~80% of the ED vertex (7× bias cut); evolutionary selected-CI trial beats adaptive/CASCI at matched determinant count. See `docs/agp_pairing_trial_report.md`, `docs/agp_vertex_plan.md`, `qmc_algo/notes/`. |

**Bottom line:** ED (oracle), DQMC-python, CP-DQMC, and CPQMC-python are all validated and
production-ready at small/medium L. The Fortran trees are validated for energy and used as
references; the active `code/src` tree has one open StepMeas estimator bug.

**NEW (2026-07-08) — momentum-resolved magnetic susceptibility χ_zz(q) in ALL engines**
(plan: `docs/PLAN_chi_spin_AHE.md`; benchmarks: `results/chi_spin_bench/`):
- **DQMC**: already had it (`code/dqmc_py/spin_susc.py::SpinDQMC` + numpy JW ED oracle).
- **CP-DQMC**: `ftcpmc.run_fb_stab(spin=True)` / CLI `--spin` — SpinDQMC Wick formula on
  the CP path (stable G(τ,0)/G(0,τ)/G(l,l), dd-aware) → `chi_spin_q`, `S_spin_q` grids.
- **CPQMC (T=0)**: `CPMC.run_bp_chi_spin()` — full-matrix generalization of the (π,π)
  `chi_block` → C_q(τ) + windowed χ_zz(q) on the whole grid, per-block errors.
- **Gates**: `pyqmc/validate_chi_spin.py` (T=0, numpy-only sector-ED Lehmann oracle —
  no QuSpin needed) and `code/dqmc_py/validate_chi_spin_ft.py` (ED vs DQMC vs
  CP-DQMC-free gated; CP-constrained bias + <sign> reported).
- **Convention (ED-verified):** χ_finite-T(β→∞) = 2 × χ_T=0^one-sided; S^z(q) identical.
- **AHE**: measurement strategy in `docs/AHE_MEASUREMENT.md` — σ_xy ≡ 0 in the current
  real-hopping model (symmetry); do the spin-splitter σ_xy^z first (no SOC needed, clones
  the χ_zz machinery); then i·t₂σ_z SOC + Chern-marker/Streda/Kubo routes.

---

## (2) PHYSICS FOUND + DATA LOCATIONS

**The verdict (robust, multi-method, sign-controlled where it matters):**

1. **d-wave is the dominant pairing channel near half-filling.** Projected vertex
   susceptibility `suscV_d ≫ suscV_s` (DQMC half ~11×, CP-DQMC ~6×, CPQMC closed-shell
   ~1.6–2.4×; grows with L, β). The correct leading-channel test — the **τ-integrated**
   pairing eigenvalue — has a d-wave leading eigenvector at tam=0 half-filling (d-overlap
   0.70 DQMC / 0.80 CP-DQMC). *NB: the equal-time eigenvalue wrongly leads ext-s; always
   use the τ-integrated one.*
2. **tam does NOT enhance d-wave SC — it SUPPRESSES it.** With sign = 1 (CP-DQMC) the
   τ-integrated d-wave character collapses with tam (d-overlap 0.58 → 0.15 by tam=0.2;
   λ → 0 by tam=0.5). The earlier DQMC "tam=0.1–0.2 enhancement" was a **sign artifact**.
3. **What tam DOES do = enhance the equal-time d-wave correlation (amplitude)** while
   suppressing the susceptibility → *"preformed pairs without coherence."* tam acts as a
   momentum-dependent pair-breaker (spin-split FS). Pseudogap-like, not a contradiction.
4. **Not useful for an SC/PDW claim:** the equal-time vertex is **short-range** (decays to
   ~0 by |R|>2, no ODLRO); susceptibility is low/falling (no instability); χ_d(q) does peak
   at finite q under tam but the peak is small and **does not grow with β** → **no PDW**.
5. **Filling dependence:** d-wave leads near half; at over-doping (n≈0.72 closed shell)
   ext-s leads. The SC-dome doping (n≈0.87, L=14/16) run was launched — **status to verify**
   (see "next steps").
6. **t1 (anisotropic t′) suppresses d-wave; tp (isotropic t′) is weak/non-monotonic.** So
   the appeal's "t′-enhancement" framing does not hold; the real story is amplitude-vs-coherence.

> **Subtlety that resolves an apparent contradiction (no bug):** the manuscript plots the
> **local R=0** vertex `N_pp(R=0)` (a k-average), which **rises** with tam; our first scan
> used the **q=0 sum** `P(q=0)=Σ_R N_pp(R)`, which **falls**. Same vertex, two projections.
> `pyqmc/unified_scan.py` now emits all four reductions (maxk / k0 / r0 / rgt) so this is
> explicit. Bond-by-bond verified `am_hopping` == the Fortran `tk`.

**Data locations:**
- `results/dqmc_scan/*.csv` — campaign data: `campA_dqmc_L{6,8,12}` (DQMC backbone),
  `campB_cpdqmc_L{6,8,10,8b10}` (CP-DQMC, sign=1), `campC_cpqmc_L{6,8,10}` +
  `campC_cpqmc_cs_L*` (CPQMC half / closed-shell), `dqmc_vertex2d_*` & `cpqmc_vertex2d_*`
  (tam×t1 vertex maps), `paireig_*` / `petau_cpdqmc_*` (pairing eigenvalues), `cpqmc_bp*` (τ-window scans).
- `results/runs.db` (force-tracked) + `docs/run_registry.md` — **provenance registry**
  (driver + git commit + cluster node + path + observables + status; ~8 campaigns / ~1000
  points). Query: `python pyqmc/run_db.py find --lx 8 --U 4 --tam 0.3 --t1 0`.
- `results/chid_scan/`, `results/un_t1.csv`, `results/un_tp.csv` — susceptibility / reduction scans (some gitignored; numbers are in the reports).
- **Cluster:** mirror `~/qmc/results` on nodes 250–258 (home NOT shared — deploy per node;
  per-node Python differs: 250–252 `~/miniconda3`, 253/258 system, 256/257 `/opt/anaconda3`;
  254/256 have a broken py3.8, skip). Manuscript Fortran data: `amax@node-253:~/run/CPQMC/dopingAM/` (READ-ONLY, different account). Access via `ssh -p <node#> profhokin@thkclusters.duckdns.org`.

---

## (3) FIGURES + MANUSCRIPT REPRODUCTION STATUS

**Key figures (`docs/*.png`) and what they show:**
- `paireig_8x8_doped.png` — pairing eigenvalue (the leading-channel test). **The headline figure.**
- `cpqmc_vertex_map_8x8.png`, `dqmc_vertex_map_6x6_b2.png` — d-wave **vertex** (corr+susc, 4
  reductions) over (t′, tam): tam ↑ peak-q equal-time corr, ↓ q=0 susc.
- `cpqmc_full_map_8x8.png` — FULL (bubble-contaminated) for contrast; `dqmc_4x4_contradiction.png` — shows the "contradiction" IS the bubble.
- `iso_vs_aniso_tprime.png` — t1 vs tp (the subleading t′ effects; t′ does not enhance).
- `xval_three_method_t1.png`, `xval_dqmc_cpqmc_{t1,tam}.png` — cross-method agreement on trends.
- `dqmc_qstar_tam_t1.png` — peak momentum q* (tam pushes the FULL peak off q=0 = bubble).
- `dqmc_chi_beta.png`, `dqmc_btrend_half_8x8.png` — β-convergence (stabilization payoff).
- `cpqmc_landscape_heatmap.png`, `dqmc_landscape_heatmap.png` — (t′, tam) landscapes.
- Manuscript-revision set: `rspace_vertex_vs_tam.png`, `pd_heatmap_Cd0.png` (n×tam phase
  diagram), `tt_heatmap_Cd0.png`, `cs_heatmaps.png`, `un_t1_{maxk,k0,r0,rgt}.png`.

**Manuscript reproduction status:**
- The PRL appeal manuscript reports the d-wave pairing **vertex** `N_pp(R=0)` *increasing*
  with tam (closed-shell n≈0.984, U=4, Fortran CPQMC; data on `amax@node-253`).
- **Reproduced:** `pyqmc` recovers the manuscript's local `C_d(R=0)` increase with tam once
  the observable (R=0, not q=0-sum), U=4, and closed-shell filling are matched
  (6×6 N17: 0.10 → 0.15; U=0 gate = 0).
- **Honest verdict for the reply (see `docs/APPEAL_FINDINGS.md`):** the defensible claim is
  **corrective/negative** — *d-wave is the dominant channel near half-filling, but the
  altermagnet enhances only short-range, incoherent, non-growing local d-wave amplitude and
  suppresses the coherent (susceptibility) channel; net it works against d-wave SC.* This is
  PRB/Comment-level (QMC adding the beyond-mean-field coherence suppression), **not** a
  positive PRL discovery. The remaining quantitative gap is the d-wave-**vertex** CP bias at
  half-filling (being attacked with the AGP / evolutionary trials).

---

## (4) FURTHER DEVELOPMENT SUGGESTIONS

**Physics (make-or-break first):**
1. **Finish the SC-dome doping test** — n≈0.87, L=14/16, CPQMC τ-integrated eigenvalue:
   does d-wave lead in BOTH the projected susceptibility AND the eigenvalue *where SC actually
   lives* (away from half-filling)? Launched on nodes 250/251 — **pull
   `~/qmc/results/cpqmc_nearhalf_L{14,16}.csv` and conclude.** This decides any positive claim.
   → **Analysis script ready:** `python pyqmc/analyze_dome.py results/cpqmc_nearhalf_L14.csv results/cpqmc_nearhalf_L16.csv`
2. **CP-DQMC pairing eigenvalue at the dome doping (sign = 1)** — the cleanest positive-or-
   negative SC statement in the doped regime DQMC cannot reach.
   → **Campaign script ready:** `python code/ftcpqmc_py/cpdqmc_dome_scan.py --lx 8 --ly 8 --betas 4,5,6,8 --nmeas 300`
   → **L=4 demo (2026-06-30):** sign=1 throughout; d-wave leads (d_ov≈0.4–0.8) at tam=0,
     falls to d_ov≈0.2–0.4 at tam=0.2 but still leads ext-s. λ_τ grows with β at tam=0 → coherent d-SC signal. Need L=8/10 for definitive conclusion.
3. **PDW check:** track the χ_d(q) finite-q peak vs β at fixed tam > 0 — it must **grow**,
   not merely be finite-q. So far it does not (PDW route closed at L ≤ 8); need larger L / lower T.
   → **Check script ready:** `python pyqmc/pdw_check.py --lx 8 --ly 8 --tam 0.2 --betas 3,4,5,6,8`
   → **L=4 demo (2026-06-30):** at β=2 peak is at q=0; at β=3 peak shifts to finite q=(π/4,π/4) with ratio=1.43. Finite-size artifact likely — needs L≥8 to confirm.

**Method:**
4. **Close the d-wave-vertex CP bias** — push the BP-AGP bra (`pyqmc/agp_bp_vertex.py`) and the
   evolutionary multidet trial (`qmc_algo/`) to the 4×2 target with a **sparse full-Fock ED**
   reference (dense 4⁸ ≈ 34 GB). A stable variance-controlled AGP/BCS importance function is the goal.
5. **Build the Fortran at manuscript scale** (gfortran INTEGER-kind shims) for an independent
   cross-check of the Python ports; fix the `code/src` StepMeas mixed-estimator bug.
6. **Proper BSE vertex** (Γ = χ₀⁻¹ − χ⁻¹) for the pairing eigenvalue instead of the
   connected-pair-matrix proxy — removes the full-bubble cancellation entirely.
7. **Distributed driver** (use the clusterq broker instead of manual per-node ssh; home is not
   shared, so deploy per node and mind per-node Python paths).
8. **CPQMC T=0 susceptibility** is variance-limited at the default short τ-window — an
   importance-sampled / longer-projection BP, or the AGP open-shell trial, would let CPQMC
   reach the true T=0 susceptibility at half-filling.

---

## Pointers / related work

- **Deep log & every number:** `docs/HANDOFF.md`. Appeal narrative: `docs/APPEAL_FINDINGS.md`.
  Validation source-of-truth: `docs/VALIDATION.md`. Plans: `docs/PLAN_dqmc_cpqmc.md`,
  `docs/PLAN_phase6.md`. Stabilization port: `docs/STABILIZATION.md`. Theory: `docs/theory/`.
- **Per-topic reports:** `docs/single_band_altermagnet_report.md`,
  `docs/pairing_susceptibility_report.md`, `docs/anisotropy_dwave_report.md`,
  `docs/agp_pairing_trial_report.md`, `docs/large_lattice_free_vs_agp.md`, `docs/run_registry.md`.
- **Engine STATUS files:** `code/dqmc_py/STATUS.md`, `code/ftcpqmc_py/STATUS.md`.
- **Web platform** (ED-vs-QMC compare UI): `platform/` — `source tools/env.sh && tools/serve.sh` → http://localhost:5180.
- **Sister repo** `thkaitools-cell/pph-qmc` — a separate Nambu / partial-particle-hole (PPH)
  CPQMC effort (Fortran + Python) targeting the *same* altermagnet d-wave physics via a pairing
  source and PPH-transformed propagation. The AGP anomalous-Wick work here referenced its §6
  Nambu bottleneck. Its development worktree (Fortran `data/code/source/*.f90`, PRB manuscript,
  stripe/science figures) is at `D:\Lee\QMC\NNN-LxLy\worktrees\2026-06-03_alt_request_clean`.
