# Plan — DQMC (finite-T) vs CPQMC (T=0 and finite-T) cross-validation

**Date:** 2026-06-27. **Goal:** independently validate the T=0 CPQMC d-wave pairing
susceptibility (the `docs/un_t1_maxk.png` quantity = peak-k connected d-wave vertex
susceptibility) against **finite-T DQMC** and a **finite-T constrained-path AFQMC**, on
small lattices in the **manageable-sign** regime. This is the method-independent and
temperature-independent check a referee will want for the appeal.

## Decisions (locked 2026-06-27)
- **Item 3 = finite-T CONSTRAINED-PATH AFQMC** (a genuine thermal method with a trial
  density matrix controlling the sign; the finite-T analog of ground-state CPMC).
- **Item 2 = FOCUSED numpy DQMC** — only the observables the comparison needs
  (energy, density, k-resolved d-wave pairing susceptibility), heavily ED-gated; not a
  1:1 port of all of BSS.f90.
- **Start with the STANDARD Hubbard (tam=t1=0), then turn on the altermagnet anisotropy.**

## Reference Fortran (already has the altermagnet model)
`code/record/benchmark/DQMC/`: BSS (Blankenbecler-Scalapino-Sugar) finite-T determinant
QMC. `BSS.f90` (4828 lines, modules mbss+link), `AppBSS.f90` (MPI driver), `BSSOutput.f90`.
Grand-canonical (μ, β=NT·dt), discrete HS, per-slice B-matrices, UDV stabilization,
Metropolis sweeps, sign tracking. Has `tam`/`ttp`/`ttn` (spin-dependent NN hopping),
`th1` (diagonal t'), and measures k-space & r-space pairing + dynamic susceptibilities
(time-displaced Green's `gt1=<c(τ)c†>`, `gt2=<c†(τ)c>`). Input `fort.501`:
`NT / dt / U,mu / warms,runs,sweeps / lamda,tam / h,hx`.

## Phases

### Phase 0 — Document the DQMC Fortran (item 1)
Read BSS.f90 model setup (sysdef/sysinit/cnfinit), the sweep (cnfmake), Green's
(makeg/makeb/udvb stabilization), and measurements (FTSF/FTSus/FTSsupBar/FTGt/FTPC →
PairKSpace/PairRSpace). **Deliverable:** `code/record/benchmark/DQMC/README.md` — model,
algorithm, parameters & input format, observables + output files, build/run, and an
explicit MAP to the CPQMC observables (PairKSpace ↔ unified_scan `P(q)`).
Optional: build the Fortran (gfortran shims for INTEGER kinds) for reference numbers.

### Phase 1 — Focused Python DQMC + benchmark (item 2)
`code/dqmc_py/dqmc.py`: numpy finite-T DQMC for the (altermagnet) Hubbard. Reuse cpqmc
discrete-HS / B-matrices. Stabilized equal-time `G=(1+B_{L}...B_1)^{-1}` (QR/UDV) +
time-displaced G for the dynamic pairing susceptibility; Metropolis HS sweeps; **track
<sign>**. Measure energy, n, and the k-resolved d-wave pairing susceptibility (match
`unified_scan` P(q)/reductions).
**Gates:** (1) U=0 = exact free-fermion finite-T; (2) 2×2/4×4 half-filling vs **finite-T
ED** (`Tr e^{-βH}O / Tr e^{-βH}`, QuSpin) for energy/n/structure-factor/pairing-χ;
(3) vs Fortran BSS on the same input; (4) **β→∞ → CPQMC T=0**.

### Phase 2 — Finite-T constrained-path AFQMC (item 3)
`code/ftcpqmc_py/ftcpmc.py`: finite-T AFQMC with a constrained-path / trial-density-matrix
constraint to tame the sign at finite T (Zhang finite-T CPMC). Gates: U=0 exact; small
ED; reduces to free-projection (= DQMC) when the constraint is released; <sign> ~ 1 in
the controlled regime. Reuse Phase-1 B-matrix machinery.

### Phase 3 — Cross-comparison (the target figure)
Manageable-sign regime: 4×4 (and 6×6), U=4, filling where <sign> stays usable (d-wave SC
is at finite doping → the sign crux; document <sign>), moderate β scanned. Compute the
**peak-k connected d-wave pairing susceptibility** with: T=0 CPQMC (unified_scan, have),
finite-T DQMC (Phase 1) at several β, finite-T CP-AFQMC (Phase 2). **Figure:** maxk-χ_d
vs tam (and vs T/β), overlaying the three, showing DQMC→CPQMC as β→∞. Provenance in run_db.

## Honest scope / risks
- **Sign problem is the crux** — d-wave pairing needs finite doping (repulsive U away from
  half), where finite-T QMC signs decay. Validation lives where <sign> is still usable
  (small lattice, moderate β, modest doping). Always report <sign>.
- **Canonical (CPQMC, fixed N) vs grand-canonical (DQMC, μ):** match filling by tuning μ
  or use canonical projection; compare at matched <n>.
- Finite-T CP-AFQMC (Phase 2) is the hardest and most novel piece; gate it hard.
