# DQMC (BSS finite-temperature determinant QMC) — code summary

Finite-temperature auxiliary-field determinant QMC (Blankenbecler–Scalapino–Sugar) for
the single-band **altermagnet Hubbard** model — the *same* Hamiltonian as the
ground-state CPQMC (`code/record/benchmark/CPQMC`, `pyqmc/`), but at **finite T**. Used
here as the independent finite-T benchmark for the d-wave pairing susceptibility (see
`docs/PLAN_dqmc_cpqmc.md`). Fortran + MPI; ~5600 lines.

## 1. Model
```
H = - Σ_{<ij>,σ} t^σ_{ij} c†_{iσ} c_{jσ}   (spin-dependent NN hopping)
    - t'  Σ_{<<ij>>,σ} (diagonal)            (th1)
    + U Σ_i n_{i↑} n_{i↓}  -  μ Σ_i n_i
```
- **Altermagnet anisotropy** (`BSS.f90:845-857`): `ttp=1+tam` (strong), `ttn=1-tam`
  (weak). Spin UP hops `-ttn` on the x-bonds and `-ttp` on the y-bonds; spin DOWN is
  rotated 90° (`-ttp` x, `-ttn` y) → the d_{x²−y²} spin splitting. (`K = -t`
  convention, identical to pyqmc `am_hopping`.) `th1` = diagonal t'.
- **Grand canonical**: chemical potential `μ` (vs CPQMC's fixed particle number) — match
  fillings by tuning μ.

## 2. Algorithm (BSS finite-T AFQMC)
- Imaginary time `β = NT·dt` discretised into `NT` slices.
- Trotter: `e^{-βH} = Π_l e^{-dt K} e^{-dt V}` (`emh/eph` = e^{∓dt K}, `emv/epv` =
  e^{∓dt V}).
- **Discrete Hubbard–Stratonovich** of the U term → Ising fields on every (site, slice).
- Per-slice one-body propagators `B_l`; the equal-time Green's function
  `G = (1 + B_{NT}…B_1)^{-1}` (`makeg`/`makeb`), with **UDV/QR stabilization** (`udvb`,
  `udvbt`) for long β.
- **Metropolis sweeps** over the HS fields (`cnfmake`), with **sign tracking** (`negs`,
  `rsign`) — the sign problem appears away from half-filling.
- **Time-displaced** Green's `gt1=<c_i(τ)c†_j(0)>`, `gt2=<c†_i(τ)c_j(0)>` (`FTGt`,
  `makeg(it)`) → dynamic (τ-integrated) susceptibilities.
- Flow (`AppBSS.f90`): `sysdef`→`sysinit`→`cnfinit` (setup) → `sysequil(nwarms)`
  (warmup) → `nmeas` bins of `sysmeas(nswps)` → `measana` → `analysisOutput`.
  MPI = independent Markov chains (one per rank, different seed) averaged.

## 3. Observables (`BSSOutput.f90`, `FT*` routines)
- energy (kinetic `enH`, `enMu`, `enU`), density `<n>`, double occupancy, local moment.
- **structure factors** `FTSF`: CDW, SDW (sdwz, sdwxy).
- **susceptibilities** `FTSus`/`sumSus`: τ-integrated spin/charge.
- **PAIRING** `FTSsupBar` (s- and d-wave, from the time-displaced Green's), `FTPC` (pair
  correlation), `FTDF` (pmdf/dsf). Written **k-resolved → `PairKSpace/`** and
  **r-resolved → `PairRSpace/`** — directly comparable to the CPQMC `Vertex_*`/`P(q)`.

## 4. Input / output
**Input `fort.501`** (unit 500), one value per line:
```
NT                 (# imaginary-time slices)
dt                 (β = NT·dt)
U  mu              (interaction, chemical potential)
warms runs sweeps  (warmup steps, measurement bins, sweeps/bin)
lamda  tam         (pair source field?, altermagnet anisotropy)
h  hx              (pinning fields)
```
e.g. `NT=24, dt=0.125 → β=3`.
**Output**: `fort.200+` (params), `700+` (`<n>` etc.), `800+` (all quantities), `900+`
(run log), `PairKSpace/`, `PairRSpace/`.

## 5. Build / run
`Makefile` (MPI Fortran; note: gfortran rejects some INTEGER(8)/(4) kind mismatches —
shims needed, see handoff). `job.sh` / `paralle.pbs` for the cluster.

## 6. Map to the CPQMC / pyqmc comparison (the target)
| DQMC (this code, finite-T) | CPQMC / pyqmc (T=0) |
|---|---|
| `PairKSpace` d-wave pairing P(q) | `unified_scan.py` `P(q)` (`maxk/k0`) |
| `PairRSpace` N_pp(R) | `unified_scan` `r0`/`rgt` |
| τ-integrated pairing susceptibility | `unified_scan` `suscV` |
| β = NT·dt (temperature knob) | ground state (β→∞ limit) |
| grand-canonical (μ) | canonical (fixed N) |

The cross-check: as **β→∞** the DQMC pairing susceptibility should approach the CPQMC
T=0 value (in the manageable-sign regime); the Python DQMC port (Phase 1) is gated
against finite-T ED and this Fortran, then both are compared to the finite-T CP-AFQMC
(Phase 2). See `docs/PLAN_dqmc_cpqmc.md`.
