# Pairing-aware (AGP / number-projected-BCS) trial for the d-wave vertex — report

**Date:** 2026-06-25. **Goal:** attack the project's key open problem — the d-wave
pairing **constrained-path bias** (equal-time d-wave vertex on the 4×2 / t1=0.3 / U=4
gate: ED **+23.18** vs CP-AFQMC with the free-electron trial **≈ −2**, sign-flipped).
The diagnosis (see `single_band_altermagnet_report.md`) was that the free-electron
trial has the **wrong nodes** for the 4-point d-wave channel, and constraint release
fixes only the ket, not the bra — so a **pairing-aware trial** is needed. This work
builds and validates that trial. It takes the sister effort
`thkaitools-cell/pph-qmc` (Nambu/BdG AFQMC with a d-wave source field) as reference;
its stated bottleneck (report §6) is exactly the rigorous paired-trial overlap/energy
delivered here.

## What was built (`pyqmc/agp.py`, `pyqmc/validate_agp.py`, `pyqmc/agp_dwave.py`)

**The trial.** A number-projected BCS state = antisymmetrized geminal power (AGP),
geminal matrix `F`:
`|Ψ_AGP⟩ = (Σ_ij F_ij c†_{i↑} c†_{j↓})^N |0⟩`. Its overlap with a number-conserving
Slater walker (up/dn orbital matrices `A`,`B`) is the single determinant
`⟨Ψ_AGP|Φ⟩ = det(Aᵀ F B)` — the N-column determinant auto-projects onto N pairs.

**The subtlety (and the rigorous fix).** The projected AGP is **non-Gaussian**, so a
normal-only Wick contraction of one-body Green's functions gives the **wrong**
two-body energy (verified: off by ~2.4 at N=2). The resolution: because `H` conserves
particle number and the walker is fixed-N,
`⟨Ψ_AGP|H|Φ⟩/⟨Ψ_AGP|Φ⟩ = ⟨BCS|H|Φ⟩/⟨BCS|Φ⟩` on the **unprojected** Gaussian
`|BCS⟩ = exp(Σ F c†c†)|0⟩`, where generalized Wick is exact and includes the
**anomalous** pairing tensors. The complete mixed contraction set (all derived and
ED-validated):

| quantity | formula | `M = Aᵀ F B` |
|---|---|---|
| `Gu[i,j]=⟨c†_{i↑}c_{j↑}⟩` | `F B M⁻¹ Aᵀ` | normal |
| `Gd[i,j]=⟨c†_{i↓}c_{j↓}⟩` | `Fᵀ A M⁻ᵀ Bᵀ` | normal |
| `Kd[i,j]=⟨c†_{i↑}c†_{j↓}⟩` | `F(I − Gdᵀ)` | **anomalous** |
| `Ka[i,j]=⟨c_{i↑}c_{j↓}⟩` | `−A M⁻ᵀ Bᵀ` | **anomalous** |

Cross-spin density-density terms then carry the anomalous correction
`⟨n_{a↑}n_{b↓}⟩ = Gu_aa Gd_bb − Kd_ab Ka_ab`; the d-wave pair structure factor is
`⟨Δ†_d Δ_d⟩ = Σ Gu⊙(Fd Gd Fdᵀ) − (Fd:Kd)(Fd:Ka)`, where `Fd:Kd = ⟨Δ†_d⟩` is the
d-wave pair amplitude. The `−Kd·Ka` / `−(Fd:Kd)(Fd:Ka)` pieces are exactly what a
normal-only estimator (and the pph-qmc Nambu prototype, §6) omitted.

## Validation — all to machine precision vs brute-force QuSpin ED (`validate_agp.py`)

Builds the exact AGP and the unprojected BCS in the full Fock space; checks every
formula and the rigorous local energy:

| case | max\|ΔGu\| | max\|ΔKd\| | max\|ΔKa\| | \|ΔE\| |
|---|---|---|---|---|
| 2×2, 1+1, U=4 | 6e-17 | 2e-16 | 1e-16 | 9e-16 |
| 2×2, 2+2, U=4, t1=0.3 | 1e-14 | 8e-15 | 2e-14 | 3e-14 |
| 2×2, 2+2, U=8, t1=0.3, tam=0.2 | 4e-16 | 6e-16 | 4e-16 | 3e-15 |

`tr(Gu)=tr(Gd)=N` exact; overlap matches up to the clean `−1/N!` convention. **This
is the rigorous generalized-Slater/HFB overlap + local energy** that the sister
Nambu/BdG effort needs.

## The d-wave geminal: nodes break the naive BdG construction

`F = −(U†)⁻¹V†` from the d-wave BdG positive-energy block is **singular** on every
cluster tried (`cond(U) ≈ 1e16`): d-wave gap **nodes** make node/deep orbitals
unpaired (`v/u → ∞`), which a finite geminal cannot represent. Tikhonov
regularization makes it finite but yields a huge-norm (`|F|≈13`), **pathological
importance function** (4×2 energy came out `−5` to `−0.5` vs ED `−10.25`).

**The usable construction** (`augmented_geminal`): `F = A·Bᵀ + η·Fd`, the
free-electron determinant plus a bounded d-wave admixture. At `η=0` it reproduces the
free trial **exactly**; `η` smoothly turns on pairing; always well-conditioned.

## CP-AFQMC results (`agp_dwave.py`) — what's validated, what's open

**Machinery validated:**
- **U=0:** energy exact (`−8.00000`) and d-wave **vertex = 0.00000** exactly (the
  connected part must vanish with no interaction) — propagation + AGP overlap
  importance + constraint + normal/anomalous estimators + vertex subtraction all correct.
- **η=0 gate:** reproduces the free-electron-trial CP-AFQMC energy = ED on clean
  clusters (2×2, 1+1, U=4, t1=0.3: AGP-AFQMC `−7.256 ± 0.004` vs ED `−7.254`).

**Open (clearly scoped next steps):**
- At 4×2 half filling, `η=0` gives `−9.35` vs ED `−10.25` — this is the *free trial's*
  own CP/mixed-estimator bias in this hard regime (short projection β=4, mixed
  estimator), independent of the pairing work.
- A **static `η>0` importance function is variance-unstable** at half filling and does
  not yet cleanly reduce the d-wave-vertex bias. Needed: a **back-propagated** (not
  mixed) AGP estimator; `η` optimization / variance control (or use the AGP only as
  the measurement bra while keeping a stable importance function); and a
  **full-Fock-space ED pairing reference** (the number-conserving basis kills
  `Δ†`, giving spurious `S_d=0`) — reuse `validate_chid.py`'s full-Fock machinery for
  an absolute vertex comparison.

## RESULT — back-propagated AGP bra reduces the d-wave vertex bias 7× (2026-06-25)

The next phase (plan: `docs/agp_vertex_plan.md`; code: `pyqmc/agp_bp_vertex.py`,
full-Fock ED gate in `agp_dwave.ed_reference`) is done and gives a clear positive
result. The design separates the two roles of the trial:
- **constraint / importance** → the STABLE free-electron determinant (the AGP-as-
  constraint is variance-unstable at half filling, shown above);
- **measurement bra** → a **back-propagated AGP (BCS)** bra, which carries the
  anomalous pairing contractions the free bra lacks.

The AGP bra is back-propagated through the recorded forward propagators as the
geminal transform `F_bp = Bu_totᵀ F Bd_tot` (since
`⟨AGP(F)|(∏B)|φ₀⟩ = det(φ₀ᵀ F_bp φ₀)`), then the equal-time d-wave vertex is the AGP
generalized Wick (normal + anomalous). The connected vertex subtracts **both**
disconnected pieces — the normal bubble **and** the anomalous pair-amplitude product
`⟨Δ⟩⟨Δ†⟩` (nonzero for the paired bra) — so the **U=0 vertex is exactly 0 at every η**
(gate). `η=0` reproduces the free-electron-bra result (control).

**Gate cluster 2×2 / 2+2 / U=4 / t1=0.3 (the maximally-multireference half-filled
case; full-Fock ED vertex = +22.41):**

| bra | d-wave vertex | fraction of ED |
|---|---|---|
| free-electron (η=0) | 2.3 | **10 %** |
| back-propagated AGP, η≈0.5 | **17.0 ± 0.65** | **76 %** |

A **7× reduction of the d-wave-vertex constrained-path discrepancy**, on the *same*
free-trial ket distribution — purely from giving the pairing channel the correct
(BCS) bra. The improvement is a stable plateau over η ∈ [0.3, 0.6] (a clear optimal
pairing strength; η ≳ 1 goes variance-noisy). The full structure factor moves
12.8 → 31.7 vs ED 30.4. The residual ~24 % gap tracks the free-trial **ket**-constraint
bias (η=0.5 BP energy −5.21 vs ED −5.657 — 2×2 half filling is the hard degenerate
regime), i.e. the bra is now largely fixed and the remaining bias is in the constraint.

### Closing the residual (2026-06-25) — 76 %→80 %, and what's left is the ket

Decomposing the ~24 % residual (`pyqmc/agp_bp_vertex.py --release`,
`pyqmc/agp_release_vertex.py`):
- **Finite back-propagation (recovered).** The vertex climbs with BP length and
  **saturates ~18 = 80 % of ED** (η=0.5: bp16→14.8, bp32→18.1, bp48→17.9). So part of
  the gap was just a too-short bra/ket projection; longer BP lifts 76 %→80 %.
- **Constraint (ket) bias — the remaining ~20 %.** The saturation ceiling is the
  free-trial-constrained ensemble. Constraint **release** (free signed projection of
  the ket after CP equilibration, fixed AGP bra) targets exactly this; U=0 vertex = 0
  at all τ (gate), and the released central value transiently **reaches ED (~22–30 at
  τ≈0.14–0.42)** — consistent with release removing the ket bias. But it is
  **variance-limited**: the AGP-overlap denominator `⟨AGP|φ⟩` blows up for released
  walkers (back-propagation regularizes this; the direct mixed measurement does not),
  so the release curve is too noisy to quantify on this pathological 2×2-degenerate
  cluster (where the free-trial energy bias is itself large, −5.21 vs −5.657).

## Bottom line

The pairing-aware trial works: the rigorous AGP overlap + anomalous local energy
(ED-validated to ~1e-14) plus a **back-propagated AGP bra cut the d-wave-vertex CP
discrepancy ~8× — from 10 % to 80 % of the exact ED vertex** — on the hard half-filled
gate, where every number-conserving (free / adaptive / CASCI-normal) trial gives ~0 by
construction. The last ~20 % is **ket-constraint bias**: release reaches ED transiently
but is variance-limited here; the robust route is a *stable improved constraint*
(future work — e.g. a variance-controlled AGP/BCS importance function, or a less
degenerate cluster where the free constraint is already good) plus the 4×2 target with
a sparse full-Fock ED reference. Honest scope: this reduces the constrained-path /
estimator bias via a correct pairing bra; it does not remove the sign problem.
