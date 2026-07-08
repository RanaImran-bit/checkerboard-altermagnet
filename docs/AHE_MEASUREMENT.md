# Measuring the anomalous Hall effect (AHE) with the platform's QMC methods

**Companion to:** `docs/PLAN_chi_spin_AHE.md` (work item W4). Suggestions + physics
constraints; ranked by feasibility on the existing engines (ED / DQMC / CP-DQMC / CPQMC).

---

## 0. The symmetry constraint you must respect first

**The current single-band altermagnet model (`am_hopping`: tam, t1, tp all real) has
σ_xy ≡ 0 identically — no QMC needed to see this.**

Two independent reasons:

1. **Reality (effective TRS per spin sector).** Both K_up and K_dn are real symmetric
   matrices, and the Hubbard/HS decomposition keeps every propagator real. A real
   Hamiltonian is invariant under complex conjugation ("spinless time reversal"), which
   forces the Berry curvature Ω_σ(k) → −Ω_σ(−k) with Ω real → σ_xy^σ = 0 for **each spin
   sector separately**, configuration by configuration. There is nothing for interactions
   to rescue.
2. **Altermagnet symmetry.** The altermagnet spin splitting is compensated: K_dn is K_up
   rotated by 90° ([C₄||C₂spin] symmetry). Even with net collinear order, the combined
   symmetry {C₄ᶻ ∘ spin-flip} maps σ_xy^up ↔ σ_xy^dn while the crystal Hall response is
   C₄-odd — so the **charge** Hall conductivities cancel between the spin sectors, while
   **spin** transverse responses survive (this is exactly the "spin-splitter" physics).

**Consequence:** to get a nonzero AHE you must extend the model with a T-breaking,
complex, spin-orbit-type term, e.g.:
- Rashba/altermagnet SOC: λ Σ_k (sin k_y σ_x − sin k_x σ_y) — couples the spin sectors
  (breaks the spin-diagonal structure all engines currently exploit), or
- a spin-diagonal complex NNN hopping i·t₂ σ_z (Haldane-like / "altermagnet + t₂ᶻ"),
  which **keeps the model spin-diagonal** — K_σ becomes complex Hermitian but the
  up/dn-factorized determinant structure of all four engines survives untouched.

The **i·t₂ σ_z route is the recommended minimal extension**: with the altermagnet
exchange it realizes the known altermagnet-AHE mechanism (weak SOC lifts the compensation),
and the only code change is `complex` dtype in the K matrices / propagators. HS fields
stay real; determinants become complex → a **phase problem** appears away from
symmetric points (see per-method notes below).

**What is measurable *without* SOC (and is the natural first target):** the
**spin-splitter / transverse spin conductivity** σ_xy^{z,spin} (a spin current j_y^z in
response to E_x). It is T-even, allowed by the altermagnet symmetry, large (it's the
defining transport signature of altermagnetism), and sign-clean at half-filling.

---

## 1. Ranked measurement routes

### Route A — Transverse **spin** conductivity (spin-splitter), no SOC needed. ★ do first
Kubo: σ_xy^{z} = (1/ω) Im Λ_{j_x^z j_y}(q=0, ω→0), from the imaginary-time correlator

    Λ_xy^{z}(τ) = < j_x^z(τ) j_y(0) >,   j_x = i Σ K_{ij}(σ) (c†_i c_j − h.c.) x-bonds,
    j_x^z = same with a σ^z weight (up − dn).

Implementation is a **direct clone of the new χ_zz(q) machinery** (`spin_susc.py`,
`ftcpmc spin=True`): the current operators are one-body bilinears, so the Wick contraction
is the same disconnected + same-spin-exchange structure — just replace the diagonal vertex
δ_ij S^z with the bond-current vertex matrices J^σ_x, J^σ_y (built from K_σ, spin-resolved,
including the tam/t1 bonds — **the spin-dependent hopping is what makes Λ_{xy} ≠ 0**).
- DQMC: exact at half-filling (sign = 1 at tam = 0; near 1 at small tam).
- CP-DQMC: sign = 1 at half-filling for all tam → the production tool, same argument as
  for χ_zz(q).
- ED gate: Lehmann with the same current operators on 2×2/3×2 — same harness as
  `validate_chi_spin_ft.py`.
- Deliverable: the τ-integrated Kubo weight ∫dτ Λ(τ) (thermodynamic; no analytic
  continuation), plus Λ(τ) curves for a later MaxEnt/Nevanlinna σ(ω) if wanted.

### Route B — Chern marker / topological Hamiltonian from the QMC Green's function (with SOC)
Cheapest *interacting* AHE diagnostic once i·t₂ σ_z is added:
1. **Topological Hamiltonian:** h_top(k) = −[G(k, iω→0)]⁻¹ (Wang–Zhang). The
   time-displaced G(τ) machinery already exists in all engines →
   G(iω_n) = ∫₀^β dτ e^{iω_n τ} G(τ); extrapolate the two lowest Matsubara points to
   ω → 0. Diagonalize h_top(k), compute Berry curvature / Chern number on the (lx, ly)
   grid (FHS lattice-gauge formula — robust even on 8×8 grids).
2. **Local Chern marker** C(r) = −2π Im Tr{P [x̂, P][ŷ, P]} from the equal-time
   one-body density matrix P = <c†c> that `corr_block`/`_accum` already measure —
   zero-frequency-pole approximation, but purely equal-time (cheapest, error bars easy).
Caveat: both identify the *topological* (intrinsic, quantized-part) contribution, not the
full σ_xy(T); they are the right tool for "does tam + SOC drive a Chern transition?".

### Route C — Streda formula via Peierls flux: σ_xy = ∂n/∂B |_μ
Thread one flux quantum 2π/(L_x L_y) (Peierls phases on the bonds; magnetic unit cell =
whole torus) and measure the density shift at fixed μ:
    σ_xy ≈ e² [n(Φ) − n(0)] / Φ.
- Propagators become complex → DQMC weight complex (phase problem); **CPQMC/CP-DQMC with
  the phaseless/cosine projection** are the natural engines (the constrained-path
  infrastructure generalizes; this is the standard AFQMC answer to complex propagators).
- Best signal in a gapped/insulating regime (half-filling with SOC gap): ∂n/∂B is then
  quantized → a *very* forgiving observable (integer × flux), ideal first demonstration.
- ED gate: trivial — same Peierls H on 2×2/4×2, exact ∂n/∂B.

### Route D — Antisymmetric current-current correlator (full Kubo σ_xy)
The honest transport number: σ_xy = lim_{ω→0} (1/ω) Im Λ^A_{xy}(ω) with
Λ^A = ½(Λ_xy − Λ_yx) from <j_x(τ) j_y(0)> (charge currents, model with SOC).
- In imaginary time, the antisymmetric part is odd in ω_n: estimate
  σ_xy ≈ Λ^A(iω_1)/ω_1 (first bosonic Matsubara) — no MaxEnt for a first pass; full
  σ_xy(ω) needs analytic continuation of a *matrix* correlator (hard; do last).
- Sign/phase: only viable where <sign> ≈ 1 (half-filling, or CP-DQMC/phaseless).
- This is also the frequency-resolved extension of Route A (same measurement code path).

### Route E — Cross-checks / free extras
- **Equal-time diagnostics of the AHE precondition:** the momentum-resolved spin
  polarization <n_up(k) − n_dn(k)> (d-wave pattern = the altermagnet splitting) is already
  available from G_σ(k) — verify the splitting survives U before hunting σ_xy.
- **ED σ_xy exactly** on 2×2/3×2 with SOC (full Kubo, Lehmann, all ω) — the arbiter for
  Routes B/C/D, same role ED plays everywhere in this platform.

---

## 2. Suggested execution order

1. **Route A** (spin-splitter σ_xy^z, no model change) — reuses this round's χ_zz(q)
   Wick/bond machinery; ED-gate on 2×2; production CP-DQMC at half-filling vs tam.
   *This is the physically meaningful "Hall-like" measurement for a pure altermagnet.*
2. Add **i·t₂ σ_z SOC** (complex spin-diagonal K) behind a flag; re-gate energy/density
   vs ED (engines stay spin-factorized; only dtype changes).
3. **Route C (Streda)** at half-filling in the SOC gap (quantized, forgiving) with
   phaseless CPQMC; **Route B** (h_top / Chern marker) from the same runs for free.
4. **Route D** last, only where 1–3 make a quantitative σ_xy worth the continuation
   effort.

Per-method summary:

| Route | ED | DQMC | CP-DQMC | CPQMC (T=0) | model change |
|---|---|---|---|---|---|
| A spin-splitter σ_xy^z | gate | half-filling | **production (sign=1)** | BP window | none |
| B h_top / Chern marker | gate | ok | ok | ok | SOC (complex K) |
| C Streda ∂n/∂B | gate | phase problem | phaseless | **phaseless** | SOC + Peierls |
| D full Kubo σ_xy | gate (small) | half-filling only | phaseless | hard (window) | SOC |
