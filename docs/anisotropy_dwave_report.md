# Does altermagnetic anisotropy enhance d-wave pairing? — an exact (ED) study

**Date:** 2026-06-25  **Model:** two-orbital altermagnet (PRL LE20050), single-band
density–density Hubbard interaction U, anisotropic intra-orbital hopping.
**Method:** exact diagonalization (QuSpin) on small clusters — no constrained-path
bias. **Code:** `pyqmc/ed_aniso_scan.py`, estimator `pyqmc/validate_pairmag_am.py`.

## Question
The appeal claims *enhancement of d-wave pairing* under altermagnetic anisotropy,
yet we had found anisotropy *suppresses the long-range d-wave correlation*. Are
these consistent? We scan the anisotropy and resolve the d-wave pair correlation
by distance.

## Definitions
- **Anisotropy** α = 1 − t2/t1 (t1 = intra-orbital strong-direction hopping fixed at
  −1; t2 = weak direction). α = 0 isotropic; α → 1 strongly anisotropic (the
  altermagnetic limit).
- **d-wave pair operator** (intra-orbital, NN): Δ_d† = Σ_orb Σ_δ f_d(δ)
  c†_{m,↑}c†_{m+δ,↓}, f_d = +1 on x bonds, −1 on y bonds. The form factor is
  geometric, so anisotropy enters only through the ground state.
- **Pairing structure factor** S_d = ⟨Δ_d Δ_d†⟩ = Σ_R P_d(R), decomposed by the
  distance R between the two pair centres. S_d = total pairing strength; P_d(R→max)
  = long-range pair coherence (the SC order parameter probe).

## Validation
The two-orbital estimator is exact at U=0 (s-wave 32.000 vs ED 32.000; d-wave ~0).
**Important caveat:** constrained-path QMC *underestimates d-wave pairing badly*
(ED 0.349 vs CPMC 0.068, ~5×, at U=4) because the free-electron node poorly
represents off-diagonal pairing order. The s-wave is only ~2% off. **So the trend
below is taken from exact ED**, not from the (d-wave-biased) QMC.

## Results (U=4, t1=−1)

**Quarter filling (4 electrons = 2↑+2↓ on 8 orbital-sites, i.e. 0.5 e/site; half
filling = 1 e/site = 8 electrons) — anisotropy ENHANCES d-wave:**
| α | S_d (2×2) | S_d (4×2) | P_d(longest), 4×2 |
|------|-----------|-----------|--------------------|
| 0.0  | 0.000     | 0.118     | 0.015              |
| 0.4  | 0.163     | 0.265     | 0.033              |
| 0.8  | **3.392** | **1.732** | **0.216**          |

d-wave grows ~15–∞× with anisotropy, and **every** range grows, including the
longest — enhancement is uniform across length scales.

**Half filling (4/8) — anisotropy SUPPRESSES d-wave, long range most:**
| α | S_d | P_d(R=0) | P_d(R=1) | P_d(R=2, longest) |
|------|---------|----------|----------|--------------------|
| 0.0  | 60.82   | 25.28    | 29.88    | 5.66               |
| 0.4  | 45.49   | 21.56    | 22.01    | 1.92               |
| 0.8  | 39.00   | 19.94    | 18.71    | **0.36** (−94%)    |
Total S_d falls ~36%; the on-site part falls ~21%; the **long-range part collapses
~94%** — exactly the "long-range suppression" found before.

**Intermediate filling (3/8) — abrupt collapse + sign change:**
S_d = 30.4 at α=0 drops to ~0.5 for any α>0, and the long-range P_d(R=2) flips from
+2.83 to **−4.67** (anti-correlated). The isotropic point appears to be a
high-symmetry/near-degenerate state with anomalously strong pairing that any
anisotropy destroys (likely a small-cluster level crossing — interpret with care).

## Conclusion
**The effect of anisotropy on d-wave pairing is filling-dependent, and the two
earlier observations are consistent, not contradictory:**

1. **Doped / low filling (away from half):** altermagnetic anisotropy *enhances*
   d-wave pairing at **all** length scales — this is the regime where the appeal's
   "enhancement of d-wave pairing" holds.
2. **At / near half filling:** anisotropy *suppresses* d-wave pairing, and it
   suppresses the **long-range** coherence far more strongly than the local pair
   amplitude (≈−94% vs ≈−21% at α=0.8) — this is the "long-range suppression" we
   had already seen.
3. **Intermediate filling** shows an abrupt isotropic-point instability (sign-
   changing long-range correlation), a small-cluster level-crossing signature.

So "anisotropy increases d-wave pairing" and "anisotropy suppresses long-range
d-wave correlation" describe **different filling regimes** of the same model: pairing
*strength* is enhanced when doped, while long-range *coherence* is what collapses as
half filling is approached.

## Caveats / next steps
- Exact but **small clusters** (2×2, 4×2): "long range" spans only R≲3; the
  enhancement/suppression *trends* are exact but the thermodynamic limit needs
  larger sizes.
- Larger sizes require QMC, but **CPMC underestimates d-wave ~5×** — to get
  trustworthy large-size d-wave one should use the validated CASCI/adaptive trial
  (which reduces the node bias) or constraint-release, or DQMC; that is the natural
  follow-up to confirm the doped-regime enhancement at scale.
- The isotropic-point degeneracy warrants a symmetry analysis and a t3/t4 (NNN /
  inter-orbital) sweep to map the phase boundary.
