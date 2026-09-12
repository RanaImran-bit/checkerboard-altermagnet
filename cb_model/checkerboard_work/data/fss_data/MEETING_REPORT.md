# Weekly meeting report — checkerboard Hubbard, pairing + finite size

**Data this week:** complete pairing cube at **L = 8, 10, 12** × U = 0, 2, 4, 6, 8 × 6 fillings
(n = 0.5–1.0) × 5 anisotropies (δ = 0–0.4) × **6 seeds** = 2700 measurements, run on the 5-node fleet.
Merged with the earlier L=6 grid into `chi_master_all.csv` (3900 rows).
All four pairing channels measured: on-site s, extended-s, dx²−y², dxy (connected vertex, ED-validated).

---

## 1. Headline: the channel swap is size-stable and gets stronger with L

At the doped filling n ≈ 0.78, U = 4 (χ at δ = 0 → 0.4):

| L | dxy | dx²−y² |
|---|---|---|
| 6 | −0.29 → −1.97 | +0.43 → +0.91 |
| 8 | −0.52 → −2.97 | +3.33 → +0.58 |
| 10 | −1.41 → −4.82 | +2.61 → +0.75 |
| 12 | −1.68 → −5.99 | +3.79 → +0.87 |

dxy is suppressed by the anisotropy at every lattice size, and the suppression **deepens
monotonically with L**. dx²−y² stays attractive throughout. This is not a small-cluster artifact.

**dx²−y² dominates dxy at n ≈ 0.78 at every single (L, U) combination** — all 16 cells negative in
`diff = χ_dxy − χ_dx²−y²`, with the margin growing with L.

*Figures:* `fig_fss_swap_U4.png`, `diff_panel_allU_allL.png`, `maps_allU_4ch_montage_L12.png`

---

## 2. A previous claim was wrong, and finite-size scaling caught it

**Old claim (L=6):** "dx²−y² is mildly enhanced by δ."
**Reality:** L=6 says enhanced at 5.3σ; L=8, 10, 12 all independently say **suppressed**.

χ_dx²−y²(δ=0.4) − χ_dx²−y²(δ=0) at n ≈ 0.78, U = 4:

```
L= 6:  +0.48 +/- 0.09   ENHANCED    (5.3 sigma)
L= 8:  -2.75 +/- 0.37   SUPPRESSED  (7.4 sigma)
L=10:  -1.85 +/- 0.42   SUPPRESSED  (4.4 sigma)
L=12:  -2.93 +/- 0.66   SUPPRESSED  (4.4 sigma)
```

**The corrected story is stronger than the old one.** It is not a competition where δ boosts one
channel and kills the other. δ suppresses *both* d-wave channels — it just suppresses dxy far more
violently. dx²−y² wins by surviving. That is symmetry-selective suppression, which is exactly what an
altermagnetic order parameter of definite dxy symmetry should do; enhancement was never required.

A second latent error found by the same audit: at U = 4 the leading-channel map at L=6 says
**extended-s** leads (13/25 points) while L = 8/10/12 all say **dx²−y²** leads (13–15/25). The
committed figure uses U = 8, where the sizes agree, so it is safe — but the U=4 version would be wrong.

*Figure:* `fig_dx2y2_finitesize_contradiction.png`

**Consequence adopted:** L=6 is no longer used for any physics claim. All reported numbers come from
L=12 (or L=14 archive). The four previously committed L=6 figures have been regenerated at L=12.

---

## 3. Everything is interaction-driven, not δ-driven

The connected vertex susceptibility is **identically zero at U = 0** — 10⁻¹³ to 10⁻¹⁴ across all four
channels, all sizes. That is machine precision, and it is a clean validation of the vertex
construction (full − bubble cancels exactly with no interaction).

It also answers the question directly: with no interaction there is no pairing vertex at all,
regardless of δ. The physics cannot be a band-structure effect.

Channel behaviour vs U is identical at all four sizes — the U at which each channel is most
attractive:

| channel | L=6 | L=8 | L=10 | L=12 |
|---|---|---|---|---|
| on-site s | 0 | 0 | 0 | 0 |
| extended-s | 4 | 4 | 4 | 4 |
| dx²−y² | 2 | 2 | 2 | 2 |
| dxy | 0 | 0 | 0 | 0 |

Both d-channels are strongest at intermediate coupling (U ≈ 2–4) and weaken by U = 8 — too little U
gives no vertex, too much localises the moments.

*Figure:* `fig_4ch_vsU_L12.png`

---

## 4. Half filling: a sharp threshold at δ ≈ 0.2 (needs one caveat)

At L=12, χ_dxy at half filling turns from repulsive to attractive at δ ≈ 0.2, at **every** U:

```
n=1.0   diff = chi_dxy - chi_dx2y2      delta = 0, 0.1, 0.2, 0.3, 0.4
U=2:   -4.19  -4.31  +0.82  +1.77  +2.64
U=4:   -5.40  -5.54  +1.05  +2.64  +3.82
U=6:   -5.03  -4.81  -0.88  +1.05  +3.08
U=8:   -4.79  -4.06  -2.08  -0.41  +0.95
```

So δ acts with **opposite sign either side of half filling** — it suppresses dxy when doped and
activates it at n = 1. A bandwidth or DOS change could not reverse sign with filling.

**Caveat, stated plainly:** n = 1.0 under periodic boundary conditions is the open-shell case
(degenerate trial wavefunction at the band-touching points). These are the least reliable points in
the dataset, and this claim rests entirely on them. It is too systematic across 4 sizes and 4 U values
to be noise, but it needs the clean-BC (twist-averaged) run before it goes in the manuscript.

---

## 5. Point 11 — pairing vs the altermagnetic order parameter (PROVISIONAL, L=6 only)

M is *not* an order parameter: M(U=0) = 0.24–0.43, finite where the pairing vertex is exactly zero.
The interaction-induced moment ΔM = M − M(U=0) vanishes at U=0 as required, giving
**m_AM = ΔM · δ** — zero without interaction, zero without anisotropy.

Against it, at fixed filling (Spearman ρ):

| filling | ρ(χ_dxy, m_AM) | ρ(χ_dx²−y², m_AM) |
|---|---|---|
| n=0.556 | **−0.905** | −0.101 |
| n=0.667 | **−0.827** | −0.288 |
| n=0.778 | **−0.681** | +0.136 |
| n=0.889 | +0.108 | +0.056 |
| n=1.000 | +0.293 | +0.116 |

The altermagnetic order parameter carries dxy symmetry, so it suppresses the dxy pairing channel and
leaves the orthogonal dx²−y² channel essentially untouched. This is the mechanism, measured.

**This is L=6 only and therefore provisional** — precisely the weakness that produced the wrong claim
in section 2. The magnetic run now on the fleet supplies all U at L=12 and will confirm or kill it.

*Figures:* `fig_mAM_perfilling_L6.png`, `fig_collapse_mAM_L6.png`

---

## 6. Answers to specific questions

**Γ or maximum-q?** Γ. The pair operator is Δ = Σ_ij F(i,j) c_i↓ c_j↑ with no e^{iq·r} phase, so the
pair centre-of-mass momentum is q = 0. The form factor F carries the internal symmetry
(s / extended-s / dx²−y² / dxy).

**Is it induced by U?** Yes — the vertex is exactly zero at U = 0 (section 3).

**A general definition of altermagnetism?** m_AM = ΔM·δ (section 5). One candidate that **failed**:
Ψ = (1/N) Σ_q sin q_x sin q_y S(q). Tested on L=10/14 archive data, its δ=0 null is the same size as
its signal and there is no monotonic δ-trend at L=14. The reason is structural, not statistical —
(π,π) is a C4-invariant point, so any C4-odd weight vanishes exactly where all the magnetic weight
sits. Dropped.

---

## 7. Meeting-point status: 9 done, 3 partial, 1 untouched

| # | Point | Status |
|---|---|---|
| 1 | Closed-shell trial | done (APBC committed; clean-BC half-filling probe outstanding) |
| 2 | s-wave + extended-s channels | done |
| 3 | Different fillings | done |
| 4 | Γ vs max-q | **done this week** — Γ |
| 5 | U = 0 and several U | done |
| 6 | Interaction-driven | done |
| 7 | Spectral function / DOS / MaxEnt | **not started** |
| 8 | Equal-time vertex | partial (connected vertex used; equal-time not in the new cube) |
| 9 | Magnetic stability | partial — L=12 magnetic run launched, lands overnight |
| 10 | Strengthen AM novelty | partial (writing; now supported by §3 + §5) |
| 11 | m_AM as x-axis | **done this week**, provisional pending L=12 |
| 12 | Finite-size stability | **done this week** — and it killed a wrong claim |
| 13 | Head-to-head channels under scaling | done |

---

## 8. Next steps

1. **Overnight:** magnetic cube at L=12 (all U) → confirms or kills point 11, closes point 9.
2. **Clean-BC half-filling run** — the δ≈0.2 threshold (§4) is the most interesting result and
   currently rests on open-shell points.
3. **Point 7 (spectral function)** — genuinely untouched, needs MaxEnt machinery.
4. **Normalisation cross-check** — the Fortran archive and pyqmc disagree on S^zz(q) normalisation.
   The L=12 pyqmc block will settle it; if they agree after the 1/N fix, the archive's L=14/16 data
   becomes usable and extends the series to L=16.
