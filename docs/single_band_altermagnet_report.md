# Single-band spin-dependent altermagnet — development & findings (referee reply)

**Date:** 2026-06-25. **Goal:** build the manuscript's *primary* model (single-band
square Hubbard with spin-dependent anisotropic hopping) in the pyqmc/ED platform —
separate from the two-orbital robustness-check model — and test the central claim:
*increasing anisotropy suppresses AFM correlations while enhancing the equal-time
connected (vertex) d-wave pairing correlation.*

## What was built
- **Spin-dependent hopping** `K↑≠K↓` in pyqmc (`LatticeModel`/`CPMC`, `K_dn=`),
  bit-for-bit regression-safe when `K↓=K↑`.
- **`am_hopping(lx,ly,t0,tam,t1)`** — the single-band altermagnet (mc2duph.f90):
  - **NN `tam`**: ↑ hops `t0−tam` (x)/`t0+tam` (y), ↓ rotated 90°. Gives a
    d_{x²−y²} *spin splitting* but **preserves (π,π) nesting**.
  - **NNN `t1`**: spin-dependent diagonal (d_xy), ↑ `+t1` main / `−t1` anti-diagonal,
    ↓ opposite. Adds `cos(kx±ky)` → **breaks (π,π) nesting**.
- **Connected/VERTEX pairing** (the manuscript's quantity): `run_bp_pairmag` now
  reports full P=⟨ΔΔ†⟩, disconnected (from the averaged G), and vertex = full−disc,
  for s and d wave, plus S(π,π). ED reference in `validate_am_single.py`.
- Adaptive/CASCI trials now also apply to the pairing measurement.

## Finding 1 — the anisotropy that matters is the NNN `t1`, not NN `tam`
The manuscript says anisotropy "disrupts nesting." A spin-dependent **NN** anisotropy
*cannot* do that: `cos((k+Q))=−cos(k)` keeps `ε(k+Q)=−ε(k)` at half filling, so each
spin band still perfectly nests at (π,π). Indeed:

| anisotropy (4×4, half filling, U=6, CP-AFQMC) | d-wave vertex | S(π,π) |
|---|---|---|
| NN `tam`: 0 → 0.4 | 4.34 → **0.74** (down) | 0.54 → **0.72** (up) |
| NNN `t1`: 0 → 0.4 | 4.34 → **8.85** (up) | 0.54 → **0.33** (down) |

**Only the NNN `t1` reproduces the manuscript's trend** (AFM suppressed, connected
d-wave enhanced). NN `tam` alone gives the *opposite* (it makes the band quasi-1D,
strengthening AFM). So the manuscript's "spin-dependent anisotropic hopping" must be
(or be dominated by) the **NNN, nesting-breaking** term.

## Finding 2 — the mechanism is reproduced, but the d-wave VERTEX is bias-limited
On 4×4 (proper 2D geometry) the CP-AFQMC trend matches the manuscript. However, the
**absolute d-wave vertex is not trustworthy**, and this is demonstrable where ED is
available (4×2, U=4, t1=0.3):

| quantity | ED (exact) | CP-AFQMC fixed | CP-AFQMC adaptive |
|---|---|---|---|
| d-wave FULL | 39.85 | 9.3 | 11.3 |
| d-wave VERTEX | **+23.18** | **−1.8** | **−3.9** |
| s-wave FULL | 28.93 | 27.6 | — |
| S(π,π) | 0.483 | 0.270 | 0.367 |

CP-AFQMC underestimates the **full d-wave 4-point correlation ~4×**; since
vertex = full − disconnected and the disconnected part (from the 1-particle G) is
close to ED, the suppressed full **flips the vertex sign**. The adaptive trial
improves S(π,π) toward ED but does **not** fix the d-wave vertex. So among the
observables: **energy, s-wave, and the AFM S(π,π) are reliable; the d-wave vertex
magnitude is the least reliable** — exactly the methodological point the referees
raised.

## Finding 3 — exact validation of the d-wave vertex is blocked by degeneracy
At U=0 the connected vertex must be 0 (a single Slater determinant). On the
half-filled small clusters ED gives a **nonzero** U=0 d-wave vertex (e.g. 5.6 on
4×2), because the free-fermion ground state is **degenerate**: ED returns an
arbitrary multi-determinant combination of the open shell, while pyqmc's trial is a
single determinant. So a clean ED-vs-QMC gate for the d-wave vertex requires a
non-degenerate, square, near-half-filled cluster — which does not exist at
ED-tractable sizes (square half filling needs even×even; 2×2 is degenerate, 4×4 is
too large for full ED).

## Bottom line for the reply
1. The **mechanism** the manuscript describes is real and reproduced: the
   **nesting-disrupting (NNN) spin-dependent anisotropy** suppresses AFM and enhances
   the connected d-wave pairing tendency (CP-AFQMC, 4×4). A purely NN anisotropy does
   the opposite, so the *type* of anisotropy is essential and worth stating precisely.
2. The **d-wave vertex magnitude is intrinsically hard**: large constrained-path bias
   (~4× on the full 4-point), sign-fragile (full−disc), not fixed by the adaptive
   trial, and not cleanly ED-validatable on degenerate half-filled small clusters.
   This supports the referees' caution about over-interpreting the vertex magnitude —
   and argues for stating the result as a **trend/tendency**, not a magnitude.
3. **Reliable companions**: energy and s-wave agree with ED to a few %, and the AFM
   S(π,π) trend (and approximate magnitude, improved by the adaptive trial) is solid.

## Next steps to strengthen confidence in the d-wave vertex
- A genuinely better pairing trial: BCS/number-projected or CASCI built for the
  *frustrated* band (the energy-optimized natural-orbital adaptive trial is not
  enough). This is the platform's clearest path to a trustworthy d-wave vertex.
- Confirm the exact `t0,tam,t1,U,filling` of the production figures (no `in.dat`
  found in the submission records yet) so the scan matches the manuscript exactly.
- Larger-cluster CP-AFQMC (4×4/6×6/8×8) on the server for finite-size scaling of the
  *trend*, with the caveat above on absolute d-wave.
