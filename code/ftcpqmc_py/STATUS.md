# FT-CPMC Phase 2 — status (2026-06-27): FORCE-BIAS + CONSTRAINT VALIDATED

## UPDATE: force-biased importance sampling WORKS (--fb); constraint validated at doping
`run_fb` / `_slice_forcebias`: per-site heat-bath sampling of the slice HS fields using
the importance Green's Phi^s=(I+A^s C^s)^{-1} (reference d=1, matching I_{l-1}), both field
values summed (Sherman-Morrison), walker weight = product of per-site normalizations,
constraint zeros the negative-weight branch. BUG FOUND+FIXED: reference must be d=1 (no-HS,
= I_{l-1}), NOT all-spin-up -- the wrong reference dropped a det factor and collapsed S_d
(monotonically wrong vs beta). After fix:
- U=0 exact.
- U=4 half, force-bias FREE: S_d TRACKS ED across beta (b1 9.07/8.92, b2 11.81/11.15,
  b3 12.75/13.58); Trotter-convergent (b2: dt 0.125->0.0625 gives E -5.38->-5.25 vs ED
  -5.223, S_d 11.81->11.08 vs 11.15). Force bias killed the free-projection variance.
- U=4 DOPED mu=1.0 beta=2 + CONSTRAINED: <sign>=1, density 0.854 vs ED 0.854 (exact),
  energy -5.99 vs -5.97, S_d 5.88 vs ED 5.94 (~1%), S_s 17.55 vs 17.38. THE GOAL: a
  sign-controlled finite-T method matching ED at doping where bare DQMC sign decays.
## DONE since: tau-integrated chi_d + partial stabilization + 3-method overlay
- chi_d (--chi): time-displaced G_s(tau,0)=B(l,0)G_s(0) along the path (= DQMC construction).
  Gated vs ED Kubo (2x2): U=0 chi_d=24.0 exact; U=4 half b1 6.61/6.55, b2 10.12/10.17,
  b3 10.84/11.87. FT-CPMC now measures the SUSCEPTIBILITY too.
- THREE-METHOD overlay docs/xval_three_method_t1.png: FT-CPMC b2 chi_d AND equal-time C_d
  both overlap the exact DQMC b2 across the whole t1 scan; CPQMC T=0 same trends steeper.
  (susceptibility UP, correlation DOWN with t'.) 4x4 b2 FT-CPMC vs DQMC: S_d 18.46 vs 18.26.
- STABILIZATION (--stab): measurement Green's via DQMC.green (ASvQRD, brute-validated).
  EXACT no-op at b2 (matches unstab to all digits); at b6 fixes ENERGY (-5.459 vs ED -5.553)
  where the raw inverse gives garbage (-0.555). PARTIAL: only the measurement Green's is
  stabilized -- the IMPORTANCE side (Tpow=(b0)^{L-l} powers + Phi=(I+AC)^{-1} in the sampler)
  still overflows at large beta, so S_d at b6 stays corrupted (-3.1). The cross-validation
  regime beta<=4 is fully validated; beta>=6 needs importance-side UDV stabilization (the
  remaining robustness item). The constraint's CP bias also appears at b4 doping (sign=1 vs
  DQMC 0.84) -- expected.

## P0.1 DONE: importance-side stabilization (run_fb_stab, --fb --stab)
Reformulated the sampler: instead of the boundary importance with the overflowing
(b0)^{L-l}, maintain the STABLE equal-time Green's at the current slice (past = sampled B,
future = trial b0) via DQMC.green (ASvQRD), and do the per-site heat-bath with the LOCAL
DQMC ratio R_s = 1 + (1-G_l[i,i]) dv (O(1), stable) + SM update. Reference field x=0
(HS factor 1 = b0) matches I_{l-1}.
GATES: U=0 exact; 2x2 U4 beta=2 half matches unstab + ED (E -5.250, S_d 11.08 vs ED 11.15);
beta=6 half PHYSICAL now (E -5.67 vs ED -5.55; unstab gave garbage E=-0.55, S_d=-3.1);
4x4 U4 beta=2 doped matches EXACT DQMC to ~1% (S_d 18.43 vs 18.26, n 0.883 vs 0.884, sign=1).
The free-projection S_d at large beta still has variance (use CONSTRAINED in production; DQMC
anchors). COST: O(L) green() calls per slice -> O(L^2) per path, slow at large L (beta=6 L=96)
-> production needs the cluster + parallel paths (and/or a G_l->G_{l+1} wrap optimization).

NEXT: P1 6x6 overlap validation (3 methods); P2 8x8 beta-series production.

--- (prior) FOUNDATION note ---

`ftcpmc.py`: finite-T constrained-path AFQMC. Walkers build the path propagator
M_l^s = B_l..B_1; importance function I_l = prod_s det(I + T_l^s M_l^s) with trial
propagator T_l^s = (b0_s)^{L-l} (b0 = U=0 B-matrix). Telescoping prod_l r_l = W(x)/Z_T
(EXACT, verified by construction); constraint drops walkers when an incremental ratio
r_l <= 0 (Zhang finite-T CP), released (sign-carried) = free projection.

GATES vs finite-T ED (dqmc.ed_finite_T), 2x2:
- U=0 free: energy/density/S_d/S_s ALL EXACT.
- U=4 half, small beta: beta=0.5 S_d 8.18 vs ED 8.15, beta=1 S_d 8.91 vs ED 8.92 (PASS);
  energy within ~3% (free-projection variance).
- U=4 DOPED mu=1.0 beta=1: density 0.865 vs ED 0.863 (PASS), energy -5.56 vs -5.43.
- Constraint mechanism in place; a no-op where <sign>=1 (small beta / this regime) -- it
  activates only when det ratios go negative (larger beta + doping).

KNOWN LIMITATION: uniform field sampling -> free-projection VARIANCE grows with beta
(beta=3 half energy biased to -5.81). The scaffolding (propagation, importance fn,
boundary Green's G=(I+M_L)^{-1}, observables) is CORRECT (U=0 exact, small-beta U=4 exact);
precision at beta=4 (the cross-validation point) needs:
  NEXT: (1) force-bias importance sampling of the HS fields (per-site, via the running
  Green's, like CPMC/DQMC sweep) to tame variance; (2) scale-preserving QR stabilization
  of M for large beta/lattice (currently raw product, ok only at small beta); (3) then the
  constraint becomes meaningful at doping -> add the sign-controlled point to the
  DQMC/CPQMC cross-validation overlay (docs/xval_dqmc_cpqmc_t1.png).
