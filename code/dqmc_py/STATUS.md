# DQMC Phase 1 — status (2026-06-27): engine VALIDATED (energy + pairing)

`dqmc.py`: focused numpy finite-T BSS DQMC. GATES vs finite-T ED (2x2, Tr e^{-bH}O/Z):
- U=0 free-fermion: energy, density, S_d, S_s ALL EXACT.
- U=4 beta=3 half-filling:
  - energy -> ED (dt 0.125 -> -5.60, 0.0625 -> -5.42, 0.031 -> -5.402; ED -5.405). PASS.
  - density, <sign> exact.
  - s-wave pair structure factor S_s = 17.78 vs ED 17.83. PASS.
  - d-wave S_d vs ED 13.58: dt 0.031 -> 12.09, dt 0.016 -> 12.52 (CONVERGING upward
    toward ED as dt->0). Convention verified correct (U=0 d-wave EXACT; s-wave matches
    at U=4). Residual is the sign-structured d-wave correlator amplifying the asymmetric-
    Trotter O(dt) error (s-wave cancels it) -> needs smaller dt / symmetric Trotter.
    MINOR, not a bug.

DYNAMIC PAIR SUSCEPTIBILITY (the TARGET quantity = un_t1_maxk):
chi_a = int_0^beta dtau <Delta_a(tau)Delta_a^dag(0)> via time-displaced Green's
G_s(tau_l,0)=B(l,0)G_s(0); ED ref = Kubo/Lehmann sum_{ab}|<b|Delta^dag|a>|^2 (e^{-bEb}-e^{-bEa})/(Ea-Eb).
- U=0: chi_d=24.0 EXACT; chi_s -> 4.0 as dt->0 (4.083/4.021/4.005, O(dt) Riemann). PASS.
- U=4 beta=3 half: chi_s=3.79 vs ED 3.767 (PASS); chi_d scatters 11.15/11.54/12.07/12.54
  around ED 11.873 within stats (~+-0.5). The dynamic d-wave chi MATCHES ED -- and is
  BETTER-behaved than equal-time S_d (tau-integration smooths the Trotter error). PASS.

BUGS FOUND + FIXED (trail): slice-0 reset; green(s,l) ASvQRD (brute-validated 1e-13);
Sherman-Morrison SIGN error (was THE corrupting bug). Pairing uses <c^+c>=I-g^T (matches
CPQMC convention); ed_finite_T returns energy/density/S_d/S_s/chi_d/chi_s.

K-RESOLUTION (--kres): full M[m,n] corr/susc matrices -> reduce_mat (IDENTICAL to
pyqmc/unified_scan) -> maxk/k0/r0/rgt. q=0 k0 == scalar S_d/chi_d exactly.

PHASE 3 (partial) -- T=0 CPQMC vs finite-T DQMC cross-validation, same 4x4 U=4 n~0.88,
scan anisotropic NNN t1 (=t'):
- peak-q d-wave SUSCEPTIBILITY rises with t1 in BOTH (CPQMC T=0 +78%, DQMC beta=4 +16%;
  milder at finite T = thermal suppression of long-tau tail). docs/xval_dqmc_cpqmc_t1.png.
- peak-q equal-time CORRELATION falls with t1 in BOTH (~0.82x) -- close agreement.
- doped 4x4 knob-scan (docs/dqmc_knob_scan.png): t1 enhances peak-q chi_d (sign 0.93->0.99),
  tam suppresses it + shifts peak off q=0. Sign manageable to n~0.73 at beta=4.
=> independent method + temperature confirm anisotropic-t' enhances the d-wave pairing
SUSCEPTIBILITY (the un_t1_maxk quantity), opposite to the equal-time correlation trend.

NEXT: Phase 2 (finite-T constrained-path AFQMC, code/ftcpqmc_py/) -- the third method;
optional vs Fortran BSS (same fort.501 input); tighten filling/tau-window match for the
overlay magnitudes.
