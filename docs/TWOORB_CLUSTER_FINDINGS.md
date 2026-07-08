# Two-orbital altermagnet — large-lattice cluster campaign (6x6 & 8x8)

216 runs: 3 methods (CPQMC T=0 full uxx=1/uxy=1/v=0.1 t3=t4=-1; DQMC & CP-DQMC finite-T
U-only reduced, t3=t4=0, half-filling) x 2 sizes x 9 anisotropies (alpha=0..0.8) x 4 seeds.
Run on cluster boxes 251-255. Data: results/dqmc_scan/cluster_twoorb/*, aggregated in
results/dqmc_scan/twoorb_{method}_L{6,8}_agg.csv; figure twoorb_cluster_scan.png.

## HEADLINE: at converged sizes the anisotropy ENHANCEMENT DISAPPEARS
The d-wave enhancement seen at tiny clusters (2x2/3x3/4x4) is a FINITE-SIZE ARTIFACT.
At 6x6 and 8x8:

- EXACT DQMC (sign=1, U-only reduced two-orbital = decoupled anisotropic bands):
  d-wave vertex susceptibility chi_dV is MONOTONICALLY SUPPRESSED by anisotropy:
    6x6: 62.5 -> 15.2   (alpha 0 -> 0.8)
    8x8: 87.9 -> 20.3
  chi_dV grows with L at fixed alpha (real, strongest at the isotropic point; anisotropy
  kills it). This reproduces the ORIGINAL single-band result at large size.

- FULL manuscript model (CPQMC, uxx/uxy/v, t3=t4=-1): the connected d-wave vertex
  susceptibility is small and slightly NEGATIVE (~ -0.3 at 6x6, ~ -0.6 at 8x8), roughly
  flat in alpha -- NO enhancement. Full chi_d and equal-time Cd(0) are flat for alpha<=0.4
  then DECREASE; ext-s dominates the equal-time throughout (Cs(0) ~ 3-4x Cd(0)).

## Caveats
- Finite-T DQMC/CP-DQMC use the U-only reduced model (decoupled orbitals); the genuine
  two-orbital coupling (t4, uxy) lives only in the CPQMC full-model runs.
- CP-DQMC 8x8 returned NaN (walker/stabilization blow-up at 128 sites, beta=5) -- needs a
  smaller dt / better pop control rerun. CP-DQMC 6x6 is noisy but decreasing.
- ALL runs are HALF-FILLING. The manuscript's enhancement DOME is at FINITE DOPING; this
  campaign tests only the anisotropy axis at half filling, where no susceptibility
  enhancement is found.

## Implication for the appeal
Supports the claim-calibrated strategy: do NOT claim a tau-integrated (susceptibility)
d-wave enhancement. At converged sizes and half filling the susceptibility is suppressed
(reduced model) or flat/slightly negative (full model). Consistent with Referee B.
