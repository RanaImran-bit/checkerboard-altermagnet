# v-scan: does the d-wave vertex rise, kink, or fall when S_AM switches on?  (2026-07-02)

Fixed manuscript bands (t1=1,t2=1.75,t3=0.85,t4=0.65) + uxx=4, uxy=1; scan the neighbour
orbital-staggering interaction v=0..1 (the density-density knob that couples directly to the
altermagnetic channel). Data: results/dqmc_scan/vscan/. 4 methods, 2x2..8x8.

## VERDICT: it FALLS. No cooperation between interaction-induced AM and d-wave.
1. ED 2x2 (exact T=0, benchmark): S_AM onsets at v_c<0.2 (degenerate-manifold artifact of the
   tiny cluster) and S_d COLLAPSES ~200x (10.17 -> 0.06). Where AM order exists, d-wave dies.
2. DQMC 4x4 (exact finite-T, sign>=0.54 for v<=0.6): S_AM flat, but d-wave falls monotonically
   anyway: SdV 1.22 -> -1.01, chi_dV +0.42 -> -2.19. v is d-wave-pair-breaking BEFORE any AM
   order forms. ext-s weakens but stays leading.
3. CPQMC 6x6 (T=0 ground state, full v range): NO S_AM onset up to v=1 (flat ~0.085 then
   decreasing). pair_dwave_vertex falls 3.7 -> negative (~ -5) by v>=0.6; ext-s falls too.
4. DQMC 6x6/8x8: v catastrophically worsens the sign problem (<sign> 0.89 -> 0.18 at v=0.2,
   dead beyond) -- only v=0/0.2 usable; consistent with the small-lattice trends. CP-DQMC NaN
   for v>0 (stabilization); v=0 point cross-checks CPQMC.

## Physics conclusion
In the two-orbital altermagnet Hubbard model, the interaction that would induce the
altermagnetic order suppresses d-wave pairing at every accessible point:
- before the onset (4x4 exact): monotone fall, vertex turns repulsive;
- at the onset (2x2 exact): collapse;
- in the ground state at 6x6: no spontaneous AM up to v=1, d-wave vertex negative anyway.
Combined with the U-scan (S_AM flat to U=8; ext-s dominant; imposed-anisotropy scan suppressing
chi_d), the QMC picture is coherent: the mean-field scenario "interaction-induced altermagnetism
enhances d-wave" does NOT survive fluctuations in the microscopic model. The t-J mean-field
regime (J~1-5) is unreachable from Hubbard (J=4t^2/U<=1 for U>=4).

## Use for the appeal
- Do NOT lean on the two-orbital mechanism for d-wave enhancement.
- USABLE: the demonstrated sign-problem quantification (sign 0.89->0.18 at v=0.2, 6x6) backs the
  manuscript's "tau-integrated susceptibility computationally inaccessible" statement concretely.
- The strongest affirmative evidence remains the single-band doped-RVB dome (main-text Fig.5).

## Caveats
- CPQMC constraint bias could suppress an incipient order (mitigated: exact 4x4 agrees).
- Finite-T DQMC densities sit at n~0.87 for mu=U/2 (non-bipartite bands), consistent across v.
- The true leeb2024 driver is Heisenberg S.S exchange (not density-density); v is the closest
  density-density proxy. A genuine J-term QMC would need a new (likely sign-problematic) HS.
