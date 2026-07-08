# Altermagnet d-wave pairing — QMC findings for the appeal (2026-06-29)

Three independent methods (DQMC finite-T, CP-DQMC sign-controlled finite-T, CPQMC T=0),
Fortran-referenced stabilization (UDV/ASvQRD, Phases 0-4), cluster campaign over full
tam=0->0.8, L=6,8,10,12, beta up to 10, at half-filling and the closed shell nearest half.

## Result 1 — d-wave is the dominant pairing channel (ROBUST)
- d-wave pair SUSCEPTIBILITY >> extended-s, every clean run: DQMC half b10 ~11x; CP-DQMC
  half ~6x; CPQMC closed-shell 1.6-2.4x. Grows with L and beta.
- The CORRECT leading-channel test (tau-integrated pairing eigenvalue, U=0-gated): at tam=0
  half-filling (sign=1) the leading eigenvector is D-WAVE (DQMC d-overlap 0.70 / CP-DQMC 0.80;
  ext-s ~0; lambda ~0.11-0.19). NB: the EQUAL-TIME eigenvalue instead leads ext-s -- that was
  the wrong observable (see Result 3).

## Result 2 — tam does NOT enhance the d-wave SC susceptibility; it SUPPRESSES it
- suscV_d falls with tam in every clean run; the tau-integrated d-wave eigenvalue falls
  (DQMC b6: lambda 0.38->0.004, d-overlap 0.46->0.01 over tam=0->0.8).
- The only window where it *looked* enhanced (tam=0.1-0.2) is where the DQMC SIGN cratered
  (0.39-0.63) -> being checked with sign=1 CP-DQMC.
- What tam DOES enhance is the EQUAL-TIME d-wave correlation (corrV_d rises, peaks ~tam0.4-0.5).

## Result 3 — equal-time UP while susceptibility DOWN is standard physics (amplitude vs coherence)
- Equal-time <Delta^dag Delta> = int dw B(w) (total pair amplitude, all energies).
- Susceptibility int dtau <Delta(tau)Delta^dag> = int dw B(w)/w (weights LOW energy / long tau).
- tam (the altermagnet up/down Fermi-surface splitting) acts like a PAIR-BREAKING field: it
  hardens the d-wave mode (shifts weight to higher energy / shorter tau-lifetime). So local
  d-wave pairs proliferate (equal-time up) but lose phase coherence (susceptibility down) --
  "preformed pairs without condensation", pseudogap-like. Not a contradiction.

## Does this contradict the manuscript's mean-field?
Not necessarily -- it depends on what the mean field computed:
- If MF gave the pairing AMPLITUDE / gap / DOS-driven bare susceptibility, the QMC
  EQUAL-TIME enhancement is CONSISTENT (correlations do strengthen with tam).
- If MF claimed enhanced SC Tc / full pair susceptibility, the QMC is in TENSION: the
  beyond-mean-field coherence suppression (pair-breaking + fluctuations) that QMC captures
  turns the amplitude enhancement into a susceptibility SUPPRESSION. MF, by construction,
  overestimates the SC tendency by ignoring this.
- Honest framing: QMC confirms the altermagnet strengthens d-wave pairing CORRELATIONS
  (supporting the MF amplitude picture) but the actual condensation tendency (susceptibility)
  is suppressed -- a refinement, not a flat contradiction.

## Appeal recommendation
Lead with: **d-wave is the dominant pairing channel (>> extended-s), robust across method,
size, temperature.** Do NOT claim a tam-driven SC (susceptibility/Tc) enhancement -- present
the tam effect as enhanced pairing AMPLITUDE with suppressed coherence (the equal-time/
susceptibility split), and note it refines rather than confirms a mean-field Tc enhancement.

Code: code/dqmc_py/{dqmc,pair_eig}.py, code/ftcpqmc_py/ftcpmc.py, pyqmc/unified_scan.py.
Data: results/dqmc_scan/{campA,campB,campC,petau,paireig}_*.csv. Cluster: ~/qmc on
profhokin@thkclusters:250-258 (home not shared; per-node python: 250-252 ~/miniconda3,
253/258 system, 256/257 /opt/anaconda3).
