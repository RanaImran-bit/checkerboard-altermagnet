# Meeting summary — checkerboard altermagnet (magnetic + pairing)

Slide-ready results summary. Date: 2026-07-29.

## Slide 1 — Model & question
- Single-band Hubbard on the **checkerboard lattice**, spin-*independent* hopping `t_pm = t' +/- delta`
  (t=1, t'=-0.3, U=4). delta = geometric anisotropy.
- Question: does an altermagnet — and its pairing — **emerge** from correlations + geometry alone, without
  imposing it via spin-dependent hopping (as prior work does)?

## Slide 2 — Magnetic result (solid)
- Interactions drive **Neel order near half-filling**; doping -> incommensurate (spiral/diagonal). [Figs 3-5]
- The geometry gives the moment a **dxy spin-splitting**, `m ~ M * delta * sin kx sin ky` -> an
  **emergent altermagnet**. [Fig 1]
- chi_zz(q) peaks at (pi,pi) at half-filling, weakens/shifts on doping. [Fig 11]
- *Novelty: neither prior paper (1 & 2) has emergent AM — both impose it via spin-dependent hopping.*

## Slide 3 — Pairing: the method trap (rigor)
- First used the **equal-time pairing vertex** (as in the PRL): showed a dx2-y2 -> dxy crossover. [Figs 6-8]
- **Benchmarked against exact diagonalization** -> the equal-time d-wave vertex is **constrained-path
  sign-biased** (ED -3.9 vs QMC +0.1). Same issue documented in our PRL appeal.
- The **tau-integrated susceptibility is ED-sign-faithful** -> switched to it as the correct observable.

## Slide 4 — Pairing result: the channel swap (headline)   [6-seed confirmed]
- Sound susceptibility (6x6, ED-validated). Two clean statements:
  - **dxy is SUPPRESSED by the anisotropy** — robust, monotonic, tight bars: chi_dxy_vtx at n=0.78 goes
    -0.39 -> -2.77 as delta 0->0.4 (Fig 9). This is the strongest result.
  - **dx2-y2 stays ATTRACTIVE and is mildly enhanced** at larger delta (+0.47 -> +0.86, growth at delta>=0.3;
    ~2sigma). Present but modest.
- So the enhanced/surviving pairing is the **conventional dx2-y2 channel, orthogonal** to the dxy symmetry of
  the altermagnet. [Figs 9-10]
- *Novelty: imposed-AM models -> pairing matches the AM symmetry; here the **emergent AM decouples pairing
  symmetry from magnetic symmetry**.*

## Slide 4b — WHY dx2-y2 and not dxy (the mechanism)   [Figs 12-13]
- chi_zz(pi,pi) (staggered magnetic susceptibility) RISES toward half filling, peaks ~n=0.89 (Fig 12).
- chi_dx2-y2 turns attractive in the SAME doped window (n~0.67-0.89) -- it TRACKS the (pi,pi) AFM
  susceptibility (Fig 13 overlay). Classic AFM-fluctuation-mediated d-wave (cuprate-like).
- dxy does not couple to (pi,pi) fluctuations and is suppressed. So: the (pi,pi) AFM correlations mediate
  dx2-y2 pairing, while the geometry's dxy channel is killed -> a coherent physical picture, not just a null.

## Slide 5 — Takeaways
1. First **unbiased** demonstration of **emergent, spin-independent** altermagnetism (checkerboard, on-site U +
   geometry).
2. Caught and corrected a **constrained-path artifact** in the pairing vertex via ED — the susceptibility is
   the sound observable.
3. The anisotropy enhances **dx2-y2** pairing while suppressing **dxy** -> magnetic and pairing symmetries
   **decouple**, distinguishing emergent from imposed altermagnetism.

## Suggested figure lineup for the deck
Fig 1 (model + Delta n) -> Fig 3 (magnetic phase diagram) -> Fig 11 (chi_zz(q) maps) ->
Fig 6 (equal-time vertex, "what we first got") -> [ED gate: equal-time is CP-biased] ->
Fig 9 (susceptibility vs delta, "what's correct": dxy suppressed, dx2-y2 enhanced) ->
Fig 10 (susceptibility vs n) -> Fig 13 (magnetism-vs-pairing overlay: dx2-y2 tracks chi_zz(pi,pi)).
Fig 12 (chi_zz(pi,pi) vs n) optional. All in cb_model/figures/{9,10,11,12,13}_*.png.

## Status / caveats to state honestly
- Susceptibility figures are 6x6 (the size where trends resolve); connected-vertex SIGN + TREND are
  ED-validated, absolute MAGNITUDE is constrained-path limited (report sign/trend, not precise magnitude).
- Fig 9 (n=0.778) being upgraded from 2 -> 6 seeds for tighter error bars. 8x8 finite-size confirmation TODO.
- Data/code: qmc-platform CP-AFQMC engine; drivers checkerboard_chi_scan.py / checkerboard_chispin.py; plots
  plot_cb_fig9/10/11.py (all committed, synced to 251).

## Next-week action items (weekly meeting, 2026-07-30)

Points 1-10 from the meeting; points 11-13 added by HoKin over chat (2026-07-30).

1. Use closed-shell (not open-shell) filling for the trial wavefunction.
2. Add susceptibility analysis for the s-wave and extended-s channels.
3. Do the susceptibility analysis for different fillings n.
4. Clarify whether I use the Gamma point or the maximum-q point in the susceptibility.
5. Plot everything for U = 0 as well, and for several U values.
6. Study different U values to show the physics is interaction-driven, not delta-driven.
7. Learn the spectral function: energy gap, DOS vs U, ground/excited states
   (refs: PRB 110, 155120 (2024); the graphene dynamic-response paper).
8. Keep quantum fluctuations in mind for the equal-time vertex pairing.
9. Check whether the magnetic order is stable, and whether the altermagnetism enhances or destroys it.
10. Strengthen the interaction-driven / altermagnetic novelty so it does not read as just another
    spin-dependent (delta-tuned) model.

### Added by HoKin (chat, 2026-07-30)
11. [PRIORITY / headline framing] Quantify how U and delta each drive the altermagnetic order parameter
    m_AM, then use the altermagnetism itself as the x-axis: plot the pairing susceptibility against m_AM
    (instead of against U and delta separately). Directly shows whether pairing tracks the altermagnetism.
12. Finite-size stability of the dx2-y2 susceptibility enhancement: check it is stable against lattice size
    (6x6 -> 8x8 -> 10x10), i.e. a real trend, not a small-cluster artifact.
13. Head-to-head channel comparison under that scaling: benchmark dx2-y2 against s-wave and dxy-wave; does
    dx2-y2 stay the dominant/growing channel as L increases? (Extends point 2.)

### Status (as of 2026-07-30)
- Points 5 & 6 IN PROGRESS: U-scan running on 251 (checkerboard_chi_grid.py, US=[0,2,4,6,8], 6x6).
  U=0 gives chi_d = chi_dxy ~ 3e-13 (numerically ZERO) -> pairing vertices vanish with no interaction;
  U=4 reproduces the channel swap. This is the direct interaction-driven proof for point 6/10.
- Fig 14 = susceptibility phase diagram over (n,delta) at fixed U (per-U CSV). Fig 15 = chi vs U (the
  interaction-driven figure). Both scripts in cb_model/code, driver = checkerboard_chi_grid.py.
