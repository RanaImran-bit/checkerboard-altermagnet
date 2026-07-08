# Findings 01 — Algorithm A(a): evolutionary selected-CI trial WORKS

Date: 2026-06-25. Code: `src/evo_trial.py`, gate `tests/validate_evo_trial.py`.
Design: `00_design.md`. Status: **both ED gates PASS on the hardest small case.**

## Result (2x2, nup=ndn=2, U=4, t=1; ED E0 = -5.656854; 36 configs)

This is the *maximally multireference* small case — the HF reference carries only
**0.049** weight in the true ground state; the dominant configs are double
excitations (weight 0.68). It is exactly where a single-determinant CP trial is
worst, so it is the right stress test.

**Gate 1 (exactness):** evolving the archive to the full space reproduces ED to
`|E - ED| = 1.6e-14`. The selected-CI engine + GA are correct.

**Gate 2 (constrained-path bias vs evolved-trial size k):**

| k          | CP bias  | \|<T\|GS>\|^2 |
|------------|----------|---------------|
| free det   | **+1.4945** | —          |
| 1          | +0.1261  | 0.469         |
| 2          | +0.0108  | 0.937         |
| 4          | +0.0097  | 0.970         |
| 8          | −0.0032  | 0.987         |
| 16         | −0.0000  | 1.000         |

The evolved multideterminant trial — built with **NO ED** — drives the CP bias
**monotonically to zero**: a ~140x reduction by k=2, statistically exact by k=8,
machine-exact by k=16. The GA's selected-CI energies match the ED-ranked oracle
truncation at every k, i.e. evolution *recovers the optimal ordering* on its own.

## Three bugs found + fixed (each instructive)

1. **Candidate-generation livelock.** In a nearly-exhausted config space, random
   offspring collide with `seen`; the size-based loop guard never advanced ->
   infinite loop. Fix: bound by attempts, not pool size.
2. **Sign-flipped Epstein-Nesbet selection.** The CIPSI score `coupling^2/(E_A−Hcc)`
   assumes the reference dominates so every external config has `Hcc > E_A`. In a
   multireference state the important doubles have `Hcc < E_A`, flipping the sign;
   a `score<0` filter then **rejects exactly the configs we need** (plateau at the
   bare HF energy −4.0). Fix: select by the *magnitude* of the EN correction
   (regularized denominator), + an immigration step for the exactness gate.
   LESSON: standard single-reference selection criteria fail in the regime this
   method targets — magnitude selection is essential.
3. **Basis-orthogonality overlap collapse (nan).** The trial lives in the
   K-eigenbasis, the *same* basis as the default free-electron walkers, so their
   overlaps are 0/1 and a degenerate multideterminant trial collapses the initial
   overlap -> nan local energy. Fix: re-initialize walkers on the trial's leading
   determinant (standard CPMC practice) so the starting overlap is finite.

## Design choices validated

- **Evolve once to k_max, then truncate by |CI coef|** -> nested truncations ->
  monotone CP bias (independent per-k evolution gave non-nested archives and a
  non-monotone bias).
- **Let the HF reference be dropped** at small k (the multireference state needs
  it gone) — achieved by truncating the full evolution by coefficient.
- **Two-population firewall** (design §3) respected: estimators come from the
  ordinary CPMC walkers; the GA only defines the trial/constraint.

## Honest caveats

- 2x2 is *enumerable*, so "GA matches the oracle" is necessary-not-sufficient: it
  proves correctness, not scaling. The compactness/scaling claim needs a larger
  cluster where full enumeration is infeasible (in progress: 2x4 half filling).
- Single-band square lattice so far. The manuscript target is the d-wave vertex
  at half filling + anisotropy — Algorithm A(b) (walker injection + BCS-seeded
  gene pool) is the next step toward that.

## Head-to-head vs adaptive and CASCI (2x2, U=4; `results/compare_2x2_U4.txt`)

Matched determinant count, same lattice/U/basis/CPMC settings (`tests/compare_trials.py`):

| method      | n_det | CP bias |
|-------------|-------|---------|
| free det    | 1     | +1.5071 |
| adaptive    | 1     | +0.5484 |
| CASCI(0,2)  | 1     | +1.5266 |
| CASCI(1,2)  | 4     | +0.2599 |
| CASCI(0,3)  | 9     | +0.1838 |
| CASCI(0,4)  | 36    |  0.0000 |
| **evolution** | 1   | **+0.1148** |
| **evolution** | 4   | **+0.0166** |
| **evolution** | 9   | **+0.0002** |
| evolution   | 36    |  0.0000 |

**Evolution wins at every matched determinant count:** ~5x lower bias than adaptive
at 1 determinant (0.115 vs 0.548, outside error bars), ~15x vs CASCI at 4, ~900x
vs CASCI at 9 (essentially exact). Mechanism: selected-CI spends each determinant
on the individually most-important config; CASCI fills a whole active window with
many low-weight configs; adaptive is capped at a single determinant. CAVEAT: 2x2 is
enumerable so the GA matches the ED oracle; the matched-k gap is real, but the
SCALING and D-WAVE advantages are A(b)'s job to show.

## Next

1. Larger-cluster scaling gate (GA explores << full space, still beats free det).
2. A(b): inject live CPMC walkers as genes + seed a projected-d-wave reference;
   test on the anisotropic half-filled case that defeated the adaptive trial.
3. Wire `evo_trial` behind a clean trial-builder API so any CPMC run can request
   an evolved trial.
