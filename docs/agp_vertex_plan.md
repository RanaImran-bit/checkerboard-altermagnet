# Plan — reduce the d-wave VERTEX discrepancy with the BCS (AGP) trial

**Date:** 2026-06-25. **Target:** the equal-time d-wave pairing vertex CP bias on the
4×2 / t1=0.3 / U=4 gate (ED **+23.18**, free-electron-trial CP-AFQMC **≈ −2**,
sign-flipped). This is *my* BCS-trial thread (`pyqmc/agp.py`, `agp_dwave.py`,
`docs/agp_pairing_trial_report.md`); the `qmc_algo/` evolutionary effort is separate
and left untouched.

## Diagnosis of the bias (where it enters, and which the BCS trial fixes)
The mixed/BP pairing estimator is `⟨ψ_T| Δ†_d Δ_d |Φ⟩/⟨ψ_T|Φ⟩`. Two sources of bias:
1. **the ket** — the constrained-path walker distribution (set by the constraint /
   importance function);
2. **the bra** — back-propagation pushes `ψ_T` toward the constrained GS, but with a
   *free-electron* `ψ_T` the 4-point pairing channel is mis-nodal even after BP
   (HANDOFF: BP free-bra vertex is still sign-flipped).

The BCS/AGP trial attacks **(2)**: its bra carries the **anomalous** contractions
`Kd=⟨c†↑c†↓⟩`, `Ka=⟨c↑c↓⟩` that directly populate the pairing channel — the d-wave
structure factor with an AGP bra is `Σ Gu⊙(Fd Gd Fdᵀ) − (Fd:Kd)(Fd:Ka)`, and the
`−(Fd:Kd)(Fd:Ka)` term is absent for any number-conserving (free/adaptive/CASCI-normal)
trial. That term is exactly the piece the equal-time free-bra vertex was missing.

## Phases (each ED-gated)

**Phase 1 — full-Fock ED equal-time vertex reference.** The number-conserving basis
kills `Δ†` (spurious `S_d=0`); reuse `validate_chid.py`'s full-Fock Lehmann machinery
at τ=0 to get the exact `⟨Δ_d Δ_d†⟩` full + connected vertex on 2×2 (cheap gate) and
4×2 (the target). Gate quantity for everything below.

**Phase 2 — back-propagated AGP-bra estimator.** Back-propagate the AGP bra: under the
recorded forward one-body propagators the geminal transforms `F → Bu·F·Bdᵀ` (the
AGP analogue of `Estimators.bp_bra`). Measure the equal-time d-wave vertex from the
back-propagated AGP bra via the normal+anomalous generalized Wick. Gates: U=0 →
vertex = 0 exactly; `energy_bp` matches ED on 2×2.

**Phase 3 — decouple constraint from bra.** The AGP-as-importance-function is
variance-unstable at half filling (report). So keep the **stable free-electron
determinant as the CP constraint/importance** (validated, well-behaved ket) and use
the **BP-AGP bra only for the measurement**. This isolates the "better bra" benefit
from the importance-function instability and is the cleanest decisive test: does a
pairing bra alone move the vertex toward ED on the free-trial ket distribution?

**Phase 4 — scan η, compare to ED, record.** On 4×2 half-filled U=4 t1=0.3 scan the
augmented-geminal pairing strength `η` (F = ABᵀ + η Fd); measure the BP-AGP-bra
d-wave vertex vs the Phase-1 ED reference. Record the bias vs η. Success = the vertex
sign/ magnitude move toward ED where the free-bra fails.

**Later (if Phase 3 helps but is bra-limited):** promote the AGP to the *constraint*
too, with variance control (η ramp, partial-projection / hybrid importance, or
release from the AGP-constrained ensemble). Also: self-consistent η (optimize the
trial's variational vertex), and the BdG frozen-core geminal for a node-faithful
trial beyond the augmented form.

## Honest framing
This reduces the **constrained-path / estimator bias** by giving the pairing channel
a correct bra; it does **not** remove the sign problem. Report the d-wave vertex as a
bias-reduced estimate with the ED gate, not an exact number.
