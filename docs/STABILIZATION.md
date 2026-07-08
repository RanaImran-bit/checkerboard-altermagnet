# QMC stabilization schemes (ported from the Fortran BSS/CPMC reference)

Reference Fortran (frozen, read-only): `code/record/benchmark/DQMC/BSS.f90`,
`code/record/benchmark/CPQMC/cp.f90`.

## 1. Equal-time Green's (DQMC) -- ASvQRD / UDV  [already in Python `DQMC.green`]
Fortran: `makeg(it)` -> `makeb` (UDV accumulation) -> `makeipb` (I+UDV factored) ->
`matinv` (stable inverse).
- `makeb`/`udvb`: accumulate the chain `B_{l}...B_1` as `U diag(D) T` (U orthogonal from
  QR, D the signed scales, T unit-upper-triangular), re-orthonormalizing every `north`
  (=`maxrolls`) slices.
- `makeipb`: `I + U D T = U (U^dag T^{-1} + D) T`, re-UDV'd -> `(U Ui) Di (Vi T)`.
- `matinv`: invert `(U Ui) Di (Vi T)` -> `[I/(Vi T)] (I/Di) (U Ui)^dag`; the `1/D` is the
  big/small split.
Python `DQMC.green(s,l0)` does the equivalent: chain QR with `Db=where(|D|>1,D,1)`,
`Ds=where(|D|>1,1,D)`, `inner = U^T/Db + Ds*T`, `G = solve(inner, U^T/Db)`. Validated 1e-13.

## 2. Time-displaced Green's (DQMC)  [Phase 1 -- the fix for suscV]
Physical: `G(tau,0) = <c(tau)c^dag(0)> = B(tau,0) (I + B(beta,0))^{-1}` for tau>0.
Naive `Gt = B_l Gt` accumulates the full dynamic range of `B(tau,0)` -> diverges at large
NT (the suscV=13095 garbage at beta=6).
STABLE identity (derived): with `1+B(beta,0) = 1+B(beta,tau)B(tau,0)`,
```
  G(tau,0) = (B(tau,0)^{-1} + B(beta,tau))^{-1}
```
Let `B(tau,0)=U1 D1 V1`, `B(beta,tau)=U2 D2 V2` (each via the same QR chain). Split
`D1=D1b D1s` (`D1b=max(|D1|,1)`, `D1s=min(|D1|,1)`), likewise `D2`. Then
```
  INNER = D1b^{-1} (U1^dag V2^{-1}) D2b^{-1}  +  D1s (V1 U2) D2s     (well-conditioned)
  G(tau,0) = V2^{-1} D2b^{-1} INNER^{-1} D1s V1
```
Every factor outside INNER is <=1, so no huge product is ever formed. Implemented as
`DQMC._green_tau(s,l)`; the susceptibility re-stabilizes to it every `nstab` slices.

## 3. Walker reorthonormalization (CPQMC / CP-DQMC)  [Phase 3]
Fortran `modgs` (modified Gram-Schmidt): periodically re-orthonormalize the walker
orbitals `Phi` (and carry the log-scale) so the back-propagation imaginary-time chain
stays well-conditioned, allowing a longer BP window (tau) without variance/overflow.

## 3b. CPQMC BP -- already stabilized (Phase 3 finding)
`agp_bp_vertex.py:135` already QR-reorthonormalizes the walkers every BP step
(`pu[i]=qr(pu[i])[0]`) + signed-weight rescaling (line 176) -- the modgs equivalent. So the
BP chain is NOT unstable. The bp-convergence scan (cpqmc_bpconv_8x8.csv) shows suscV does
NOT converge with the tau-window (bp 16->96: 2.5/3.1/4.9/3.1/4.2 -- bouncing), and even the
EQUAL-TIME corrV varies with bp (18.8->29->27) though it must not. => the limit is the
back-propagation ESTIMATOR variance/bias growing with bp, which stabilization cannot fix.
The T=0 connected pair susceptibility via BP is intrinsically noisy; use equal-time corrV.

## Phase 1 payoff: DQMC suscV no longer diverges at large beta
2x2/4x4 beta=6: suscV was 13095 (raw) -> 12.05/7.94 (stabilized), finite & physical.

## Fortran -> Python symbol map
| Fortran | Python |
|---|---|
| `makeg/makeb/udvb` | `DQMC.green` (ASvQRD chain) |
| `makeipb`+`matinv` | the `inner`/`solve` in `DQMC.green` |
| `makeg(it)` (displaced start) | `DQMC.green(s, l0)` |
| time-displaced `gtime`/`FTGt` | `DQMC._green_tau(s, l)` (Phase 1) |
| `modgs` | walker reorthonormalization in CPMC (Phase 3) |
