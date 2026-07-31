# Closed-shell fillings — checkerboard Hubbard trial wavefunction

Meeting point #1: use closed-shell (non-degenerate) fillings for the free-electron trial.

**Method.** Diagonalize the non-interacting checkerboard hopping `K = checkerboard_hopping(L,L,-1,0.3,-δ)`.
A filling `nup` is **closed-shell** if there is a gap between the highest occupied level (index nup-1)
and the lowest unoccupied level (index nup): `E[nup] - E[nup-1] > 0`. Then the free-electron ground
state is unique (non-degenerate) and the trial is well defined. This depends on both L and δ.
Filling `n = 2*nup/(L*L)`; `n=1.0` = half-filling. (Lattice is NOT particle-hole symmetric — diagonal
bonds — so `n>1` is a genuine electron-doped side, not a mirror of `n<1`.)

## Closed-shell fillings at ALL δ ∈ {0,0.1,0.2,0.3,0.4}, restricted to n ≤ 1.0

| L | nup | n | Fermi gap (δ=0 … 0.4) | note |
|---|-----|-----|----------------------|------|
| 6  | 9  | 0.500 | 0.10 … 0.43 | dilute-ish |
| 6  | 13 | 0.722 | 0.40 … 0.10 | **good doped point** |
| 8  | 25 | 0.781 | 0.21 … (open at δ=0.4) | use nup=21 (n=0.656) at δ=0.4 |
| 10 | 25 | 0.500 | 0.40 … 0.04 | gap collapses at large δ |
| 10 | 37 | 0.740 | 0.10 … 0.08 | **good doped point** |
| 14 | 49 | 0.500 | 0.07 … 0.07 | marginal (small gap) |
| 14 | 57 | 0.582 | 0.02 … 0.05 | small gap at δ=0 |
| 14 | 61 | 0.622 | 0.22 … 0.08 | **best doped point (largest min gap); = Fig-6 filling** |
| 14 | 73 | 0.745 | 0.03 … 0.08 | good |

## Key facts

1. **Half-filling (n=1.0) is OPEN-shell at every L and every δ** — a degenerate manifold sits at the
   Fermi level (van Hove). The half-filling results (dxy enhancement) therefore rest on an open-shell,
   degenerate trial. Report half-filling as supporting/altermagnetic context, not the quantitative
   headline; lead the clean statements with the doped closed-shell points below.
2. **n ≈ 0.72–0.75 is closed-shell across sizes** (L=6,10,14, all δ) — a natural finite-size series.
3. **Avoid n=0.5 at large L** — its gap collapses (0.01 at L=14) → nearly degenerate.

## Recommended fillings

- **Finite-size study (#12), n ≈ 0.73 series (all closed-shell):**
  `nup = 13 (L=6), 25 (L=8), 37 (L=10), 73 (L=14)`. n drifts 0.72→0.78; at L=8 δ=0.4 is open (drop it
  or use nup=21).
- **L=14 main-text doped points:** `n = 0.582 (nup=57), 0.622 (nup=61)` — closed-shell at all δ;
  n=0.622 matches the earlier equal-time-vertex work (Fig 6).
- Driver: run `FILLINGS=closed python checkerboard_chi_grid.py` to use the closed-shell preset per L.
