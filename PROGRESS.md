# Progress tracker — checkerboard altermagnet

Living status of the 13 meeting points (2026-07-30). Legend: ✅ done · 🔄 in progress · ⬜ todo

## Runs in flight
- 🔄 **4-channel 6-seed U-scan** on 251 (`checkerboard_chi_grid.py`, NSEED=6):
  US={0,2,4,6,8} × 8 fillings × 5 δ, channels = on-site s / extended-s / dx2-y2 / dxy.
  Started 2026-07-30 ~16:xx, ETA ~5-6 h. Output: `chi_grid_U*.csv` + `chi_grid_all.csv`.
  U=0 validated: all 4 channels = 0 (~1e-13).

## Meeting points
| # | Point | Status | Notes |
|---|-------|--------|-------|
| 1 | Closed-shell filling for trial WF | ⬜ | need closed-shell filling list per L (differs by size) |
| 2 | s-wave + extended-s channels | 🔄 | driver updated; running now |
| 3 | Susceptibility vs filling n | ✅ | 8 fillings in the (n,δ) grid |
| 4 | Clarify Γ vs max-q | ⬜ | we compute q=0 (Γ); state it explicitly, maybe add max-q |
| 5 | Plot U=0 + several U | ✅ | U-scan 0,2,4,6,8; U=0 blank |
| 6 | Different U → interaction-driven | ✅ | U=0 → χ=0, grows with U (Figs 15/16) |
| 7 | Spectral function (gap, DOS vs U) | ⬜ | learning; refs PRB 110 155120 + graphene paper |
| 8 | Quantum fluctuations in equal-time vertex | 🔄 | ED gate showed equal-time vertex is CP-biased |
| 9 | Magnetic order stability; AM enhance/destroy | ⬜ | have S^z(π,π), χ_zz(q); need stability analysis |
| 10 | Strengthen interaction-driven novelty | 🔄 | U-scan is the core evidence; #11 is the framing |
| 11 | **AM order parameter as x-axis** (HoKin) | ⬜ | PRIORITY: m_AM(U,δ), then χ vs m_AM |
| 12 | Finite-size stability of dx2-y2 (HoKin) | ⬜ | run timing test → size scan 6→8→10 |
| 13 | Compare dx2-y2 vs s-wave vs dxy (HoKin) | 🔄 | s-wave added; compare once run lands |

## Next actions (priority order)
1. When U-scan lands → 4-channel figures → commit driver + result.
2. #12: run `checkerboard_timing.py 6 8` (on node 256, not 251) → size-scan plan.
3. #11: extract m_AM vs U from existing half-filling Fortran folders (U=0,2,3,3.5,4,4.5,5).
4. #1: pick closed-shell fillings per L.
5. #7: read the two spectral-function reference papers.

## Data locations
- Canonical repo: `~/Desktop/checkerboard-altermagnet/` (branch `checkerboard`).
- Cluster runs: 251 `~/qmc-platform-master/pyqmc/`.
- Magnetic U-series (existing): 251 `~/Checkerboard_Model/L14n1.000u*tt0.3N98` (half-filling, δ=0.1-0.3).
