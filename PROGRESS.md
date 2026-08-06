# Progress tracker — checkerboard altermagnet

Living status of the 13 meeting points (2026-07-30). Legend: ✅ done · 🔄 in progress · ⬜ todo

## Runs completed
- ✅ **4-channel 6-seed U-scan** on 251 (`checkerboard_chi_grid.py`, NSEED=6), done 2026-07-31:
  US={0,2,4,6,8} × 8 fillings × 5 δ, channels = on-site s / extended-s / dx2-y2 / dxy.
  U=0 = 0 in all channels (~1e-13). Output committed: `data/chi_grid_all.csv`,
  figures `figures/fig_4ch_*.png`, script `code/plot_cb_fig16_channels.py`.

## Meeting points
| # | Point | Status | Notes |
|---|-------|--------|-------|
| 1 | Closed-shell filling for trial WF | 🔄 | fillings computed + `FILLINGS=closed` preset added (`docs/closed_shell_fillings.md`); closed-shell run pending |
| 2 | s-wave + extended-s channels | ✅ | added + validated; on-site s negative (U kills it), ext-s largest positive |
| 3 | Susceptibility vs filling n | ✅ | 8 fillings in the (n,δ) grid |
| 4 | Clarify Γ vs max-q | ⬜ | we compute q=0 (Γ); state it explicitly, maybe add max-q |
| 5 | Plot U=0 + several U | ✅ | U-scan 0,2,4,6,8; U=0 blank |
| 6 | Different U → interaction-driven | ✅ | U=0 → χ=0, grows with U |
| 7 | Spectral function (gap, DOS vs U) | 🔄 | learning; ALF MaxEnt code on 250/251/253/255; refs Jarrell MaxEnt + SAC 2202.09870 |
| 8 | Quantum fluctuations in equal-time vertex | 🔄 | ED gate showed equal-time vertex is CP-biased |
| 9 | Magnetic order stability; AM enhance/destroy | ⬜ | have S^z(π,π), χ_zz(q); need stability analysis |
| 10 | Strengthen interaction-driven novelty | 🔄 | U-scan is the core evidence; #11 is the framing |
| 11 | **AM order parameter as x-axis** | ✅ | M(U,δ) scanned; NO clean collapse (was small-N/U-confound). Honest result = per-filling δ-dependence (fig_pol_chi_vs_delta): doped δ suppresses dxy / dx2-y2 δ-blind; half-filling δ enhances dxy. Half-filling open-shell caveat RESOLVED via anti-periodic BC (fig_apbc_vs_periodic): dxy-vs-δ clean at all U (corr +0.97) → physical, not trial artifact |
| 12 | Finite-size stability of dx2-y2 | 🔄 | driver checkerboard_fss.py ready; timing: ~N^1.5 (L=8 ~15min, L=10 ~30, L=12 ~55). L=6→8 preview: dxy suppression stable/strengthens. Run L=8/10/12 on 251 |
| 12 | Finite-size stability of dx2-y2 | ⬜ | run timing test → size scan 6→8→10 |
| 13 | Compare dx2-y2 vs s-wave vs dxy | 🔄 | comparison done at 6×6 (dx2-y2 leads doped, dxy at half-filling); size-scaling pending (#12) |

## Honest caveats locked in
- Robust δ-effect is **dxy suppression**, NOT dx2-y2 enhancement (dx2-y2 vs δ is non-monotonic, falls at U=8).
- Extended-s is the largest positive channel but conventional/AFM-tied — not the SC instability.
- At low filling both unconventional channels are repulsive; "leading" there = less-repulsive, not real pairing.
- Half-filling is OPEN-shell (degenerate trial) — feature the doped closed-shell points for quantitative claims.

## Next actions (priority order)
1. #12: run `checkerboard_timing.py 6 8` (on node 256, not 251) → size-scan plan.
2. #11: extract m_AM vs U from existing half-filling Fortran folders (U=0,2,3,3.5,4,4.5,5).
3. #7: read ALF MaxEnt code + the two reference papers.
4. #1: run a closed-shell scan (`FILLINGS=closed`).
5. Push repo to GitHub (RanaImran-bit/checkerboard-altermagnet), add thkaitool.

## Data locations
- Canonical repo: `~/Desktop/checkerboard-altermagnet/` (branch `checkerboard`).
- Cluster runs: 251 `~/qmc-platform-master/pyqmc/`.
- Magnetic U-series (existing): 251 `~/Checkerboard_Model/L14n1.000u*tt0.3N98` (half-filling, δ=0.1-0.3).
