# cb_model — Interaction-driven altermagnetism in the checkerboard Hubbard model

Working folder for the paper. Everything paper-related lives here.

## Paper in one line
Single-band Hubbard model on the anisotropic checkerboard lattice, spin-independent hopping
`t± = t′ ± δ` (t=1, t′=−0.3, U=4). Altermagnetism is **emergent** (from U + geometry), with **dxy**
symmetry, treated by **unbiased CPQMC** (benchmarked by ED). Magnetic **+** pairing paper.

## Folder layout
- `code/`   — figure scripts (regenerate everything)
  - `fig1_model_dnk.py`  — Fig 1: model schematic + mean-field Δn(k) dxy maps (runs anywhere)
  - `fig2_fs_dos.py`     — Fig 2: non-interacting Fermi surface (δ=0, 0.4) + VHS + DOS (runs anywhere)
  - `fig3_phase_diagram.py` — Fig 3: magnetic (n, δ) phase diagram from CPQMC S^z(k). **RUN ON 251** (needs raw data)
  - `fig4_rp_pointwise.py`  — R_p sharpness, point-wise scatter over (n, δ). **RUN ON 251**
  - `fig5_rp_colormap.py`   — R_p sharpness, interpolated colormap over (n, δ). **RUN ON 251**
  - `fig_dos_standalone.py` — DOS panel on its own
- `figures/` — the three paper figures at dpi=600 (names match the LaTeX `\includegraphics`):
  - `1_dnk.png` (Fig 1, from `fig1_model_dnk.py`)
  - `2_vhs_dos.png` (Fig 2, from `fig2_fs_dos.py`)
  - `3_phase_diagram.png` (Fig 3, from `fig3_phase_diagram.py`, run on 251)
  - `dos.png` (standalone DOS component)
- `docs/`    — research plan (docx/pdf), figure captions, reading list, writing skill
- `data/`    — processed phase-diagram CSV

## Data locations (on cluster 251, user phd25imran)
- Full L=14 archive: `~/L14_archive/` (51 folders = 46 doped + 5 half-filling), U=4, t′=−0.3, δ=0–0.4.
- Each folder has: `dir-kVals/sdwz.dat` (S^z(k), used for the phase diagram), `dir-kVals/n_up.dat`,
  `n_dn.dat` (spin-resolved occupations → Δn(k), Δtot), and `dir-kVals/Vertex_*` / `Unpair_*`
  (pairing correlators, k-space, with error bars).
- Pairing channel map (verified from mc2duph.f90): **dx2-y2 = `dwave`**, **dxy = `dd12wave`**
  (df/ddf form factors (+1,+1,−1,−1) on NN / diagonal bonds).

## Status
Done: Fig 1 (model + Δn), Fig 2 (FS + VHS + DOS), Fig 3 (magnetic phase diagram, spiral+diagonal share q/π
colorbar), Fig 4 (R_p sharpness, point-wise), Fig 5 (SDW k-space + real-space, 4-panel wide; δ=0.20, n=0.929/
0.847/0.745 = Néel/spiral/diagonal). Captions for Figs 1–5 in docs/title_abstract_captions.md.
Fig 6 (pairing vertex crossover: d_{x2-y2} -> d_{xy} vs δ, crossover δ≈0.11–0.22, robust at n=0.622 & 0.582).
Captions for Figs 1–6 in docs/title_abstract_captions.md. Research plan written.
Pairing section = 3 figures (equal-time vertex, from L14 data): Fig 6 line cut (all channels + d-wave zoom,
crossover at fixed n, 6_pairing.png), Fig 7 (n,δ) pairing phase diagram (dx2-y2 vs dxy + difference,
7_pairing_pd.png), Fig 8 bare-bubble(Unpair) vs connected-vertex (was "Fig A", 8_d_crossover.png). The k-space
vertex maps figure is CUT from the main set (code fig8_pairing_kspace.py kept for a possible supplement).
Captions Figs 1–8 in docs/. NOTE: these are EQUAL-TIME VERTEX (CP-bias-prone); the dynamic SUSCEPTIBILITY
cross-check (platform + ED gate) is in progress and will decide whether these stay or get replaced.
HEADLINE PAIRING RESULT: among the unconventional channels, checkerboard anisotropy δ selects d_{xy} over
d_{x2-y2} (crossover). NOT global dxy dominance (ext-s larger = conventional/AFM-tied).
POSITIONING (decided w/ user 2026-07-28): LEAD on emergent spin-independent AM (novelty vs paper 1 & 2 which
IMPOSE AM via spin-dependent hopping); pairing crossover = CONFIRMATION, cite paper 2 (unpublished, no scoop).
paper1 = "Enhancement of d-wave Pairing in Strongly Correlated Altermagnet" (imposed dx2-y2 AM, square).
paper2 = user's OWN unfinished "Superconducting pairing symmetry in altermagnetic square Hubbard" (imposed
dx2-y2+dxy via spin-dep tA+t', shows SAME dx2-y2->dxy crossover driven by t'; codes in ~/Desktop/plot_summary.ipynb).
Next: L=16 confirmation; 4-method benchmark; optional pair-field susceptibility to demote ext-s.

## Writing
Apply the `academic-humanizer` skill (docs/SKILL.md): no em-dashes, no semicolons, claims tied to evidence,
channels written as dxy / dx2-y2 (not B2g/B1g).

---

## Code & reproducibility (added 2026-07-30)

`code/` is **self-contained** — the CP-AFQMC engine `cpqmc.py` depends only on numpy, so the
whole susceptibility pipeline runs without the external qmc-platform install.

**Engine + model**
- `cpqmc.py` — constrained-path AFQMC engine (numpy-only).
- `checkerboard.py` — checkerboard model: `checkerboard_hopping` (spin-independent K),
  `nn_bond_factors` / `diag_bond_factors` (s / dx2-y2 / dxy form factors),
  `run_bp_chid_cb` (multi-channel dynamic pairing susceptibility), `run_bp_chi_spin` (chi_zz(q)).

**Drivers (run on the cluster)**
- `checkerboard_chi_grid.py` — (n, delta, U) grid scan of the connected-vertex pairing
  susceptibility. `NSEED` env var sets seeds/point (quick=3, production=6+). Writes per-U CSVs
  + `chi_grid_all.csv`.
  - e.g. `NPROC=32 NSEED=6 python checkerboard_chi_grid.py`
- `checkerboard_chispin.py` — chi_zz(q) staggered magnetic susceptibility.
- `checkerboard_timing.py` — finite-size wall-clock probe (L=6,8,10) for the size-scaling study.

**Plots** — `plot_cb_fig9..15_*.py`, each takes a CSV path as argv and has editable style vars
at the top. Fig 14 = susceptibility phase diagram (per-U CSV); Fig 15 = chi vs U (interaction-driven).

**Result headline.** At U=0 both pairing vertices vanish (~1e-13) despite full anisotropy delta;
they turn on only for U>0 -> the pairing is interaction-driven, not a delta band-structure effect.
dx2-y2 is AFM-mediated (grows toward half-filling); dxy is the altermagnetic channel (half-filling
only, suppressed by doping) -> the channel swap.

**Not included here** (live on qmc48): the exact-diagonalization modules `checkerboard_ed.py` /
`checkerboard_ed_np.py` that `validate_checkerboard_chid*.py` import. The validation scripts are
kept as a record of the ED sign gate (which showed the equal-time vertex is CP-biased while the
tau-integrated susceptibility is sign-faithful).
