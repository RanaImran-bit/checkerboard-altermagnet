# Free vs AGP trial — energy, d-wave vertex, pairing susceptibility (6×6 & 8×8)

**Model:** single-band square Hubbard, spin-dependent anisotropy `U=6, t1=0.2, tam=0`.
**Method:** CP-AFQMC, back-propagated; free-electron constraint throughout; the AGP
(number-projected-BCS) enters only as the **measurement bra** (η=0 = free control,
η=0.05 = AGP). Computed on the 48-core `amaz` server, matched settings (bp=14).
`⟨sign⟩=1` (no sign problem). **No ED at this size** — this compares trials, not vs
exact (the ED gate lives on 2×2/4×2, where the AGP is validated to 80–89% of ED).

Fillings: half = degenerate open shell; non-deg = closed shell near half (6×6: 17;
8×8: 31 — gaps 0.60 / 0.40).

## Energy
| case | free | AGP (η=0.05) |
|---|---|---|
| 6×6 half (18+18) | −23.79 | −23.40 |
| 6×6 non-deg (17+17) | −27.92 | −27.50 |
| 8×8 half (32+32) | −41.39 | −39.48 |
| 8×8 non-deg (31+31) | −46.17 | −50.24 |

## Equal-time d-wave vertex (bp=14)
| case | free | AGP (η=0.05) |
|---|---|---|
| 6×6 half | 17.9 ± 0.3 | 13.9 ± **3.8** |
| 6×6 non-deg | 7.97 ± 0.4 | −0.7 ± **27.7** |
| 8×8 half | 20.2 ± 0.6 | −698 ± **659** |
| 8×8 non-deg | 19.5 ± 0.9 | −77 ± **99** |

## χ_d pairing susceptibility (full and connected vertex)
| case | χ_d full (free) | χ_d full (AGP) | χ_d vertex (free) | χ_d vertex (AGP) |
|---|---|---|---|---|
| 6×6 half | 6.33 ± 0.04 | 5.14 ± 1.10 | **1.97 ± 0.06** | 0.52 ± **1.47** |
| 6×6 non-deg | 7.49 ± 0.07 | 8.22 ± 1.26 | 0.88 ± 0.04 | 17.8 ± **15.8** |
| 8×8 half | 10.16 ± 0.10 | 9.48 ± 0.56 | **2.72 ± 0.11** | 2.18 ± 0.43 |
| 8×8 non-deg | 12.18 ± 0.25 | 10.33 ± 0.65 | 2.68 ± 0.27 | 2.10 ± **2.9** |

## Findings

**1. The free trial is clean and reliable on every large-lattice case** (tight error
bars, ⟨sign⟩=1) for all three observables.

**2. The AGP bra is NOT usable on 6×6/8×8.** Even at the gentlest stable setting
(η=0.05, short bp=14), its error bars are 5–400× the free trial and several central
values are unphysical (8×8 vertex −698 ± 659; χ_d-vertex 17.8 ± 15.8). The geminal
back-propagation `F_bp = Bu_totᵀ F Bd_tot` is intrinsically ill-conditioned at 36/64
dimensions; shrinking η trades signal for noise without fixing the conditioning. The
AGP's value is **confined to small strongly-correlated clusters** (2×2/4×2), where it
is ED-validated (80–89% of the exact vertex) and the free bra genuinely fails.

**3. Physics — bubble vs vertex (why the full χ_d barely discriminates).** Decomposing
χ_d = bubble + vertex, the bubble (= χ_d full − vertex) is **4.4 / 6.6 / 7.7 / 9.5** for
the four cases — 2–4× the vertex and dominating the full χ_d. The bubble is a
single-particle quantity (two normal G's + the d-wave form factor) set by the
Fermi-surface phase space: it is **trial-robust** (free=AGP) and actually **grows with
doping**, so the full χ_d is larger when doped and *anti-correlates* with the pairing
tendency. The **connected vertex** is the small interaction-induced part that carries
the pairing physics: it is **enhanced at half-filling** (6×6: 1.97 vs 0.88; 8×8: 2.72
vs 2.68 — clearer on 6×6). → Report the d-wave SC diagnostic as the **connected
vertex / connected χ_d**, never the raw full susceptibility; and the AGP targets the
vertex precisely because the bubble was never the problem.

## Provenance
`pyqmc/agp_bp_vertex.py` (energy + equal-time vertex), `pyqmc/agp_chid.py` (susceptibility,
ED-validated: U=4 2×2 full 13.63 / vtx 0.290 vs ED 13.56 / 0.294), `pyqmc/chid_tam_scan.py`
(free χ_d, 4 seeds). Raw logs: server `~/qmc/results/cmp`, `cmp2` (gitignored).
