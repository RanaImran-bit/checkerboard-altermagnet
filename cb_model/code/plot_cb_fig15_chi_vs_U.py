#!/usr/bin/env python3
"""Fig 15: interaction dependence of the pairing susceptibility (the interaction-driven test).
Reads chi_grid_all.csv (columns: U,nup,delta,seed,chi_d,chi_dxy,n) from checkerboard_chi_grid.py.
At a fixed delta, plots chi_dx2-y2 and chi_dxy vs U for a few representative fillings; seeds are
averaged (mean +/- sem). Both channels should vanish at U=0 and grow with U if interaction-driven.
  python plot_cb_fig15_chi_vs_U.py ../data/chi_grid_all.csv
"""
import sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- EDIT HERE: selection + style ----
DELTA     = 0.4                 # anisotropy to hold fixed
FILLINGS  = [0.778, 1.0]        # representative n: doped vs half-filling (matched with tolerance)
NTOL      = 0.02
FS_LABEL, FS_TICK, FS_TITLE, FS_LEG, MS = 18, 14, 17, 14, 9
COLORS    = ['#1f77b4', '#d62728', '#2ca02c', '#9467bd']
OUT       = "Fig15_chi_vs_U.png"
# --------------------------------------

CSV = sys.argv[1] if len(sys.argv) > 1 else "../data/chi_grid_all.csv"
d = pd.read_csv(CSV)
d = d[np.isclose(d.delta, DELTA)].copy()

fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.2), facecolor='white', sharex=True)
CH = [('chi_d',   r'$\chi_{d_{x^2-y^2}}$'),
      ('chi_dxy', r'$\chi_{d_{xy}}$')]

for ax, (col, sym) in zip(axes, CH):
    for c, nsel in zip(COLORS, FILLINGS):
        sub = d[np.isclose(d.n, nsel, atol=NTOL)]
        g = (sub.groupby('U')[col].agg(['mean', 'sem']).reset_index().sort_values('U'))
        ax.errorbar(g.U, g['mean'], yerr=g['sem'].fillna(0), marker='o', ms=MS, lw=2,
                    capsize=4, color=c, label=fr'$n={nsel:g}$')
    ax.axhline(0, color='gray', ls='--', lw=0.8)
    ax.set_xlabel(r'Interaction $U$', fontsize=FS_LABEL)
    ax.set_ylabel(sym, fontsize=FS_LABEL)
    ax.set_title(sym + r'  vs $U$', fontsize=FS_TITLE)
    ax.tick_params(labelsize=FS_TICK)
    ax.legend(fontsize=FS_LEG, frameon=False)

fig.suptitle(fr'Interaction-driven test:  connected-vertex $\chi$ vs $U$   '
             fr'($6\times6$, $\delta={DELTA:g}$)', fontsize=FS_TITLE, y=1.02)
plt.tight_layout()
plt.savefig(OUT, dpi=300, bbox_inches='tight', facecolor='white')
print("saved", OUT, " U values:", sorted(d.U.unique()),
      " fillings matched:", FILLINGS)
