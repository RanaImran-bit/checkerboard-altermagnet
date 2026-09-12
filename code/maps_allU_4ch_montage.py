#!/usr/bin/env python3
# ALL U in one figure, 4 channels: on-site s / extended-s / dx2-y2 / dxy, rows U=0..8.
import sys, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

CSV = sys.argv[1] if len(sys.argv) > 1 else "chi_grid_all.csv"
OUT = sys.argv[2] if len(sys.argv) > 2 else "maps_allU_4ch_montage.png"
NMIN = 0.55
US   = [0.0, 2.0, 4.0, 6.0, 8.0]
CH = [('chi_son', r'$\chi_{s}$'), ('chi_sext', r'$\chi_{s\text{-}ext}$'),
      ('chi_d', r'$\chi_{d_{x^2-y^2}}$'), ('chi_dxy', r'$\chi_{d_{xy}}$')]

d = pd.read_csv(CSV); d = d[d.n >= NMIN]
cols = [c for c, _ in CH]
g = d.groupby(['U', 'n', 'delta'])[cols].mean().reset_index()
lims = {c: (np.nanpercentile(g[c], 5), np.nanpercentile(g[c], 95)) for c in cols}

fig, axes = plt.subplots(len(US), 4, figsize=(19, 3.4 * len(US)), facecolor='white')
for row, U in enumerate(US):
    sub = g[np.isclose(g.U, U)]
    for col_i, (col, sym) in enumerate(CH):
        a = axes[row, col_i]; vmin, vmax = lims[col]
        sc = a.scatter(sub.n, sub.delta, c=sub[col], cmap='viridis', vmin=vmin, vmax=vmax,
                       s=240, marker='s', edgecolors='k', linewidths=0.4)
        a.set_title(fr'{sym}   $U={U:g}$', fontsize=14)
        a.set_xlim(0.5, 1.05); a.set_ylim(-0.05, 0.45); a.tick_params(labelsize=10)
        if row == len(US) - 1:
            a.set_xlabel('Filling $n$', fontsize=12)
        if col_i == 0:
            a.set_ylabel(r'$\delta$', fontsize=18)
        fig.colorbar(sc, ax=a, fraction=0.046, pad=0.04)
plt.tight_layout()
plt.savefig(OUT, dpi=200, bbox_inches='tight', facecolor='white')
print("saved", OUT)
