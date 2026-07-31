#!/usr/bin/env python3
"""U=0 vs U=8 pairing-susceptibility maps (the interaction-driven contrast).
2x2: rows = U (0, 8), cols = channel (dx2-y2, dxy), shared color scale per channel so
U=0 reads as blank. Demonstrates: no interaction -> no pairing; turn on U -> structure.
  python plot_cb_U0_vs_U8.py chi_grid_all.csv
"""
import sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

CSV = sys.argv[1] if len(sys.argv) > 1 else "chi_grid_all.csv"
NMIN = 0.55
d = pd.read_csv(CSV); d = d[d.n >= NMIN]
g = d.groupby(['U', 'n', 'delta'])[['chi_d', 'chi_dxy']].mean().reset_index()
CH = [('chi_d', r'$\chi_{d_{x^2-y^2}}$'), ('chi_dxy', r'$\chi_{d_{xy}}$')]

fig, axes = plt.subplots(2, 2, figsize=(11, 9), facecolor='white')
for col_i, (col, sym) in enumerate(CH):
    ref = g[np.isclose(g.U, 8.0)][col].values
    vmin, vmax = np.nanpercentile(ref, 5), np.nanpercentile(ref, 95)
    for row_i, U in enumerate([0.0, 8.0]):
        ax = axes[row_i, col_i]; sub = g[np.isclose(g.U, U)]
        sc = ax.scatter(sub.n, sub.delta, c=sub[col], cmap='viridis', vmin=vmin, vmax=vmax,
                        s=320, marker='s', edgecolors='k', linewidths=0.4)
        ax.set_title(fr'{sym}   $U={U:g}$', fontsize=16)
        ax.set_xlabel('Filling $n$', fontsize=13); ax.set_ylabel(r'$\delta$', fontsize=13)
        ax.set_xlim(0.5, 1.05); ax.set_ylim(-0.05, 0.45); ax.tick_params(labelsize=11)
        fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
fig.suptitle(r'No interaction $\to$ no pairing:  $U=0$ vs $U=8$', fontsize=16, y=1.00)
plt.tight_layout()
plt.savefig('fig_U0_vs_U8.png', dpi=200, bbox_inches='tight', facecolor='white')
print("saved fig_U0_vs_U8.png")
