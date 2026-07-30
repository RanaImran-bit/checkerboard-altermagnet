#!/usr/bin/env python3
"""Fig 14: connected-vertex pairing susceptibility phase diagram over (n, delta).
Reads chi_grid.csv (columns: nup,delta,seed,chi_d,chi_dxy,n), averages seeds per (n,delta).
Panels: (a) chi_dx2-y2, (b) chi_dxy, (c) difference (red = dxy favored, blue = dx2-y2 favored).
  python plot_cb_fig14_chi_phasediagram.py ../data/chi_grid.csv
"""
import sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- EDIT HERE: style ----
NMIN = 0.55         # drop the dilute limit (noise-dominated, off-story, hijacks color scale)
FS_LABEL = 20
FS_TITLE = 19
FS_TICK  = 16
FS_CBAR  = 14
S        = 340      # marker size
OUT      = "Fig14_chi_phasediagram.png"
# --------------------------

CSV = sys.argv[1] if len(sys.argv) > 1 else "../data/chi_grid.csv"

d = pd.read_csv(CSV)                          # nup,delta,seed,chi_d,chi_dxy,n
d = d[d.n >= NMIN].copy()
g = (d.groupby(['n', 'delta'])
       .agg(chi_d=('chi_d', 'mean'), chi_dxy=('chi_dxy', 'mean'))
       .reset_index())
g['diff'] = g.chi_dxy - g.chi_d

fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.0), facecolor='white')

def panel(ax, col, cmap, title, sym=False):
    v = g[col].values
    if sym:
        c = np.nanpercentile(np.abs(v), 92); vmin, vmax = -c, c
    else:
        vmin, vmax = np.nanpercentile(v, 5), np.nanpercentile(v, 95)
    sc = ax.scatter(g.n, g.delta, c=v, cmap=cmap, vmin=vmin, vmax=vmax,
                    s=S, marker='s', edgecolors='k', linewidths=0.4)
    ax.set_xlabel(r'Filling $n$', fontsize=FS_LABEL)
    ax.set_ylabel(r'Anisotropy $\delta$', fontsize=FS_LABEL)
    ax.set_title(title, fontsize=FS_TITLE, pad=10)
    ax.tick_params(labelsize=FS_TICK)
    ax.set_xlim(NMIN - 0.05, 1.05); ax.set_ylim(-0.05, 0.45)
    cb = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
    cb.ax.tick_params(labelsize=FS_CBAR)

panel(axes[0], 'chi_d',   'viridis',  r'(a)  $\chi_{d_{x^2-y^2}}$')
panel(axes[1], 'chi_dxy', 'viridis',  r'(b)  $\chi_{d_{xy}}$')
panel(axes[2], 'diff',    'coolwarm',
      r'(c)  $\chi_{d_{xy}} - \chi_{d_{x^2-y^2}}$', sym=True)

plt.tight_layout()
plt.savefig(OUT, dpi=300, bbox_inches='tight', facecolor='white')
print("saved", OUT, " points:", len(g),
      " fillings:", sorted(g.n.round(3).unique()))
