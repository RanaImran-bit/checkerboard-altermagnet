#!/usr/bin/env python3
"""Fig 13: magnetism-vs-pairing overlay (twin axis). Left axis = chi_zz(pi,pi) (magnetic);
right axis = chi_{dx2-y2} connected-vertex (pairing), both vs filling n at fixed delta.
  python plot_cb_fig13_mag_vs_pair.py ../data/chizz_vs_n.csv ../data/chi_n_d04.csv
"""
import sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- EDIT HERE: style ----
FS_LABEL = 18
FS_TICK  = 12
FS_TITLE = 15
FS_LEG   = 13
MS       = 9
C_MAG    = '#8c2d8c'    # magnetic (left axis)
C_PAIR   = '#1f77b4'    # dx2-y2 pairing (right axis)
TITLE    = r'magnetism vs pairing, $\delta=0.4$'
OUT      = "Fig13_mag_vs_pair.png"
# --------------------------

CZ  = sys.argv[1] if len(sys.argv) > 1 else "../data/chizz_vs_n.csv"
CHI = sys.argv[2] if len(sys.argv) > 2 else "../data/chi_n_d04.csv"

cz = pd.read_csv(CZ).sort_values('n')
d = pd.read_csv(CHI, sep=r'\s+', engine='python', on_bad_lines='skip')
d = d[pd.to_numeric(d['delta'], errors='coerce').notna()].copy()
for c in d.columns:
    d[c] = d[c].astype(float)
d['n'] = 2 * d['nup'] / (d['lx'] * d['ly'])
gp = d.groupby('n').agg(dd=('chi_d_vtx', 'mean'), dde=('chi_d_vtx', 'sem')).reset_index().sort_values('n')

fig, ax = plt.subplots(figsize=(7.5, 5.5), facecolor='white')
l1 = ax.errorbar(cz.n, cz.chizz_pipi, yerr=cz.err, marker='s', ms=MS, lw=2, capsize=4,
                 color=C_MAG, label=r'$\chi_{zz}(\pi,\pi)$ (magnetic)')
ax.set_xlabel(r'Filling $n$', fontsize=FS_LABEL)
ax.set_ylabel(r'$\chi_{zz}(\pi,\pi)$', fontsize=FS_LABEL - 2, color=C_MAG)
ax.tick_params(axis='y', labelcolor=C_MAG, labelsize=FS_TICK)
ax.tick_params(axis='x', labelsize=FS_TICK + 1)

ax2 = ax.twinx()
l2 = ax2.errorbar(gp.n, gp.dd, yerr=gp.dde.fillna(0), marker='o', ms=MS, lw=2, capsize=4,
                  color=C_PAIR, label=r'$\chi_{d_{x^2-y^2}}$ (pairing)')
ax2.axhline(0, color='gray', ls='--', lw=0.8)
ax2.set_ylabel(r'$\chi^{\mathrm{vertex}}_{d_{x^2-y^2}}$', fontsize=FS_LABEL - 2, color=C_PAIR)
ax2.tick_params(axis='y', labelcolor=C_PAIR, labelsize=FS_TICK)

ax.set_title(TITLE, fontsize=FS_TITLE)
lns = [l1, l2]
ax.legend(lns, [x.get_label() for x in lns], fontsize=FS_LEG, frameon=False, loc='upper center')
plt.tight_layout()
plt.savefig(OUT, dpi=300, bbox_inches='tight', facecolor='white')
print("saved", OUT)
