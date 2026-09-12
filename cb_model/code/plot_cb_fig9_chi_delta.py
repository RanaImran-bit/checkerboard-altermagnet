#!/usr/bin/env python3
"""Fig 9: pairing-susceptibility connected vertex chi_d, chi_dxy vs anisotropy delta
(fixed filling). Reads a checkerboard_chi_scan.py CSV (default cb_chi_6x6.csv).
  python pyqmc/plot_cb_fig9_chi_delta.py cb_chi_6x6.csv
"""
import sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

CSV = sys.argv[1] if len(sys.argv) > 1 else "cb_chi_6x6.csv"
d = pd.read_csv(CSV, sep=r'\s+', engine='python', on_bad_lines='skip')
d = d[pd.to_numeric(d['delta'], errors='coerce').notna()].copy()
for c in d.columns:
    d[c] = d[c].astype(float)
n = 2 * d['nup'].iloc[0] / (d['lx'].iloc[0] * d['ly'].iloc[0])
L = int(d['lx'].iloc[0])
g = d.groupby('delta').agg(d=('chi_d_vtx', 'mean'), de=('chi_d_vtx', 'std'),
                           x=('chi_dxy_vtx', 'mean'), xe=('chi_dxy_vtx', 'std')
                           ).reset_index().sort_values('delta')
fig, ax = plt.subplots(figsize=(7, 5.5), facecolor='white')
ax.errorbar(g.delta, g.d, yerr=g.de.fillna(0), marker='o', ms=10, lw=2, capsize=4,
            color='#1f77b4', label=r'$d_{x^2-y^2}$')
ax.errorbar(g.delta, g.x, yerr=g.xe.fillna(0), marker='^', ms=10, lw=2, capsize=4,
            color='#2ca02c', label=r'$d_{xy}$')
ax.axhline(0, color='gray', ls='--', lw=1)
ax.set_xlabel(r'Anisotropy $\delta$', fontsize=18)
ax.set_ylabel(r'$\chi^{\mathrm{vertex}}_{\zeta}$ (pairing susceptibility)', fontsize=16)
ax.set_title(rf'${L}\times{L}$,  $n={n:.2f}$,  $U=4$', fontsize=16)
ax.tick_params(labelsize=13)
ax.legend(fontsize=15, frameon=False)
plt.tight_layout()
plt.savefig("Fig9_chi_vs_delta.png", dpi=300, bbox_inches='tight', facecolor='white')
print("saved Fig9_chi_vs_delta.png")
