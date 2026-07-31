#!/usr/bin/env python3
"""APBC-vs-periodic half-filling comparison (the open-shell resolution).
Overlays chi_dxy, chi_dx2-y2, and the staggered moment M vs delta at half-filling for
anti-periodic BC (APBC-x, trial-clean) vs periodic BC (open-shell), one line per U.

Result: the dxy enhancement with delta is CLEAN and monotonic at all U under APBC (the
trial-clean trial), confirming it is physical; the periodic (open-shell) version shows the
same trend but is noisy at weak U and has contaminated magnitudes -> open-shell caveat removed.

  python plot_cb_apbc_compare.py apbc_half_APX-1.csv apbc_half_APX1.csv
"""
import sys
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

APBC = sys.argv[1] if len(sys.argv) > 1 else "apbc_half_APX-1.csv"   # trial-clean
PER  = sys.argv[2] if len(sys.argv) > 2 else "apbc_half_APX1.csv"    # open-shell
gA = pd.read_csv(APBC).groupby(['U', 'delta'])[['chi_dxy', 'chi_d', 'M']].mean().reset_index()
gP = pd.read_csv(PER ).groupby(['U', 'delta'])[['chi_dxy', 'chi_d', 'M']].mean().reset_index()

cols = [('chi_dxy', r'$\chi_{d_{xy}}$'), ('chi_d', r'$\chi_{d_{x^2-y^2}}$'), ('M', 'staggered moment $M$')]
colors = {2: '#ff7f0e', 4: '#2ca02c', 6: '#d62728', 8: '#9467bd'}
fig, ax = plt.subplots(1, 3, figsize=(17, 5), facecolor='white')
for a, (c, lab) in zip(ax, cols):
    for U in [2, 4, 6, 8]:
        sa = gA[gA.U == U].sort_values('delta'); sp = gP[gP.U == U].sort_values('delta')
        a.plot(sa.delta, sa[c], '-o', color=colors[U], lw=2, label=f'U={U} APBC')
        a.plot(sp.delta, sp[c], '--s', color=colors[U], lw=1.5, alpha=0.6, label=f'U={U} periodic')
    a.axhline(0, color='gray', ls=':', lw=0.8)
    a.set_xlabel(r'anisotropy $\delta$', fontsize=14); a.set_ylabel(lab, fontsize=14)
    a.set_title(f'{lab} vs δ (n=1)', fontsize=13)
ax[0].legend(fontsize=8, ncol=2, frameon=False)
fig.suptitle('Half-filling: APBC (solid, trial-clean) vs periodic (dashed, open-shell)', fontsize=14, y=1.02)
plt.tight_layout()
plt.savefig('fig_apbc_vs_periodic.png', dpi=200, bbox_inches='tight', facecolor='white')
print("saved fig_apbc_vs_periodic.png")
