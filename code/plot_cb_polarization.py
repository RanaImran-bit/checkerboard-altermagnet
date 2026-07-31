#!/usr/bin/env python3
"""Polarization / channel-vs-delta analysis (meeting point #11).
Merges pol_grid_all.csv (staggered moment M) with chi_grid_all.csv (pairing susceptibility),
and shows the HONEST per-filling delta-dependence -- no pooled correlations.

Main figure fig_pol_chi_vs_delta.png (2x2): rows = filling (half-filling n=1, doped n=0.778),
cols = channel (dxy, dx2-y2), one line per U vs delta. Reads:
  chi_grid_all.csv (U,nup,delta,seed,chi_son,chi_sext,chi_d,chi_dxy,n)
  pol_grid_all.csv (U,nup,delta,seed,Spipi,M,n)
Usage:
  python plot_cb_polarization.py pol_grid_all.csv chi_grid_all.csv

Result (honest, no collapse):
  doped n=0.778   : delta SUPPRESSES chi_dxy (all U); chi_dx2-y2 is delta-blind  -> decoupling
  half-filling n=1: delta ENHANCES chi_dxy (clear at U=6,8); SUPPRESSES chi_dx2-y2 -> channel swap
  CAVEAT: half-filling is open-shell (degenerate free-electron trial) -> trial-dependent; U=2,4 noisy.
"""
import sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

POL = sys.argv[1] if len(sys.argv) > 1 else "pol_grid_all.csv"
CHI = sys.argv[2] if len(sys.argv) > 2 else "chi_grid_all.csv"
L = 6; N = L * L
FILLINGS = [1.0, 0.778]                       # rows: half-filling, doped
US = [0, 2, 4, 6, 8]

pol = pd.read_csv(POL); chi = pd.read_csv(CHI)
P = pol.groupby(['U', 'nup', 'delta'])['M'].mean().reset_index()
C = chi.groupby(['U', 'nup', 'delta'])[['chi_d', 'chi_dxy']].mean().reset_index()
m = P.merge(C, on=['U', 'nup', 'delta']); m['n'] = 2 * m['nup'] / N

CH = [('chi_dxy', r'$\chi_{d_{xy}}$'), ('chi_d', r'$\chi_{d_{x^2-y^2}}$')]
fig, axes = plt.subplots(len(FILLINGS), 2, figsize=(13, 5 * len(FILLINGS)), facecolor='white')
for row, nsel in enumerate(FILLINGS):
    for col, (ch, lab) in enumerate(CH):
        a = axes[row, col]
        for U in US:
            su = m[(m.U == U) & (np.isclose(m.n, nsel, atol=0.02))].sort_values('delta')
            a.plot(su.delta, su[ch], marker='o', lw=2, label=f'U={U}')
        a.axhline(0, color='gray', ls='--', lw=0.8)
        a.set_xlabel(r'anisotropy $\delta$', fontsize=14); a.set_ylabel(lab, fontsize=14)
        a.set_title(fr'{lab} vs $\delta$   (n={nsel:g})', fontsize=14); a.legend(fontsize=10)
plt.tight_layout()
plt.savefig('fig_pol_chi_vs_delta.png', dpi=300, bbox_inches='tight', facecolor='white')
print("saved fig_pol_chi_vs_delta.png")
