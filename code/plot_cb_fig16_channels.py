#!/usr/bin/env python3
"""Fig 16: four-channel pairing susceptibility analysis (on-site s / extended-s / dx2-y2 / dxy).
Reads chi_grid_all.csv (cols: U,nup,delta,seed,chi_son,chi_sext,chi_d,chi_dxy,n) and writes:
  fig_4ch_maps_U{U}.png          - (n,delta) phase maps, all 4 channels, at U=UMAP
  fig_4ch_vsU.png                - chi vs U line graphs, 4 channels, at 2 fillings
  fig_dominant_channel_U{U}.png  - leading (most attractive) channel incl. extended-s
  fig_unconventional_channel_U{U}.png - leading unconventional channel (dx2-y2 vs dxy)
Usage:
  python plot_cb_fig16_channels.py chi_grid_all.csv
"""
import sys
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# ---- EDIT HERE ----
CSV  = sys.argv[1] if len(sys.argv) > 1 else "chi_grid_all.csv"
NMIN = 0.55            # trim dilute limit
UMAP = 8.0             # U for the (n,delta) map + dominant-channel panels
DLINE = 0.4            # delta for the vs-U line graphs
FILLINGS_LINE = [0.778, 1.0]
S = 360
COL = ['#7f7f7f', '#2ca02c', '#1f77b4', '#d62728']   # on-site s, ext-s, dx2-y2, dxy
CH  = [('chi_son', 'on-site s'), ('chi_sext', 'extended-s'),
       ('chi_d', r'$d_{x^2-y^2}$'), ('chi_dxy', r'$d_{xy}$')]
# -------------------

d = pd.read_csv(CSV); d = d[d.n >= NMIN]
g = d.groupby(['U', 'n', 'delta'])[['chi_son', 'chi_sext', 'chi_d', 'chi_dxy']].mean().reset_index()
g8 = g[np.isclose(g.U, UMAP)]

# --- Fig A: 4-channel (n,delta) phase maps ---
fig, ax = plt.subplots(1, 4, figsize=(20, 4.7), facecolor='white')
for a, (col, lab) in zip(ax, CH):
    v = g8[col].values; vmin, vmax = np.nanpercentile(v, 5), np.nanpercentile(v, 95)
    sc = a.scatter(g8.n, g8.delta, c=v, cmap='viridis', vmin=vmin, vmax=vmax,
                   s=S, marker='s', edgecolors='k', linewidths=0.4)
    a.set_title(lab, fontsize=18); a.set_xlabel('Filling $n$', fontsize=15)
    a.set_ylabel(r'$\delta$', fontsize=15); a.set_xlim(0.5, 1.05); a.set_ylim(-0.05, 0.45)
    a.tick_params(labelsize=12); fig.colorbar(sc, ax=a, fraction=0.046, pad=0.04)
fig.suptitle(fr'Connected-vertex pairing susceptibility over $(n,\delta)$ at $U={UMAP:g}$',
             fontsize=17, y=1.02)
plt.tight_layout(); plt.savefig(f'fig_4ch_maps_U{UMAP:g}.png', dpi=600, bbox_inches='tight', facecolor='white'); plt.close()

# --- Fig B: chi vs U line graphs ---
dd = d[np.isclose(d.delta, DLINE)]
fig, ax = plt.subplots(1, len(FILLINGS_LINE), figsize=(6.5*len(FILLINGS_LINE), 5.2), facecolor='white')
for a, nsel in zip(np.atleast_1d(ax), FILLINGS_LINE):
    sub = dd[np.isclose(dd.n, nsel, atol=0.02)]
    for (col, lab), c in zip(CH, COL):
        s = sub.groupby('U')[col].agg(['mean', 'sem']).reset_index().sort_values('U')
        a.errorbar(s.U, s['mean'], yerr=s['sem'].fillna(0), marker='o', ms=8, lw=2, capsize=4, color=c, label=lab)
    a.axhline(0, color='gray', ls='--', lw=0.8); a.set_title(fr'$n={nsel:g}$, $\delta={DLINE:g}$', fontsize=15)
    a.set_xlabel('Interaction $U$', fontsize=15); a.set_ylabel(r'$\chi^{\rm vertex}$', fontsize=15)
    a.legend(fontsize=12, frameon=False)
plt.tight_layout(); plt.savefig('fig_4ch_vsU.png', dpi=600, bbox_inches='tight', facecolor='white'); plt.close()

# --- Fig C: leading channel (all 4) ---
cols = ['chi_son', 'chi_sext', 'chi_d', 'chi_dxy']; labs = ['on-site s', 'extended-s', 'dx2-y2', 'dxy']
win = np.argmax(g8[cols].values, axis=1)
fig, a = plt.subplots(figsize=(7.5, 6), facecolor='white')
for i in range(len(g8)):
    a.scatter(g8.n.values[i], g8.delta.values[i], c=COL[win[i]], s=460, marker='s', edgecolors='k', linewidths=0.5)
a.set_xlabel('Filling $n$', fontsize=16); a.set_ylabel(r'Anisotropy $\delta$', fontsize=16)
a.set_title(fr'Leading (most attractive) pairing channel, $U={UMAP:g}$', fontsize=15)
a.set_xlim(0.5, 1.05); a.set_ylim(-0.05, 0.45); a.tick_params(labelsize=13)
a.legend(handles=[Patch(facecolor=COL[i], edgecolor='k', label=labs[i]) for i in range(4)],
         fontsize=12, frameon=False, loc='center left', bbox_to_anchor=(1.02, 0.5))
plt.tight_layout(); plt.savefig(f'fig_dominant_channel_U{UMAP:g}.png', dpi=600, bbox_inches='tight', facecolor='white'); plt.close()

# --- Fig D: leading UNCONVENTIONAL channel (dx2-y2 vs dxy) ---
w2 = np.where(g8.chi_d.values >= g8.chi_dxy.values, 0, 1); c2 = ['#1f77b4', '#d62728']; l2 = ['dx2-y2', 'dxy']
fig, a = plt.subplots(figsize=(7.5, 6), facecolor='white')
for i in range(len(g8)):
    a.scatter(g8.n.values[i], g8.delta.values[i], c=c2[w2[i]], s=460, marker='s', edgecolors='k', linewidths=0.5)
a.set_xlabel('Filling $n$', fontsize=16); a.set_ylabel(r'Anisotropy $\delta$', fontsize=16)
a.set_title(fr'Leading unconventional channel, $U={UMAP:g}$', fontsize=15)
a.set_xlim(0.5, 1.05); a.set_ylim(-0.05, 0.45); a.tick_params(labelsize=13)
a.legend(handles=[Patch(facecolor=c2[i], edgecolor='k', label=l2[i]) for i in range(2)],
         fontsize=13, frameon=False, loc='center left', bbox_to_anchor=(1.02, 0.5))
plt.tight_layout(); plt.savefig(f'fig_unconventional_channel_U{UMAP:g}.png', dpi=600, bbox_inches='tight', facecolor='white'); plt.close()

print("saved: fig_4ch_maps_U%g.png, fig_4ch_vsU.png, fig_dominant_channel_U%g.png, fig_unconventional_channel_U%g.png"
      % (UMAP, UMAP, UMAP))
