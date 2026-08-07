#!/usr/bin/env python3
"""Point 7: the altermagnetic order parameter over the (U, delta) plane.

Separates what U contributes from what delta contributes. m_AM = Delta m * delta with
Delta m = m - m(U=0), so it vanishes if either ingredient is absent. Heatmaps rather than
scatter, per the publication-figure standard.

  python plot_mAM_UD.py [data/mAM_L12_grid.csv]
"""
import sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

CSV = sys.argv[1] if len(sys.argv) > 1 else "data/mAM_L12_grid.csv"
FILLINGS = [0.6667, 0.7778, 0.8889, 1.0]
g = pd.read_csv(CSV)
US = np.sort(g.U.unique()); DS = np.sort(g.delta.unique())

def grid(n, col):
    s = g[np.isclose(g.n, n, atol=0.01)]
    return s.pivot(index='delta', columns='U', values=col).reindex(index=DS, columns=US).values

fig, ax = plt.subplots(2, len(FILLINGS), figsize=(4.3*len(FILLINGS), 8.0), facecolor='white')
for j, n in enumerate(FILLINGS):
    for i, (col, lab, cmap) in enumerate([('dm', r'$\Delta m = m - m_{U=0}$', 'viridis'),
                                          ('mAM', r'$m_{\rm AM} = \Delta m\,\delta$', 'magma')]):
        Z = grid(n, col)
        a = ax[i, j]
        im = a.pcolormesh(US, DS, Z, cmap=cmap, shading='gouraud', vmin=0)
        a.set_xticks(US); a.set_yticks(DS)
        a.set_xlabel('interaction $U$', fontsize=13)
        if j == 0: a.set_ylabel(r'anisotropy $\delta$', fontsize=14)
        a.set_title(rf'$n={n:.3f}$', fontsize=13) if i == 0 else None
        cb = fig.colorbar(im, ax=a, fraction=0.046, pad=0.03)
        cb.ax.tick_params(labelsize=9)
        if j == len(FILLINGS)-1: cb.set_label(lab, fontsize=12)
        a.tick_params(labelsize=11)
ax[0,0].text(-0.42, 0.5, 'the moment\n(U creates it,\n$\\delta$ suppresses it)', transform=ax[0,0].transAxes,
             rotation=90, va='center', ha='center', fontsize=12, color='#1E2761')
ax[1,0].text(-0.42, 0.5, 'the AM order parameter\n(needs BOTH)', transform=ax[1,0].transAxes,
             rotation=90, va='center', ha='center', fontsize=12, color='#B85042')
fig.suptitle(r'What $U$ and $\delta$ each contribute to the altermagnetism  ($L=12$)',
             fontsize=17, y=0.98)
plt.tight_layout(rect=[0.03, 0, 1, 0.96])
plt.savefig('figures/fig_mAM_U_delta_L12.png', dpi=600, bbox_inches='tight', facecolor='white')
print("saved figures/fig_mAM_U_delta_L12.png")

# ---- the separation, as line cuts ----
fig2, ax2 = plt.subplots(1, 3, figsize=(16, 4.6), facecolor='white')
cmn = plt.cm.viridis(np.linspace(0, .85, len(FILLINGS)))
for n, c in zip(FILLINGS, cmn):
    s = g[np.isclose(g.n, n, atol=0.01)]
    d0 = s[np.isclose(s.delta, 0.0)].sort_values('U')
    ax2[0].plot(d0.U, d0.dm, '-o', ms=8, lw=2.2, color=c, label=f'n={n:.3f}')
    u8 = s[np.isclose(s.U, 8.0)].sort_values('delta')
    ax2[1].plot(u8.delta, u8.dm, '-o', ms=8, lw=2.2, color=c)
    ax2[2].plot(u8.delta, u8.mAM, '-o', ms=8, lw=2.2, color=c)
ax2[0].set_xlabel('interaction $U$', fontsize=14); ax2[0].set_ylabel(r'$\Delta m$', fontsize=15)
ax2[0].set_title(r'(a) $U$ creates the moment   ($\delta=0$)', fontsize=14); ax2[0].legend(fontsize=11, frameon=False)
ax2[1].set_xlabel(r'anisotropy $\delta$', fontsize=14); ax2[1].set_ylabel(r'$\Delta m$', fontsize=15)
ax2[1].set_title(r'(b) $\delta$ suppresses it   ($U=8$)', fontsize=14)
ax2[2].set_xlabel(r'anisotropy $\delta$', fontsize=14); ax2[2].set_ylabel(r'$m_{\rm AM}$', fontsize=15)
ax2[2].set_title(r'(c) the product still grows   ($U=8$)', fontsize=14)
for a in ax2: a.tick_params(labelsize=12)
plt.tight_layout()
plt.savefig('figures/fig_mAM_separation_L12.png', dpi=600, bbox_inches='tight', facecolor='white')
print("saved figures/fig_mAM_separation_L12.png")
