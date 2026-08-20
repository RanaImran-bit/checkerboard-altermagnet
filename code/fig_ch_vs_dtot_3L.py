"""Pairing channels against anisotropy, coloured by Delta_tot, at three sizes.

Columns are the four channels, rows are L = 8, 10, 12. In every panel:

    x = anisotropy delta
    y = pairing vertex at k = 0 for that channel
    colour = Delta_tot/N of the SAME cell

so one point is one (U, delta) cell and the colour says how large the
altermagnetic splitting is where that pairing value was measured. ONE colour
scale for all twelve panels, because colour is a single quantity throughout.

This is the only arrangement that shows Delta_tot as colour AND keeps the four
channels distinguishable. Putting Delta_tot on the colour axis of a (delta, U)
map makes all four columns identical, since Delta_tot does not depend on the
channel -- the channel has to occupy an axis.

Delta_tot/N is per site (bare Eq. 3 sum divided by L^2), which is the only form
comparable across L. The pairing vertex is raw: it is already intensive.
"""
import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt

mpl.rcParams.update({
    'font.family':'serif', 'font.serif':['DejaVu Serif'],
    'mathtext.fontset':'dejavuserif', 'axes.linewidth':1.2,
    'xtick.direction':'in', 'ytick.direction':'in',
    'xtick.top':True, 'ytick.right':True,
    'axes.labelsize':15, 'xtick.labelsize':11, 'ytick.labelsize':11})

D = '/Users/liujiaxin/Desktop/checkerboard-altermagnet/data'
CH = [('son', r'on-site $s$'), ('sext', r'extended $s$'),
      ('d', r'$d_{x^2-y^2}$'), ('dxy', r'$d_{xy}$')]
d = pd.read_csv(f'{D}/fortran_L8_L10_L12_dedup.csv')
d = d[np.isclose(d.n, 1.0) & (d.U >= 2)]
LS = sorted(d.L.unique())
norm = plt.Normalize(d.dtot_N.min(), d.dtot_N.max())

fig, ax = plt.subplots(len(LS), len(CH), figsize=(5.2*len(CH), 4.2*len(LS)),
                       facecolor='white', sharex=True)
for j, (key, lab) in enumerate(CH):
    ylo, yhi = d[key].min(), d[key].max(); pad = 0.06*(yhi-ylo)
    for i, L in enumerate(LS):
        a = ax[i, j]; s = d[d.L == L]
        sc = a.scatter(s.delta, s[key], c=s.dtot_N, cmap='jet', norm=norm,
                       s=105, edgecolors='k', linewidths=0.5, zorder=3)
        a.axhline(0, color='k', lw=0.9, ls='--')
        a.set_ylim(ylo-pad, yhi+pad)          # same y range down a column
        a.grid(alpha=0.18, lw=0.7); a.set_axisbelow(True)
        if i == 0: a.set_title(lab, fontsize=17, pad=10)
        if i == len(LS)-1: a.set_xlabel(r'anisotropy  $\delta$')
        if j == 0: a.set_ylabel(rf'$L={L}$' '\n' r'vertex at $\mathbf{k}=0$', fontsize=14)
cb = fig.colorbar(sc, ax=ax, pad=0.012, fraction=0.014)
cb.set_label(r'$\Delta_{tot}/N$', fontsize=17)
cb.ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3f'))
cb.ax.tick_params(labelsize=11)
plt.savefig(f'{D}/../figs/fig_ch_vs_dtot_3L.png', dpi=600, bbox_inches='tight',
            facecolor='white')

print(f'colour scale (Delta_tot/N): {d.dtot_N.min():.4f} .. {d.dtot_N.max():.4f}')
print('\ncorrelation of each channel with Delta_tot/N, within fixed U, averaged:')
for key, lab in CH:
    print(f'  {lab:16s} ' + '  '.join(
        f'L{L}: {np.mean([np.corrcoef(g.dtot_N, g[key])[0,1] for _, g in d[d.L==L].groupby("U")]):+.2f}'
        for L in LS))
