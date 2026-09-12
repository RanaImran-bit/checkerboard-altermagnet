"""How each pairing channel behaves with delta and with U, coloured by Delta_tot.

  (a) x = anisotropy delta   (b) x = interaction U/t
  y  = pairing vertex at k = 0, all four channels, distinguished by marker shape
  colour = Delta_tot, the altermagnetic splitting of that same cell

This is the combination the earlier attempts could not reach. A colour map
encodes exactly ONE field, so a Delta_tot colorbar over four panels that differ
from each other is impossible: nothing would vary between them. Moving the
channel onto the MARKER frees the colour for Delta_tot, and then all four
channels, both control parameters, and the splitting fit in two panels.

Read it as: follow one marker shape across a panel to see that channel's
response, and read the colour to see whether the splitting is large there.

L = 12, half filling, complete 7 x 7 grid, U = 0 excluded (every channel and
Delta_tot vanish identically there).

Delta_tot = (1/L^2) sum_k |n_up(k) - n_dn(k)|, our per-site convention.
"""
import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.lines import Line2D

mpl.rcParams.update({
    'font.family':'serif', 'font.serif':['DejaVu Serif'],
    'mathtext.fontset':'dejavuserif', 'axes.linewidth':1.2,
    'xtick.direction':'in', 'ytick.direction':'in',
    'xtick.top':True, 'ytick.right':True, 'legend.frameon':False,
    'axes.labelsize':16, 'xtick.labelsize':12, 'ytick.labelsize':12})

D = '/Users/liujiaxin/Desktop/checkerboard-altermagnet/data'
d = pd.read_csv(f'{D}/fortran_L12_all49.csv')
d = d[np.isclose(d.n, 1.0) & (d.U >= 2)]
CH = [('son', r'on-site $s$', 'v'), ('sext', r'extended $s$', 'o'),
      ('d', r'$d_{x^2-y^2}$', 's'), ('dxy', r'$d_{xy}$', 'D')]

norm = plt.Normalize(d.dtot_N.min(), d.dtot_N.max())
fig, ax = plt.subplots(1, 2, figsize=(14.4, 5.4), facecolor='white', sharey=True)

for p, (xcol, xlab) in enumerate([('delta', r'anisotropy  $\delta$'),
                                  ('U', r'$U/t$')]):
    a = ax[p]
    for key, lab, mk in CH:
        a.scatter(d[xcol], d[key], c=d.dtot_N, cmap='jet', norm=norm,
                  s=110, marker=mk, edgecolors='k', linewidths=0.6, zorder=3)
    a.axhline(0, color='k', lw=1.0, ls='--')
    a.set_xlabel(xlab)
    if p == 0: a.set_ylabel(r'pairing vertex at $\mathbf{k}=0$')
    a.set_title(f'({chr(97+p)})', loc='left', fontsize=15)
    a.grid(alpha=0.18, lw=0.7); a.set_axisbelow(True)

ax[0].legend(handles=[Line2D([], [], marker=mk, ls='none', mfc='0.85', mec='k',
                             ms=9, label=lab) for _, lab, mk in CH],
             fontsize=11, loc='center left')

sm = plt.cm.ScalarMappable(cmap='jet', norm=norm); sm.set_array([])
cb = fig.colorbar(sm, ax=ax, pad=0.015, fraction=0.03)
cb.set_label(r'$\Delta_{tot}$', fontsize=17)
cb.ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3f'))
cb.ax.tick_params(labelsize=11)
plt.savefig(f'{D}/../figs/fig_pair_vs_delta_U_dtot.png', dpi=600,
            bbox_inches='tight', facecolor='white')

print(f'{len(d)} cells, {len(CH)} channels, Delta_tot '
      f'{d.dtot_N.min():.4f} to {d.dtot_N.max():.4f}')
for key, lab, _ in CH:
    rd = np.mean([np.corrcoef(g.delta, g[key])[0,1] for _, g in d.groupby('U')])
    ru = np.mean([np.corrcoef(g.U, g[key])[0,1] for _, g in d.groupby('delta')])
    print(f'  {lab:16s} vs delta r={rd:+.2f}   vs U r={ru:+.2f}')
