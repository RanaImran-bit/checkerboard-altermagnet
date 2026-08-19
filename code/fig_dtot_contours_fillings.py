"""Delta_tot colour field with pairing contours, one ROW per filling.

Rows are fillings, columns are the four pairing channels. In every panel the
colour is Delta_tot and the black lines are that channel's vertex, so the
magnetic and the pairing response can be read against each other cell by cell,
and now also across doping.

ONLY TWO FILLINGS EXIST. The figure needs a full (delta, U) grid of BOTH
Delta_tot and all four channels at each filling, and that has been run at
n = 1.000 and n = 0.847 only. No other filling in any dataset has more than a
scattered handful of cells. Adding a third and fourth row means launching two
more 49-cell Fortran grids.

DELTA RANGE. n = 1.000 is complete to delta = 0.7. n = 0.847 is missing the
delta = 0.7 column for U >= 3, which was still running when this was written.
Both rows are therefore drawn over the COMMON complete range, delta = 0.1 to
0.6, so the two rows share axes and neither is interpolated across a hole.
Re-run once those five cells land and DMAX can go to 0.7.

Colour scale is shared across BOTH rows so the doping dependence of Delta_tot is
visible. Contour levels are per panel, since the channels differ in range.
"""
import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt
from scipy.interpolate import griddata

mpl.rcParams.update({
    'font.family':'serif', 'font.serif':['DejaVu Serif'],
    'mathtext.fontset':'dejavuserif', 'axes.linewidth':1.2,
    'xtick.direction':'out', 'ytick.direction':'out',
    'xtick.top':False, 'ytick.right':False,
    'axes.labelsize':15, 'xtick.labelsize':11, 'ytick.labelsize':11})

D = '/Users/liujiaxin/Desktop/checkerboard-altermagnet/data'
DMAX = 0.6                      # common complete range; raise to 0.7 when n=0.847 finishes
CH = [('son', r'on-site $s$'), ('sext', r'extended $s$'),
      ('d', r'$d_{x^2-y^2}$'), ('dxy', r'$d_{xy}$')]
SRC = [('fortran_L12_all49.csv', 1.000), ('fortran_L12_n0847.csv', 0.847)]

rows = []
for f, nv in SRC:
    t = pd.read_csv(f'{D}/{f}')
    t = t[(t.U >= 2) & (t.delta <= DMAX + 1e-9)]
    cov = t.pivot_table(index='U', columns='delta', values='dtot_N')
    assert not cov.isna().any().any(), f'{f}: hole in the grid at delta<={DMAX}'
    rows.append((nv, t))

VMIN = min(t.dtot_N.min() for _, t in rows)
VMAX = max(t.dtot_N.max() for _, t in rows)
LEV = np.linspace(VMIN, VMAX, 100)

fig, ax = plt.subplots(len(rows), 4, figsize=(21.0, 4.6*len(rows)),
                       facecolor='white', sharex=True, sharey=True, squeeze=False)
for r, (nv, t) in enumerate(rows):
    gx, gy = np.meshgrid(np.linspace(t.delta.min(), t.delta.max(), 300),
                         np.linspace(t.U.min(), t.U.max(), 300))
    BG = griddata((t.delta.values, t.U.values), t.dtot_N.values, (gx, gy), method='linear')
    for c, (key, lab) in enumerate(CH):
        a = ax[r, c]
        im = a.contourf(gx, gy, BG, levels=LEV, cmap='jet', extend='both')
        F = griddata((t.delta.values, t.U.values), t[key].values, (gx, gy), method='linear')
        cs = a.contour(gx, gy, F, levels=6, colors='k', linewidths=1.2)
        a.clabel(cs, inline=True, fontsize=8, fmt='%.2f')
        if r == 0: a.set_title(f'contours: {lab}', fontsize=15)
        if r == len(rows)-1: a.set_xlabel(r'anisotropy  $\delta$')
        if c == 0: a.set_ylabel(rf'$n={nv:.3f}$' '\n' r'$U/t$', fontsize=14)
    print(f'n={nv:.3f}: Delta_tot {t.dtot_N.min():.4f} to {t.dtot_N.max():.4f}, {len(t)} cells')

cb = fig.colorbar(im, ax=ax, pad=0.012, fraction=0.014)
cb.set_label(r'$\Delta_{tot}$', fontsize=16)
cb.ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3f'))
cb.ax.tick_params(labelsize=10)
plt.savefig(f'{D}/../figs/fig_dtot_contours_fillings.png', dpi=600,
            bbox_inches='tight', facecolor='white')
print(f'shared colour scale {VMIN:.4f} to {VMAX:.4f}')
