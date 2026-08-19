"""Delta_tot over (delta, U), L = 12, half filling. Complete 7 x 7 grid.

Single panel: with the pairing contours AND the measured-point markers removed
there is only one field left to show, so a four-panel row would have been four
copies of the same picture.

Nothing marks the measured cells any more, so state the grid in the caption:
7 interaction values x 7 anisotropies, every cell measured, U = 0 dropped.

The grid has no holes, so griddata interpolates between neighbours everywhere
and no convex hull appears. This is the difference from the older L = 14 figure,
where 19 of 35 cells left a white diamond. It also reaches delta = 0.7.

NORMALISATION. Delta_tot is the momentum sum divided by the number of sites:

    Delta_tot = (1/L^2) sum_k |n_up(k) - n_dn(k)|

The absolute value is inside the sum, so it stays finite although the net
magnetisation vanishes. The 1/L^2 is OUR per-site convention, chosen so that
different L are comparable. Earlier versions of these scripts attributed it to
the reference paper's Eq. (2); that attribution was not verified and should not
be repeated. Quote the source equation verbatim in the manuscript and define
this one separately.

U = 0 is excluded: it is identically zero at every delta, and the U grid is
uneven below 2, so shading would smooth across a wide unsampled stretch.
"""
import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt
from scipy.interpolate import griddata

mpl.rcParams.update({
    'font.family':'serif', 'font.serif':['DejaVu Serif'],
    'mathtext.fontset':'dejavuserif', 'axes.linewidth':1.2,
    'xtick.direction':'out', 'ytick.direction':'out',
    'xtick.top':False, 'ytick.right':False,
    'axes.labelsize':17, 'xtick.labelsize':13, 'ytick.labelsize':13})

D = '/Users/liujiaxin/Desktop/checkerboard-altermagnet/data'
d = pd.read_csv(f'{D}/fortran_L12_all49.csv')
d = d[np.isclose(d.n, 1.0) & (d.U >= 2)]

gx, gy = np.meshgrid(np.linspace(d.delta.min(), d.delta.max(), 300),
                     np.linspace(d.U.min(), d.U.max(), 300))
Z = griddata((d.delta.values, d.U.values), d.dtot_N.values, (gx, gy), method='linear')
assert not np.isnan(Z).any(), 'hull hole -- grid is not complete'

fig, a = plt.subplots(figsize=(7.6, 5.6), facecolor='white')
im = a.contourf(gx, gy, Z, levels=100, cmap='jet', extend='both')
cb = fig.colorbar(im, ax=a, pad=0.02)
cb.set_label(r'$\Delta_{tot}$', fontsize=17)
cb.ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3f'))
cb.ax.tick_params(labelsize=11)
a.set_xlabel(r'anisotropy  $\delta$'); a.set_ylabel(r'$U/t$')
plt.tight_layout()
plt.savefig(f'{D}/../figs/fig_dtot_map_L12.png', dpi=600, bbox_inches='tight',
            facecolor='white')

print(f'Delta_tot {d.dtot_N.min():.4f} to {d.dtot_N.max():.4f}, {len(d)} measured cells')
print('ridge -- delta of maximum at each U:')
for U in sorted(d.U.unique()):
    s = d[d.U == U]
    print(f'   U={U:g}: delta={s.loc[s.dtot_N.idxmax(), "delta"]:.1f}  ({s.dtot_N.max():.4f})')
