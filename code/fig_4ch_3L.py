"""All four pairing channels at L = 8, 10 and 12, over (delta, U), half filling.

Rows are lattice sizes, columns are channels. Each COLUMN shares one colour
scale across the three sizes, so reading down a column shows the finite-size
behaviour of that channel. Colours must NOT be compared between columns: the
channels differ by an order of magnitude and each has its own bar.

NORMALISATION, and it differs from Delta_tot. The Fortran vertex at k=0 is
already intensive -- extended s runs 0.186..0.500 at L=8 and 0.182..0.503 at
L=12 -- so it is plotted RAW. Dividing by L^2 here would manufacture a spurious
1/L^2 shrinkage. Delta_tot is the opposite case: it is a bare momentum sum that
grows with the number of k-points, so it must be divided by L^2 before sizes can
be compared. Two quantities, two conventions, both deliberate.

All runs are Fortran at the production settings (Delta_tau = 0.01, beta = 32,
N_walkers = 1000). The beta = 3 issue affected only the Python chi drivers.
"""
import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt
from scipy.interpolate import griddata
from gridinterp import reggrid

mpl.rcParams.update({
    'font.family':'serif', 'font.serif':['DejaVu Serif'],
    'mathtext.fontset':'dejavuserif', 'axes.linewidth':1.2,
    'xtick.direction':'out', 'ytick.direction':'out',
    'xtick.top':False, 'ytick.right':False,
    'axes.labelsize':15, 'xtick.labelsize':11, 'ytick.labelsize':11})

D = '/Users/liujiaxin/Desktop/checkerboard-altermagnet/data'
CH = [('son', r'on-site $s$'), ('sext', r'extended $s$'),
      ('d', r'$d_{x^2-y^2}$'), ('dxy', r'$d_{xy}$')]
d = pd.read_csv(f'{D}/fortran_L8_L10_L12_dedup.csv')
d = d[np.isclose(d.n, 1.0) & (d.U >= 2)]
LS = sorted(d.L.unique())

fig, ax = plt.subplots(len(LS), len(CH), figsize=(5.6*len(CH), 4.5*len(LS)),
                       facecolor='white', sharex=True, sharey=True)
gx, gy = np.meshgrid(np.linspace(d.delta.min(), d.delta.max(), 300),
                     np.linspace(d.U.min(), d.U.max(), 300))
for j, (key, lab) in enumerate(CH):
    lo, hi = d[key].min(), d[key].max()          # column-wide scale, all L
    lev = np.linspace(lo, hi, 100)
    for i, L in enumerate(LS):
        a = ax[i, j]; s = d[d.L == L]
        Z = reggrid(s, 'delta', 'U', key, gx, gy)
        assert not np.isnan(Z).any(), f'L={L} {key}: hull hole'
        im = a.contourf(gx, gy, Z, levels=lev, cmap='jet', extend='both')
        if i == 0: a.set_title(lab, fontsize=17, pad=10)
        if i == len(LS)-1: a.set_xlabel(r'anisotropy  $\delta$')
        if j == 0: a.set_ylabel(rf'$L={L}$' '\n' r'$U/t$', fontsize=15)
    cb = fig.colorbar(im, ax=ax[:, j], pad=0.015, fraction=0.020)
    cb.ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.2f'))
    cb.ax.tick_params(labelsize=10)
    print(f'{key:5s} scale {lo:+.3f}..{hi:+.3f}   ' +
          '  '.join(f'L{L}:{d[d.L==L][key].min():+.3f}..{d[d.L==L][key].max():+.3f}' for L in LS))
plt.savefig(f'{D}/../figs/fig_4ch_3L.png', dpi=600, bbox_inches='tight', facecolor='white')
print('\ndxy - dx2y2 at half filling, U=4, by L:')
for L in LS:
    s = d[(d.L == L) & (d.U == 4)].sort_values('delta')
    print(f'  L={L:2d}: ' + '  '.join(f'{r.delta:g}:{r.dxy-r.d:+.3f}' for _, r in s.iterrows()))
