"""All four pairing channels at L = 8, 10 and 12, over (delta, U), half filling.

Rows are lattice sizes, columns are channels. ONE colour scale for all twelve
panels, so magnitudes are comparable everywhere -- extended s really is the
largest channel and the picture now says so.

The cost is real and worth knowing: the scale spans -0.16..+0.51 to cover
extended s, while on-site s only occupies -0.16..-0.06 of it. That channel
therefore reads as a nearly uniform block and its internal structure is no
longer visible. dx2-y2 loses contrast for the same reason. Use fig_4ch_3L.py,
which gives each column its own bar, when the structure WITHIN a channel is
the point.

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
LO = min(d[k].min() for k, _ in CH)              # one scale, all channels, all L
HI = max(d[k].max() for k, _ in CH)
LEV = np.linspace(LO, HI, 100)
for j, (key, lab) in enumerate(CH):
    lev = LEV
    for i, L in enumerate(LS):
        a = ax[i, j]; s = d[d.L == L]
        Z = griddata((s.delta.values, s.U.values), s[key].values, (gx, gy), method='linear')
        assert not np.isnan(Z).any(), f'L={L} {key}: hull hole'
        im = a.contourf(gx, gy, Z, levels=lev, cmap='jet', extend='both')
        if i == 0: a.set_title(lab, fontsize=17, pad=10)
        if i == len(LS)-1: a.set_xlabel(r'anisotropy  $\delta$')
        if j == 0: a.set_ylabel(rf'$L={L}$' '\n' r'$U/t$', fontsize=15)
    print(f'{key:5s} ' +
          '  '.join(f'L{L}:{d[d.L==L][key].min():+.3f}..{d[d.L==L][key].max():+.3f}' for L in LS))
cb = fig.colorbar(im, ax=ax, pad=0.012, fraction=0.014)
cb.set_label(r'pairing vertex at $\mathbf{k}=0$', fontsize=16)
cb.ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.2f'))
cb.ax.tick_params(labelsize=11)
print(f'\nsingle shared scale {LO:+.3f} .. {HI:+.3f}')
plt.savefig(f'{D}/../figs/fig_4ch_3L_shared.png', dpi=600, bbox_inches='tight', facecolor='white')
print('\ndxy - dx2y2 at half filling, U=4, by L:')
for L in LS:
    s = d[(d.L == L) & (d.U == 4)].sort_values('delta')
    print(f'  L={L:2d}: ' + '  '.join(f'{r.delta:g}:{r.dxy-r.d:+.3f}' for _, r in s.iterrows()))
