"""All four pairing channels, unequal time only, over (n, delta) at L = 12.

Dropping the equal-time row frees the figure from the L = 10 eqtime files,
which reach delta = 0.7 for only 6 fillings. dense_L12_U4_all.csv is a COMPLETE
19 filling x 8 anisotropy grid at 6 seeds per cell, so the map covers the full
delta range at every filling, at L = 12 rather than L = 10, with no holes for
the interpolator to invent values across.

Normalisation. With only one time sector in the figure there is nothing to make
commensurate, so chi is divided by N alone and labelled per site. It still
carries its factor of imaginary time; that only mattered when an equal-time row
had to sit beside it.

Each panel keeps its own colour scale. The channels differ by more than an order
of magnitude, so colours must NOT be read across panels, only the numbers.
"""
import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import gaussian_filter

mpl.rcParams.update({
    'font.family':'serif', 'font.serif':['DejaVu Serif'],
    'mathtext.fontset':'dejavuserif', 'axes.linewidth':1.2,
    'xtick.direction':'out', 'ytick.direction':'out',
    'xtick.top': False, 'ytick.right': False,
    'axes.labelsize':16, 'xtick.labelsize':12, 'ytick.labelsize':12})

D = '/Users/liujiaxin/Desktop/checkerboard-altermagnet/data'
NSITES = 144                      # L = 12
BLUR = 4.0
CHAN = [('chi_son', r'on-site $s$'), ('chi_sext', r'extended $s$'),
        ('chi_d', r'$d_{x^2-y^2}$'), ('chi_dxy', r'$d_{xy}$')]

d = pd.read_csv(f'{D}/dense_L12_U4_all.csv')
NS = np.sort(d.n.unique()); DS = np.sort(d.delta.unique())
cov = d.groupby(['n', 'delta']).size().unstack()
assert not cov.isna().any().any(), 'grid has holes; interpolator would invent values'
print(f'complete grid: {len(NS)} fillings x {len(DS)} anisotropies, '
      f'delta {DS.min():g} to {DS.max():g}, {int(cov.values.min())} seeds/cell')

nf = np.linspace(NS.min(), NS.max(), 400); df = np.linspace(DS.min(), DS.max(), 400)
NG, DG = np.meshgrid(nf, df); pts = np.stack([DG.ravel(), NG.ravel()], -1)


def smooth_field(col):
    g = d.groupby(['n', 'delta'])[col].mean()
    M = g.unstack('n').values                      # (len(DS), len(NS)), the order
    F = RegularGridInterpolator((DS, NS), M)(pts).reshape(DG.shape)
    return gaussian_filter(F, BLUR) / NSITES


fig, ax = plt.subplots(1, 4, figsize=(23.0, 5.0), sharex=True, sharey=True,
                       facecolor='white')
for c, (key, sym) in enumerate(CHAN):
    a = ax[c]
    F = smooth_field(key)
    im = a.contourf(nf, df, F, levels=100, cmap='jet', extend='both')
    cb = fig.colorbar(im, ax=a, pad=0.02, shrink=0.92)
    cb.set_ticks(np.linspace(F.min(), F.max(), 5))
    cb.ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.4f'))
    cb.ax.tick_params(labelsize=10)
    a.set_title(f'({chr(97+c)}) {sym}', fontsize=16)
    a.set_xlabel(r'Filling  $n$')
    if c == 0: a.set_ylabel(r'Anisotropy  $\delta$')
    a.set_xticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    print(f'({chr(97+c)}) {key:9s} {F.min():+9.5f} to {F.max():+9.5f}')
plt.tight_layout()
plt.savefig(f'{D}/../figs/fig_4ch_uneq_L12.png', dpi=600,
            bbox_inches='tight', facecolor='white')
