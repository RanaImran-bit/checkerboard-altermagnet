"""All four pairing channels, equal time against unequal time, over (n, delta).

Extends the reference two-channel figure to all four channels. The delta range
is UNCHANGED at 0 to 0.4, because that is where the equal-time data ends. See
the note at the bottom of this file for what reaching 0.7 would take.

Units. C(0) is a correlation; chi = integral of C(tau) d tau carries an extra
factor of imaginary time and is not comparable to C(0) as it stands. Dividing
the equal-time row by N and the unequal-time row by N*tau_max puts both in the
same units, which is the only way the two rows can honestly sit in one figure.

Each panel keeps its own colour scale. The four channels differ by more than an
order of magnitude, so colours must NOT be read across panels, only the numbers.
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
NSITES, TAU_MAX = 100, 16*0.05          # L = 10; MUST match BP, DT in the driver
BLUR = 4.0
CHAN = [('son', r'on-site $s$'), ('sext', r'extended $s$'),
        ('d', r'$d_{x^2-y^2}$'), ('dxy', r'$d_{xy}$')]
SECT = [('eq', NSITES, 'equal time'), ('chi', NSITES*TAU_MAX, 'unequal time')]

d  = pd.read_csv(f'{D}/eqtime_L10_U4_all.csv')
NS = np.sort(d.n.unique()); DS = np.sort(d.delta.unique())
nf = np.linspace(NS.min(), NS.max(), 400); df = np.linspace(DS.min(), DS.max(), 400)
NG, DG = np.meshgrid(nf, df); pts = np.stack([DG.ravel(), NG.ravel()], -1)


def smooth_field(col, scale):
    g = d.groupby(['n', 'delta'])[col].mean()
    assert len(g) == len(NS)*len(DS), 'grid has holes'
    M = g.unstack('n').values                     # (len(DS), len(NS)), the order
    F = RegularGridInterpolator((DS, NS), M)(pts).reshape(DG.shape)
    return gaussian_filter(F, BLUR) / scale


fig, ax = plt.subplots(2, 4, figsize=(23.0, 9.4), sharex=True, sharey=True,
                       facecolor='white')
for r, (pre, sc, sect) in enumerate(SECT):
    for c, (key, sym) in enumerate(CHAN):
        a = ax[r, c]
        F = smooth_field(f'{pre}_{key}_vertex', sc)
        im = a.contourf(nf, df, F, levels=100, cmap='jet', extend='both')
        cb = fig.colorbar(im, ax=a, pad=0.02, shrink=0.92)
        cb.set_ticks(np.linspace(F.min(), F.max(), 5))
        cb.ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3f'))
        cb.ax.tick_params(labelsize=10)
        a.set_title(f'({chr(97 + 4*r + c)}) {sym}, {sect}', fontsize=15)
        if r == 1: a.set_xlabel(r'Filling  $n$')
        if c == 0: a.set_ylabel(r'Anisotropy  $\delta$')
        a.set_xticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
        print(f'({chr(97 + 4*r + c)}) {key:5s} {sect:13s} '
              f'{F.min():+9.4f} to {F.max():+9.4f}')
plt.tight_layout()
plt.savefig(f'{D}/../figs/fig_4ch_eq_vs_uneq.png', dpi=600,
            bbox_inches='tight', facecolor='white')

# ---------------------------------------------------------------------------
# WHY delta STOPS AT 0.4.
# eqtime_L10_U4_all.csv is the ONLY dataset carrying an equal-time measurement,
# and it runs delta = 0, 0.1, 0.2, 0.3, 0.4 at L = 10, U = 4.
# pr_all.csv looks like a candidate because it reaches delta = 0.7 with all four
# channels, but it is NOT equal time: summing its vertex over R reproduces
# pairing_master's chi exactly (half filling, U = 4: +3.98, +4.05, +2.68, +2.37,
# +1.78, +1.38, +1.19, +1.04 for d). It is the distance-resolved decomposition
# of the same tau-integrated quantity.
# The unequal-time row alone can reach 0.7 today, from dense_L12_U4_all.csv.
# A matched pair of rows to 0.7 needs new equal-time runs at delta = 0.5, 0.6, 0.7.
# ---------------------------------------------------------------------------
