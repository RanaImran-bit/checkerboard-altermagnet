"""All four channels, equal time against unequal time, over (n, delta).

WHY THE TWO ROWS ALWAYS SHARE A delta RANGE. checkerboard_eqtime.py writes
eq_*  (the tau = 0 value) and chi_*  (the integrated curve) from the SAME run,
so every cell that has one has the other. There is no data state in which the
unequal-time row reaches further in delta than the equal-time row. What varies
is which (n, delta) cells were run at all.

Coverage as of this writing:
    delta <= 0.4 : all 11 fillings
    delta 0.5-0.7: only n <= 0.82  (nup 25 to 41)
The near-half-filling cells at high delta were still running when this was
written, which is precisely the region the dxy result lives in.

So the script does not hard-code a delta range. It merges every eqtime file it
finds and picks the largest COMPLETE rectangle of (n, delta), because
RegularGridInterpolator requires a full grid and silently produces nonsense if
handed one with holes. Re-run it after more cells land and the panel simply
grows -- no edit needed.

Units. C(0) is a correlation; chi carries an extra factor of imaginary time.
Dividing the equal-time row by N and the unequal-time row by N*tau_max puts both
in the same units, which is what lets the two rows sit in one figure at all.

Each panel keeps its own colour scale: the channels differ by more than an order
of magnitude, so colours must NOT be read across panels, only the numbers.
"""
import glob, os
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

# ---- merge every eqtime file present, newest wins on duplicate cells ----
files = sorted(glob.glob(f'{D}/eqtime_L10_U4*.csv'))
d = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
d = d.drop_duplicates(subset=['nup', 'delta', 'seed'], keep='last')
print('merged:', ', '.join(os.path.basename(f) for f in files))

# ---- largest complete rectangle: keep every delta, drop fillings with holes ----
cov = d.groupby(['n', 'delta']).size().unstack()
full_n = cov.dropna().index.values                    # fillings present at ALL delta
if len(full_n) < 3:                                   # too few, fall back on delta
    cov = cov.dropna(axis=1)
    full_n = cov.dropna().index.values
NS = np.sort(full_n); DS = np.sort(cov.columns.values)
d = d[d.n.isin(NS) & d.delta.isin(DS)]
print(f'complete grid: {len(NS)} fillings x {len(DS)} anisotropies '
      f'(delta {DS.min():g} to {DS.max():g})')
print(f'  fillings: {", ".join(f"{x:.2f}" for x in NS)}')

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
plt.savefig(f'{D}/../figs/fig_4ch_eq_vs_uneq_auto.png', dpi=600,
            bbox_inches='tight', facecolor='white')
