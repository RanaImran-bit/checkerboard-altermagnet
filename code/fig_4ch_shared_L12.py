"""All four pairing channels over (delta, U), one shared colour scale.

L = 12, half filling, complete 7 x 7 grid, so griddata interpolates between
neighbours everywhere and no convex hull appears.

A SHARED scale is meaningful here because the four channels happen to lie within
one order of magnitude of each other (-0.156 to +0.503). That is not generally
true of these quantities, and it is why earlier versions of this figure gave
each panel its own scale. With one scale the panels ARE directly comparable by
colour: extended s is visibly the largest, on-site s the only negative one.

ON THE COLORBAR LABEL. It reads "pairing vertex", not Delta_tot. The field being
coloured is the pairing vertex at k = 0 for the channel named in each title.
Delta_tot is a different quantity, the momentum-space spin splitting, and
labelling this bar with it would mis-state what the colours are. If the colours
should instead BE Delta_tot, that is fig_dtot_bg_pair_contours.py, where the
identical Delta_tot field is drawn in every panel and the channel appears as
contour lines on top.
"""
import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt
from scipy.interpolate import griddata

mpl.rcParams.update({
    'font.family':'serif', 'font.serif':['DejaVu Serif'],
    'mathtext.fontset':'dejavuserif', 'axes.linewidth':1.2,
    'xtick.direction':'out', 'ytick.direction':'out',
    'xtick.top':False, 'ytick.right':False,
    'axes.labelsize':16, 'xtick.labelsize':12, 'ytick.labelsize':12})

D = '/Users/liujiaxin/Desktop/checkerboard-altermagnet/data'
d = pd.read_csv(f'{D}/fortran_L12_all49.csv')
d = d[np.isclose(d.n, 1.0) & (d.U >= 2)]
CH = [('son', r'on-site $s$'), ('sext', r'extended $s$'),
      ('d', r'$d_{x^2-y^2}$'), ('dxy', r'$d_{xy}$')]

VMIN = min(d[k].min() for k, _ in CH)
VMAX = max(d[k].max() for k, _ in CH)
LEV = np.linspace(VMIN, VMAX, 100)

gx, gy = np.meshgrid(np.linspace(d.delta.min(), d.delta.max(), 300),
                     np.linspace(d.U.min(), d.U.max(), 300))

fig, ax = plt.subplots(1, 4, figsize=(21.0, 4.7), facecolor='white',
                       sharex=True, sharey=True)
for c, (key, lab) in enumerate(CH):
    a = ax[c]
    Z = griddata((d.delta.values, d.U.values), d[key].values, (gx, gy), method='linear')
    assert not np.isnan(Z).any(), f'{key}: hull hole -- grid is not complete'
    im = a.contourf(gx, gy, Z, levels=LEV, cmap='jet', extend='both')
    a.set_title(f'({chr(97+c)}) {lab}', fontsize=17)
    a.set_xlabel(r'anisotropy  $\delta$')
    if c == 0: a.set_ylabel(r'$U/t$')
    print(f'({chr(97+c)}) {key:5s} {d[key].min():+.3f} to {d[key].max():+.3f}')

cb = fig.colorbar(im, ax=ax, pad=0.015, fraction=0.016)
cb.set_label(r'pairing vertex at $\mathbf{k}=0$', fontsize=16)
cb.ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.2f'))
cb.ax.tick_params(labelsize=11)
plt.savefig(f'{D}/../figs/fig_4ch_shared_L12.png', dpi=600, bbox_inches='tight',
            facecolor='white')
print(f'\nshared scale {VMIN:+.3f} to {VMAX:+.3f}')
