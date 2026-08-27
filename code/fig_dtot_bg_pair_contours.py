"""Delta_tot as the colour field, pairing channels as contours on the same axes.

Every panel shows the IDENTICAL Delta_tot map over (delta, U) at L = 12, half
filling, complete 7 x 7 grid. Overlaid on each is one pairing channel drawn as
labelled contour lines. So the colour answers "where is the altermagnetic
splitting large" and the lines answer "where is this channel strong", on one
set of axes, cell by cell.

That is the comparison the two-figure version could not make: the ridge in
Delta_tot and the maximum of dxy sit in visibly different places.

NORMALISATION. Delta_tot here is the raw momentum sum divided by the number of
sites, i.e. per site. This is OUR convention, chosen so different L are
comparable. It is NOT copied from the reference paper's equation, and earlier
versions of these scripts wrongly attributed it there. Quote the source
equation verbatim in the manuscript and define this one separately.

    Delta_tot / N  =  (1/L^2) sum_k |n_up(k) - n_dn(k)|

U = 0 is excluded: Delta_tot and all four channels vanish identically there, and
the U grid is uneven below 2, so shading would smooth across an unsampled gap.
"""
import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt
from scipy.interpolate import griddata
from gridinterp import reggrid

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

gx, gy = np.meshgrid(np.linspace(d.delta.min(), d.delta.max(), 300),
                     np.linspace(d.U.min(), d.U.max(), 300))
BG = reggrid(d, 'delta', 'U', 'dtot_N', gx, gy)
assert not np.isnan(BG).any(), 'hull hole -- grid is not complete'

fig, ax = plt.subplots(1, 4, figsize=(22.0, 4.9), facecolor='white',
                       sharex=True, sharey=True)
for c, (key, lab) in enumerate(CH):
    a = ax[c]
    # identical background in every panel: colour is ALWAYS Delta_tot
    im = a.contourf(gx, gy, BG, levels=100, cmap='jet',
                    vmin=d.dtot_N.min(), vmax=d.dtot_N.max(), extend='both')
    F = reggrid(d, 'delta', 'U', key, gx, gy)
    cs = a.contour(gx, gy, F, levels=7, colors='k', linewidths=1.3)
    a.clabel(cs, inline=True, fontsize=9, fmt='%.2f')
    a.plot(d.delta, d.U, 'o', ms=4, mfc='none', mec='k', mew=0.8)
    a.set_title(f'({chr(97+c)}) contours: {lab}', fontsize=16)
    a.set_xlabel(r'anisotropy  $\delta$')
    if c == 0: a.set_ylabel(r'$U/t$')
    print(f'({chr(97+c)}) {key:5s} contour range {d[key].min():+.3f} to {d[key].max():+.3f}')

cb = fig.colorbar(im, ax=ax, pad=0.015, fraction=0.016)
cb.set_label(r'$\Delta_{tot}$', fontsize=16)
cb.ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3f'))
cb.ax.tick_params(labelsize=11)
plt.savefig(f'{D}/../figs/fig_dtot_bg_pair_contours.png', dpi=600,
            bbox_inches='tight', facecolor='white')
print(f'\ncolour field identical in all panels: Delta_tot/N '
      f'{d.dtot_N.min():.4f} to {d.dtot_N.max():.4f}')
