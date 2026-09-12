"""Delta_tot over (delta, U) at L = 8, 10 and 12, half filling.

All three are Fortran runs at the production settings (Delta_tau = 0.01,
beta = 32, N_walkers = 1000, E_T = -50), identical to the published PRB. The
beta = 3 problem applied only to the Python chi drivers, never to Delta_tot.

NORMALISATION. The paper's Eq. (3) is a bare sum, sum_k |n_up - n_dn|, with no
1/L^2. That is what dtot_raw holds. But the raw sum grows with the number of
k-points, so comparing it across L would show a size effect that is pure
counting. Every panel here uses dtot_N = dtot_raw / L^2, per site, which is the
only form comparable across L. Quote Eq. (3) for absolute values and define the
per-site quantity separately in the manuscript.

DUPLICATE k. n_up.dat stores (L+1)^2 rows: the +pi zone edge repeats -pi, so a
naive sum counts the boundary up to 4 times at the corners and inflates
Delta_tot by 0.4% at L=8 rising to more with L -- biasing exactly the trend
being measured here. The extraction snaps |k| within 1e-4 of pi onto -pi, which
recovers exactly L^2 unique momenta at every size.

One shared colour scale across the three panels, otherwise each autoscales and
the finite-size trend is invisible.
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
d = pd.read_csv(f'{D}/fortran_L8_L10_L12_dedup.csv')
d = d[np.isclose(d.n, 1.0) & (d.U >= 2)]
LS = sorted(d.L.unique())
VMAX = d.dtot_N.max()
LEV = np.linspace(0, VMAX, 100)

fig, ax = plt.subplots(1, len(LS), figsize=(6.2*len(LS), 5.2), facecolor='white',
                       sharex=True, sharey=True)
for i, L in enumerate(LS):
    a = ax[i]; s = d[d.L == L]
    gx, gy = np.meshgrid(np.linspace(s.delta.min(), s.delta.max(), 300),
                         np.linspace(s.U.min(), s.U.max(), 300))
    Z = reggrid(s, 'delta', 'U', 'dtot_N', gx, gy)
    assert not np.isnan(Z).any(), f'L={L}: hull hole -- grid is not complete'
    im = a.contourf(gx, gy, Z, levels=LEV, cmap='jet', extend='both')
    a.set_title(rf'$L={L}$', fontsize=17)
    a.set_xlabel(r'anisotropy  $\delta$')
    if i == 0: a.set_ylabel(r'$U/t$')
    peak = s.loc[s.dtot_N.idxmax()]
    print(f'L={L:2d}: {len(s)} cells, dtot_N {s.dtot_N.min():.4f}..{s.dtot_N.max():.4f}'
          f'   peak at U={peak.U:g}, delta={peak.delta:g}')
    print('     ridge (delta of max at each U): ' +
          '  '.join(f'U{U:g}:{s[s.U==U].loc[s[s.U==U].dtot_N.idxmax(),"delta"]:.1f}'
                    for U in sorted(s.U.unique())))
cb = fig.colorbar(im, ax=ax, pad=0.015, fraction=0.018)
cb.set_label(r'$\Delta_{tot}/N$', fontsize=17)
cb.ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3f'))
cb.ax.tick_params(labelsize=11)
plt.savefig(f'{D}/../figs/fig_dtot_L8_L10_L12.png', dpi=600, bbox_inches='tight',
            facecolor='white')
print(f'\nshared colour scale 0 to {VMAX:.4f}')
