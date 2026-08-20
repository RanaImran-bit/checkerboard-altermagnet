%matplotlib inline
import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt
from scipy.interpolate import griddata

mpl.rcParams.update({
    'font.family':'serif', 'font.serif':['DejaVu Serif'],
    'mathtext.fontset':'dejavuserif', 'axes.linewidth':1.2,
    'xtick.direction':'out', 'ytick.direction':'out',
    'xtick.top':False, 'ytick.right':False,
    'axes.labelsize':15, 'xtick.labelsize':15, 'ytick.labelsize':15})

D = '/home/phd25imran/analysis'
CH = [('son', r'on-site $s$'), ('sext', r'extended $s$'),
      ('d', r'$d_{x^2-y^2}$'), ('dxy', r'$d_{xy}$')]
d = pd.read_csv(f'{D}/fortran_L8_L10_L12_dedup.csv')
d = d[np.isclose(d.n, 1.0) & (d.U >= 2)]
LS = sorted(d.L.unique())
LEV = np.linspace(0, d.dtot_N.max(), 100)

fig, ax = plt.subplots(len(LS), len(CH), figsize=(5.5*len(CH), 4.4*len(LS)),
                       facecolor='white', sharex=True, sharey=True)
gx, gy = np.meshgrid(np.linspace(d.delta.min(), d.delta.max(), 300),
                     np.linspace(d.U.min(), d.U.max(), 300))
for i, L in enumerate(LS):
    s = d[d.L == L]
    BG = griddata((s.delta.values, s.U.values), s.dtot_N.values, (gx, gy), method='linear')
    assert not np.isnan(BG).any(), f'L={L}: hull hole'
    for j, (key, lab) in enumerate(CH):
        a = ax[i, j]
        im = a.contourf(gx, gy, BG, levels=LEV, cmap='jet', extend='both')
        F = griddata((s.delta.values, s.U.values), s[key].values, (gx, gy), method='linear')
        cs = a.contour(gx, gy, F, levels=6, colors='k', linewidths=1.2)
        a.clabel(cs, inline=True, fontsize=8, fmt='%.2f')
        if i == 0: a.set_title(lab, fontsize=20, pad=10)
        if i == len(LS)-1: a.set_xlabel(r'anisotropy  $\delta$', fontsize=22)
        if j == 0: a.set_ylabel(rf'$L={L}$' '\n' r'$U/t$', fontsize=18)
cb = fig.colorbar(im, ax=ax, pad=0.012, fraction=0.014)
cb.set_label(r'$\Delta_{tot}/N$', fontsize=20)
cb.ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3f'))
cb.ax.tick_params(labelsize=15)
plt.savefig('/home/phd25imran/analysis/fig_dtot_colour_ch_contours_3L.png', dpi=600,
            bbox_inches='tight', facecolor='white')
plt.show()
print(f'colour = Delta_tot/N, one scale for all panels: 0 .. {d.dtot_N.max():.4f}')
for L in LS:
    s = d[d.L == L]
    print(f'  L={L:2d}: contour ranges ' + '  '.join(f'{k}:{s[k].min():+.2f}..{s[k].max():+.2f}' for k,_ in CH))
