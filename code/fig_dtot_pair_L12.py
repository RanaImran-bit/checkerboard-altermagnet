"""Delta_tot and all four pairing channels over the SAME (delta, U) plane,
plus the Delta n(k) maps, from the complete L=12 half-filling grid.

What changed from fig_dtot_prl. That figure drew L=14, where the (U, delta)
grid held 19 of 35 cells, so griddata returned NaN outside the convex hull and
the panel carried a white diamond that read as missing physics. This grid is
49 of 49 with no holes, so the maps interpolate between neighbours everywhere
and no hull appears. It also reaches delta = 0.7 instead of 0.4.

Putting the pairing channels on the same axes as Delta_tot is the point of the
first figure: the magnetic and the pairing response can be compared cell by cell
rather than across two figures with different grids.

U = 0 is excluded from every panel. Delta_tot AND all four pairing channels are
identically zero there, so the row carries no information, and the U grid is
uneven below 2 (gap of 2, then 1, then 0.5), which would make the shading smooth
across a wide unsampled stretch.

Each panel keeps its own colour scale. The quantities differ by orders of
magnitude, so read the numbers, not the colours, across panels. Sign is on the
colorbar: the pairing channels take both signs, Delta_tot is positive by
construction.
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
z = np.load(f'{D}/fortran_L12_maps.npz', allow_pickle=True)
maps, meta = z['maps'], z['meta']                 # meta: L, n, U, delta

PAN = [('dtot_N', r'$\Delta_{tot}$'), ('son', r'on-site $s$'),
       ('sext', r'extended $s$'), ('d', r'$d_{x^2-y^2}$'), ('dxy', r'$d_{xy}$')]

# ---------- FIGURE 1: same (delta, U) plane for all five quantities ----------
fig, ax = plt.subplots(1, 5, figsize=(26.0, 4.9), facecolor='white',
                       sharex=True, sharey=True)
gx, gy = np.meshgrid(np.linspace(d.delta.min(), d.delta.max(), 300),
                     np.linspace(d.U.min(), d.U.max(), 300))
for c, (key, lab) in enumerate(PAN):
    a = ax[c]
    gi = griddata((d.delta.values, d.U.values), d[key].values, (gx, gy),
                  method='linear')
    assert not np.isnan(gi).any(), f'{key}: hull hole -- grid is not complete'
    im = a.contourf(gx, gy, gi, levels=100, cmap='jet', extend='both')
    cb = fig.colorbar(im, ax=a, pad=0.02)
    cb.ax.tick_params(labelsize=10)
    cb.ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3f'))
    a.plot(d.delta, d.U, 'o', ms=4, mfc='none', mec='k', mew=0.8)
    a.set_title(f'({chr(97+c)}) {lab}', fontsize=16)
    a.set_xlabel(r'anisotropy  $\delta$')
    if c == 0: a.set_ylabel(r'$U/t$')
    print(f'({chr(97+c)}) {key:7s} {d[key].min():+8.4f} to {d[key].max():+8.4f}')
plt.tight_layout()
plt.savefig(f'{D}/../figs/fig_dtot_pair_L12.png', dpi=600, bbox_inches='tight',
            facecolor='white')

# ---------- FIGURE 2: Delta n(k) over the Brillouin zone ----------
want = [(1.0, 4.0, 0.1), (1.0, 4.0, 0.4), (1.0, 4.0, 0.7)]
sel = []
for n_, U_, dl_ in want:
    hit = np.where(np.isclose(meta[:,1], n_) & np.isclose(meta[:,2], U_) &
                   np.isclose(meta[:,3], dl_))[0]
    if len(hit): sel.append((hit[0], dl_))

fig2, ax2 = plt.subplots(1, len(sel), figsize=(5.0*len(sel), 4.8), facecolor='white')
# one shared symmetric scale, so a weak panel reads as weak instead of being
# stretched to fill its own range
vmax = max(np.abs(maps[i][:,2]).max() for i, _ in sel)
for j, (i, dl_) in enumerate(sel):
    a = ax2[j] if len(sel) > 1 else ax2
    kx, ky, dn = maps[i][:,0], maps[i][:,1], maps[i][:,2]
    # the file writes pi as 3.141592741 (float32), 8.7e-8 LARGER than float64
    # pi, so a bare > np.pi sends the +pi points to -pi and leaves a white strip
    kxf = np.where(kx > np.pi + 1e-6, kx - 2*np.pi, kx)
    kyf = np.where(ky > np.pi + 1e-6, ky - 2*np.pi, ky)
    ggx, ggy = np.meshgrid(np.linspace(-np.pi, np.pi, 260),
                           np.linspace(-np.pi, np.pi, 260))
    gz = griddata((kxf, kyf), dn, (ggx, ggy), method='linear')
    im = a.contourf(ggx, ggy, gz, levels=80, cmap='RdBu_r', vmin=-vmax, vmax=vmax)
    a.set_title(rf'$\delta = {dl_:g}$', fontsize=17)
    a.set_xlabel(r'$k_x$'); a.set_aspect('equal')
    a.set_xticks([-np.pi, 0, np.pi]); a.set_xticklabels([r'$-\pi$','0',r'$\pi$'])
    a.set_yticks([-np.pi, 0, np.pi])
    a.set_yticklabels([r'$-\pi$','0',r'$\pi$'] if j == 0 else [])
    if j == 0: a.set_ylabel(r'$k_y$')
    print(f'  dn(k) delta={dl_:g}: min {dn.min():+.4f}  max {dn.max():+.4f}')
cb = fig2.colorbar(im, ax=ax2, pad=0.02, fraction=0.03)
cb.set_label(r'$\Delta n(\mathbf{k}) = n_\uparrow(\mathbf{k}) - n_\downarrow(\mathbf{k})$',
             fontsize=15)
fig2.suptitle(r'$L=12$, half filling, $U=4$', fontsize=17, y=1.02)
plt.savefig(f'{D}/../figs/fig_dnk_L12.png', dpi=600, bbox_inches='tight',
            facecolor='white')
