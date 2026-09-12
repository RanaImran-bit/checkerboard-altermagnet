"""n_up(k) over the Brillouin zone, from the Fortran output.

Reads dir-kVals/n_up.dat directly. Companion script plots n_dn(k) on the SAME
colour scale, so the two can be laid side by side and compared honestly.

Note before comparing by eye: at L=14, U=4, delta=0.2 the spin difference is
mean |n_up - n_dn| = 0.050 while n(k) itself is of order 0.5 to 1.0. The
asymmetry is therefore ~7% of the value, and the two maps look nearly identical.
The splitting is real but only becomes visible in the DIFFERENCE.
"""
import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt
from scipy.interpolate import griddata

mpl.rcParams.update({
    'font.family':'serif', 'font.serif':['DejaVu Serif'],
    'mathtext.fontset':'dejavuserif', 'axes.linewidth':1.2,
    'xtick.direction':'out', 'ytick.direction':'out',
    'xtick.top':False, 'ytick.right':False,
    'axes.labelsize':17, 'xtick.labelsize':13, 'ytick.labelsize':13})

BASE  = '/home/phd25imran/Checkerboard_Model'
SPIN  = 'up'                       # <-- 'up' or 'dn'; the only line to change
ARROW = {'up': r'\uparrow', 'dn': r'\downarrow'}[SPIN]
RUNS  = [('L14n1.000u4.0tA0.0tt0.3N98',   0.0),
         ('L14n1.000u4.0tA-0.2tt0.3N98',  0.2),
         ('L14n1.000u4.0tA-0.4tt0.3N98',  0.4)]

# One shared colour scale across all panels and BOTH spins, otherwise each panel
# autoscales and a 7% difference is stretched to look like a large one.
VMIN, VMAX = 0.0, 1.0

fig, ax = plt.subplots(1, len(RUNS), figsize=(5.0*len(RUNS), 4.8), facecolor='white')
for j, (run, dl) in enumerate(RUNS):
    a = ax[j] if len(RUNS) > 1 else ax
    f = pd.read_csv(f'{BASE}/{run}/dir-kVals/n_{SPIN}.dat', sep=r'\s+',
                    skiprows=1, header=None, names=['kx','ky','v','e'])
    # The file already spans [-pi, pi], so no folding is needed. An earlier
    # version folded on kx > np.pi, which was wrong twice over: the convention
    # is not [0, 2pi), AND the file writes pi as 3.141592741 (float32), which is
    # 8.7e-8 LARGER than numpy's float64 pi. That sent the +pi points to -pi and
    # left the strip between 2.693 and pi with no data -- the white edge.
    # The tolerance keeps this safe if a file ever does use [0, 2pi).
    kx = np.where(f.kx > np.pi + 1e-6, f.kx - 2*np.pi, f.kx)
    ky = np.where(f.ky > np.pi + 1e-6, f.ky - 2*np.pi, f.ky)
    gx, gy = np.meshgrid(np.linspace(-np.pi, np.pi, 260),
                         np.linspace(-np.pi, np.pi, 260))
    gz = griddata((kx, ky), f.v.values, (gx, gy), method='linear')
    im = a.contourf(gx, gy, gz, levels=np.linspace(VMIN, VMAX, 80),
                    cmap='jet', extend='both')
    a.set_title(rf'$\delta = {dl:g}$', fontsize=17)
    a.set_xlabel(r'$k_x$'); a.set_aspect('equal')
    a.set_xticks([-np.pi, 0, np.pi]); a.set_xticklabels([r'$-\pi$','0',r'$\pi$'])
    a.set_yticks([-np.pi, 0, np.pi])
    a.set_yticklabels([r'$-\pi$','0',r'$\pi$'] if j == 0 else [])
    if j == 0: a.set_ylabel(r'$k_y$')
    print(f'delta={dl:g}  n_{SPIN}(k): min {f.v.min():.4f}  max {f.v.max():.4f}  '
          f'mean {f.v.mean():.4f}')
cb = fig.colorbar(im, ax=ax, pad=0.02, fraction=0.03)
cb.set_label(rf'$n_{{{ARROW}}}(\mathbf{{k}})$', fontsize=16)
fig.suptitle(rf'$n_{{{ARROW}}}(\mathbf{{k}})$   $L=14$, half filling, $U=4$',
             fontsize=17, y=1.02)
plt.savefig(f'/home/phd25imran/analysis/fig_n{SPIN}_k.png', dpi=600,
            bbox_inches='tight', facecolor='white')
plt.show()
