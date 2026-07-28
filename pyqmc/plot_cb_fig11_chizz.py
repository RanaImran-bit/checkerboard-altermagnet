#!/usr/bin/env python3
"""Fig 11: magnetic susceptibility chi_zz(q) over the Brillouin zone. Reads the
chizz_*.npy grids from checkerboard_chispin.py. Pass file/label pairs:
  python pyqmc/plot_cb_fig11_chizz.py \
     chizz_L6_n0.611_d0.0_s1.npy 'n=0.61, delta=0.0' \
     chizz_L6_n0.611_d0.4_s1.npy 'n=0.61, delta=0.4' \
     chizz_L6_n1.000_d0.4_s1.npy 'n=1.00, delta=0.4'
"""
import sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

pi = np.pi
args = sys.argv[1:]
pairs = [(args[i], args[i + 1]) for i in range(0, len(args) - 1, 2)]
if not pairs:
    print("usage: plot_cb_fig11_chizz.py f1.npy 'label1' [f2.npy 'label2' ...]"); sys.exit(1)

fig, axes = plt.subplots(1, len(pairs), figsize=(5 * len(pairs), 4.6),
                         facecolor='white', squeeze=False)
for ax, (fn, lab) in zip(axes[0], pairs):
    chi = np.fft.fftshift(np.load(fn))          # center (0,0); corners -> (pi,pi)
    im = ax.imshow(chi.T, origin='lower', cmap='magma',
                   extent=[-pi, pi, -pi, pi], aspect='equal')
    fig.colorbar(im, ax=ax, shrink=0.82)
    ax.set_xticks([-pi, 0, pi]); ax.set_yticks([-pi, 0, pi])
    ax.set_xticklabels([r'$-\pi$', '0', r'$\pi$'])
    ax.set_yticklabels([r'$-\pi$', '0', r'$\pi$'])
    ax.set_xlabel(r'$q_x$', fontsize=15)
    ax.set_title(r'$\chi_{zz}(\mathbf{q})$: ' + lab, fontsize=14)
axes[0][0].set_ylabel(r'$q_y$', fontsize=15)
plt.tight_layout()
plt.savefig("Fig11_chizz_q.png", dpi=300, bbox_inches='tight', facecolor='white')
print("saved Fig11_chizz_q.png")
