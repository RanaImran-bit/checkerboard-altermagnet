"""Fig 8: k-space pairing vertex maps, momentum-resolved d-wave crossover. RUN ON 251.
rows = delta (0.1 dx2-y2-dominant, 0.4 dxy-dominant); cols = Vertex_dx2-y2, Vertex_dxy. n=0.622.
As delta grows, dx2-y2 weakens (center even goes negative) and dxy grows into a coherent X
peaked on the zone diagonals at (+-pi/2, +-pi/2) -- the SAME momenta where the altermagnetic
spin splitting m=4*delta*sin kx sin ky is maximal. Pairing and altermagnetism share the dxy
structure (the thesis). Analog of paper-2 Fig 5. Folders: tA = -delta, tt = 0.3."""
import os, re
import numpy as np, pandas as pd
import matplotlib.pyplot as plt

BASE = "/home/phd25imran/L14_archive"
U, L, N_TARGET, TT = 4.0, 14, 0.622, 0.3
PAT = re.compile(r'L(\d+)n([\d.]+)u([\d.]+)tA(-?[\d.]+)tt([\d.-]+)N(\d+)')

def find_folder(n, delta, tol=0.02):
    target_tA = -delta
    best, bs = None, 1e9
    for f in os.listdir(BASE):
        m = PAT.match(f)
        if not m:
            continue
        if int(m.group(1)) != L or abs(float(m.group(3)) - U) > 1e-6:
            continue
        nv, tAv, tt = float(m.group(2)), float(m.group(4)), float(m.group(5))
        if abs(nv-n) > tol or abs(tAv-target_tA) > tol or abs(tt-TT) > tol:
            continue
        s = abs(nv-n) + abs(tAv-target_tA)
        if s < bs and os.path.isdir(os.path.join(BASE, f, 'dir-kVals')):
            bs, best = s, f
    return best

def kgrid(fp):
    d = pd.read_csv(fp, sep=r'\s+', header=None, skiprows=1)   # kx ky value error
    xu, yu = np.unique(d[0]), np.unique(d[1])
    z = np.full((len(yu), len(xu)), np.nan)
    xi = {v:i for i,v in enumerate(xu)}; yi = {v:i for i,v in enumerate(yu)}
    for _, r in d.iterrows():
        z[yi[r[1]], xi[r[0]]] = r[2]
    return xu, yu, z

DELTAS = [0.10, 0.40]
COLS = [('Vertex_dwave',    r'$N^{\mathrm{Vertex}}_{d_{x^2-y^2}}$'),
        ('Vertex_dd12wave', r'$N^{\mathrm{Vertex}}_{d_{xy}}$')]

fig, axes = plt.subplots(2, 2, figsize=(11, 10), facecolor='white')
pi = np.pi
p = 0
for r, delta in enumerate(DELTAS):
    folder = find_folder(N_TARGET, delta)
    for c, (fname, label) in enumerate(COLS):
        ax = axes[r][c]
        if folder is None:
            ax.text(0.5, 0.5, 'missing', transform=ax.transAxes, ha='center'); continue
        fp = os.path.join(BASE, folder, 'dir-kVals', fname + '.dat')
        xu, yu, z = kgrid(fp)
        vmax = np.nanmax(np.abs(z))
        cf = ax.contourf(xu, yu, z, levels=60, cmap='coolwarm', vmin=-vmax, vmax=vmax)
        cb = fig.colorbar(cf, ax=ax, shrink=0.85, pad=0.03); cb.ax.tick_params(labelsize=12)
        ax.set_aspect('equal')
        ax.set_xticks([-pi, 0, pi]); ax.set_yticks([-pi, 0, pi])
        ax.set_xticklabels([r'$-\pi$', '0', r'$\pi$'], fontsize=15)
        ax.set_yticklabels([r'$-\pi$', '0', r'$\pi$'], fontsize=15)
        if r == 1: ax.set_xlabel(r'$k_x$', fontsize=17)
        if c == 0: ax.set_ylabel(r'$k_y$', fontsize=17)
        ax.set_title(rf"({chr(97+p)}) {label},  $\delta={delta}$", fontsize=15)
        p += 1

plt.tight_layout()
plt.savefig("Fig8_pairing_kspace.png", dpi=600, bbox_inches='tight', facecolor='white')
print("saved Fig8_pairing_kspace.png")
