"""Fig 7: d-wave pairing phase diagram over the (n, delta) plane. RUN ON 251.
(a) N^Vertex_dx2-y2 (largest at small delta), (b) N^Vertex_dxy (largest at large delta),
(c) difference N_dxy - N_dx2-y2 with the dashed CROSSOVER line (where the two are equal).
Whole-plane version of Fig 6; analog of paper-2 Fig 2(b,c) "opposite regions". Emergent,
spin-independent checkerboard => the crossover is GEOMETRY-driven (delta), not an imposed
spin-dependent term. Vertex_* summed over k (connected). Folders: tA = -delta."""
import os, re
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

BASE = "/home/phd25imran/L14_archive"
U, L = 4.0, 14
PAT = re.compile(r'L(\d+)n([\d.]+)u([\d.]+)tA(-?[\d.]+)tt([\d.-]+)N(\d+)')

def vsum(fp):
    d = pd.read_csv(fp, sep=r'\s+', header=None, skiprows=1)   # kx ky value error
    d = d[(d[0] != 3.141592741) | (d[1] != -3.141592741)]
    return d[2].sum()

rows = []
for f in os.listdir(BASE):
    m = PAT.match(f)
    if not m:
        continue
    if int(m.group(1)) != L or abs(float(m.group(3)) - U) > 1e-6:
        continue
    n, tA = float(m.group(2)), float(m.group(4))
    kdir = os.path.join(BASE, f, 'dir-kVals')
    fx = os.path.join(kdir, 'Vertex_dwave.dat')       # dx2-y2
    fy = os.path.join(kdir, 'Vertex_dd12wave.dat')    # dxy
    if not (os.path.isfile(fx) and os.path.isfile(fy)):
        continue
    rows.append((n, -tA, vsum(fx), vsum(fy)))         # delta = -tA

df = pd.DataFrame(rows, columns=['n', 'delta', 'dx2', 'dxy'])
df['diff'] = df['dxy'] - df['dx2']
print(f"{len(df)} points, delta {sorted(df.delta.unique())}")

fig, axes = plt.subplots(1, 3, figsize=(18, 5.2), facecolor='white')
titles = [r'(a) $N^{\mathrm{Vertex}}_{d_{x^2-y^2}}$',
          r'(b) $N^{\mathrm{Vertex}}_{d_{xy}}$',
          r'(c) $N_{d_{xy}} - N_{d_{x^2-y^2}}$']
cols   = ['dx2', 'dxy', 'diff']

for ax, col, title in zip(axes, cols, titles):
    if col == 'diff':
        vmax = np.abs(df['diff']).max()
        sc = ax.scatter(df.n, df.delta, c=df[col], cmap='coolwarm', vmin=-vmax, vmax=vmax,
                        s=220, edgecolors='k', linewidths=0.6, zorder=3)
        ni = np.linspace(df.n.min(), df.n.max(), 120)
        di = np.linspace(df.delta.min(), df.delta.max(), 120)
        NI, DI = np.meshgrid(ni, di)
        G = griddata((df.n, df.delta), df['diff'], (NI, DI), method='linear')
        ax.contour(NI, DI, G, levels=[0], colors='k', linewidths=2.0, linestyles='--')
    else:
        sc = ax.scatter(df.n, df.delta, c=df[col], cmap='viridis',
                        s=220, edgecolors='k', linewidths=0.6, zorder=3)
    cb = fig.colorbar(sc, ax=ax, shrink=0.85, pad=0.02)
    cb.ax.tick_params(labelsize=12)
    ax.set_xlabel(r'Filling $n$', fontsize=17)
    if ax is axes[0]:
        ax.set_ylabel(r'Anisotropy $\delta$', fontsize=17)
    ax.set_title(title, fontsize=16)
    ax.tick_params(axis='both', labelsize=13)

plt.tight_layout()
plt.savefig("Fig7_pairing_phase_diagram.png", dpi=600, bbox_inches='tight', facecolor='white')
print("saved Fig7_pairing_phase_diagram.png")
