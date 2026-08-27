"""Appendix Fig A: uncorrelated (bubble) vs connected (vertex) d-wave pairing at n=0.622. RUN ON 251.
Answers the kinematic objection. (a) Unpair (uncorrelated/bubble) part: dx2-y2 & dxy vary WEAKLY with
delta, cross near delta~0.32. (b) Vertex (connected/interaction) part: dx2-y2 COLLAPSES ~75% and dxy more
than DOUBLES, crossover at delta~0.22. => geometry sets the DIRECTION, correlations DRIVE & SHARPEN it.
Full = Unpair + Vertex. Both files present per channel. Folders: tA = -delta, tt = 0.3."""
import os, re
import numpy as np, pandas as pd
import matplotlib.pyplot as plt

BASE = "/home/phd25imran/L14_archive"
N_TARGET, TT, U, L = 0.622, 0.3, 4.0, 14
PAT = re.compile(r'L(\d+)n([\d.]+)u([\d.]+)tA(-?[\d.]+)tt([\d.-]+)N(\d+)')

def vsum(fp):
    d = pd.read_csv(fp, sep=r'\s+', header=None, skiprows=1)   # kx ky value error
    d = d[(d[0] != 3.141592741) | (d[1] != -3.141592741)]
    return d[2].sum()

def crossing(ds, y1, y2):
    diff = np.array(y2) - np.array(y1)
    for i in range(len(diff)-1):
        if diff[i] <= 0 <= diff[i+1] and diff[i+1] != diff[i]:
            return ds[i] - diff[i]*(ds[i+1]-ds[i])/(diff[i+1]-diff[i])
    return None

data = {}   # delta -> {(prefix, channel): value}
for f in os.listdir(BASE):
    m = PAT.match(f)
    if not m:
        continue
    if int(m.group(1)) != L or abs(float(m.group(3))-U) > 1e-6:
        continue
    n, tA, tt = float(m.group(2)), float(m.group(4)), float(m.group(5))
    if abs(n-N_TARGET) > 0.02 or abs(tt-TT) > 0.02:
        continue
    k = os.path.join(BASE, f, 'dir-kVals')
    d = -tA; data[d] = {}
    for pfx in ('Unpair', 'Vertex'):
        for ch in ('dwave', 'dd12wave'):
            fp = os.path.join(k, f'{pfx}_{ch}.dat')
            if os.path.isfile(fp):
                data[d][(pfx, ch)] = vsum(fp)

ds = sorted(data)
def series(pfx, ch): return [data[d][(pfx, ch)] for d in ds]

fig, (axb, axv) = plt.subplots(1, 2, figsize=(13, 5.4), facecolor='white')
for ax, pfx, tag in [(axb, 'Unpair', '(a) uncorrelated (bubble)'),
                     (axv, 'Vertex', '(b) connected vertex')]:
    x2 = series(pfx, 'dwave'); xy = series(pfx, 'dd12wave')
    ax.plot(ds, x2, 'o-', ms=10, lw=2.2, color='#1f77b4', markerfacecolor='#1f77b4', label=r'$d_{x^2-y^2}$')
    ax.plot(ds, xy, '^-', ms=10, lw=2.2, color='#2ca02c', markerfacecolor='#2ca02c', label=r'$d_{xy}$')
    c = crossing(ds, x2, xy)
    if c is not None:
        ax.axvline(c, color='gray', ls=':', lw=1.5)
        ax.text(c, ax.get_ylim()[0], rf'  $\delta_c\approx{c:.2f}$', va='bottom', fontsize=12, color='dimgray')
    ax.set_xlabel(r'Anisotropy $\delta$', fontsize=18)
    ax.set_ylabel(r'$N_\zeta$', fontsize=18)
    ax.set_title(tag, fontsize=16)
    ax.tick_params(axis='both', labelsize=13)
    ax.legend(fontsize=15, frameon=False, loc='center left')

plt.tight_layout()
plt.savefig("FigA_bare_vs_vertex.png", dpi=600, bbox_inches='tight', facecolor='white')
print("saved FigA_bare_vs_vertex.png")
