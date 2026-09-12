"""Fig 6: pairing vertex N_zeta vs anisotropy delta. RUN ON 251 (needs Vertex_*.dat in dir-kVals).
(a) all channels at fixed doped n=0.622 (full context: on-site s repulsive, ext-s largest).
(b) zoom on the d-wave crossover for BOTH fillings n=0.622 (solid) and n=0.582 (dashed);
    d_{x2-y2} falls and d_{xy} rises with delta, crossing near delta~0.11-0.22 (shaded band).
HONEST CLAIM: among the UNCONVENTIONAL (sign-changing) channels, delta selects d_{xy} over d_{x2-y2}.
Do NOT claim d_{xy} is the globally dominant pairing (ext-s is larger; it is the conventional/AFM-tied
channel and not the SC instability). Vertex_* is the CONNECTED vertex -> sum over k directly (no Unpair
subtraction). Channel map Fortran-verified from mc2duph.f90. Folders: tA = -delta, tt = 0.3."""
import os, re
import numpy as np, pandas as pd
import matplotlib.pyplot as plt

BASE = "/home/phd25imran/L14_archive"
TT, U, L = 0.3, 4.0, 14
PAT = re.compile(r'L(\d+)n([\d.]+)u([\d.]+)tA(-?[\d.]+)tt([\d.-]+)N(\d+)')

# symmetry -> (file, label, marker, color)   [Fortran-verified channel map]
CHANNELS = [
    ("dwave",      r'$d_{x^2-y^2}$',            'o', '#1f77b4'),
    ("dd12wave",   r'$d_{xy}$',                 '^', '#2ca02c'),
    ("sowave",     r'$s$',                      '*', '#d62728'),
    ("swave",      r'ext. $s$',                 'D', '#9467bd'),
    ("puupywave",  r'$p_{\uparrow\uparrow}$',   's', '#ff7f0e'),
    ("pudpywave",  r'$p_{\uparrow\downarrow}$', 'v', '#8c564b'),
]

def vsum(fp):
    d = pd.read_csv(fp, sep=r'\s+', header=None, skiprows=1)   # kx ky value error
    d = d[(d[0] != 3.141592741) | (d[1] != -3.141592741)]
    return d[2].sum(), np.sqrt((d[3]**2).sum())

def collect(n_target):
    """delta -> {channel: (val, err)} for one filling."""
    out = {}
    for f in os.listdir(BASE):
        m = PAT.match(f)
        if not m:
            continue
        if int(m.group(1)) != L or abs(float(m.group(3)) - U) > 1e-6:
            continue
        n, tA, tt = float(m.group(2)), float(m.group(4)), float(m.group(5))
        if abs(n - n_target) > 0.02 or abs(tt - TT) > 0.02:
            continue
        kdir = os.path.join(BASE, f, 'dir-kVals')
        d = -tA
        out[d] = {}
        for fname, *_ in CHANNELS:
            fp = os.path.join(kdir, f'Vertex_{fname}.dat')
            if os.path.isfile(fp):
                out[d][fname] = vsum(fp)
    return out

def crossing(ds, y1, y2):
    """delta where y2 (dxy) overtakes y1 (dx2-y2), linear interp; None if none."""
    diff = np.array(y2) - np.array(y1)
    for i in range(len(diff) - 1):
        if diff[i] <= 0 <= diff[i+1] and diff[i+1] != diff[i]:
            return ds[i] - diff[i] * (ds[i+1] - ds[i]) / (diff[i+1] - diff[i])
    return None

data622 = collect(0.622)
data582 = collect(0.582)

fig, (axa, axb) = plt.subplots(1, 2, figsize=(14, 6), facecolor='white')

# ---- panel (a): all channels at n=0.622 ----
ds = sorted(data622)
for fname, label, mk, col in CHANNELS:
    xs = [d for d in ds if fname in data622[d]]
    ys = [data622[d][fname][0] for d in xs]
    es = [data622[d][fname][1] for d in xs]
    axa.errorbar(xs, ys, yerr=es, marker=mk, ms=10, lw=1.8, capsize=4,
                 color=col, markerfacecolor='w', markeredgecolor=col,
                 markeredgewidth=1.6, label=label, zorder=5)
axa.axhline(0, color='gray', ls='--', lw=1.0, alpha=0.6)
axa.set_xlabel(r'Anisotropy $\delta$', fontsize=20)
axa.set_ylabel(r'$N_{\zeta}^{\mathrm{Vertex}}$', fontsize=22)
axa.set_title(r'(a) all channels,  $n=0.622$', fontsize=17)
axa.tick_params(axis='both', labelsize=14)
axa.legend(fontsize=14, ncol=2, frameon=False, loc='upper left')

# ---- panel (b): d-wave crossover, both fillings ----
def dwave_curves(dat):
    ds = sorted(dat)
    x  = [d for d in ds if 'dwave' in dat[d] and 'dd12wave' in dat[d]]
    y1 = [dat[d]['dwave'][0]    for d in x]   # dx2-y2
    e1 = [dat[d]['dwave'][1]    for d in x]
    y2 = [dat[d]['dd12wave'][0] for d in x]   # dxy
    e2 = [dat[d]['dd12wave'][1] for d in x]
    return x, y1, e1, y2, e2

x6, y16, e16, y26, e26 = dwave_curves(data622)
x5, y15, e15, y25, e25 = dwave_curves(data582)

axb.errorbar(x6, y16, yerr=e16, marker='o', ms=10, lw=2.2, capsize=4, color='#1f77b4',
             markerfacecolor='#1f77b4', label=r'$d_{x^2-y^2}$, $n=0.622$', zorder=6)
axb.errorbar(x6, y26, yerr=e26, marker='^', ms=10, lw=2.2, capsize=4, color='#2ca02c',
             markerfacecolor='#2ca02c', label=r'$d_{xy}$, $n=0.622$', zorder=6)
axb.errorbar(x5, y15, yerr=e15, marker='o', ms=9, lw=1.8, ls='--', capsize=3, color='#1f77b4',
             markerfacecolor='w', alpha=0.6, label=r'$d_{x^2-y^2}$, $n=0.582$', zorder=5)
axb.errorbar(x5, y25, yerr=e25, marker='^', ms=9, lw=1.8, ls='--', capsize=3, color='#2ca02c',
             markerfacecolor='w', alpha=0.6, label=r'$d_{xy}$, $n=0.582$', zorder=5)

c6 = crossing(x6, y16, y26)
c5 = crossing(x5, y15, y25)
cs = [c for c in (c5, c6) if c is not None]
if cs:
    axb.axvspan(min(cs), max(cs), color='gold', alpha=0.18, zorder=0)
    axb.text(np.mean(cs), axb.get_ylim()[1]*0.05, 'crossover',
             ha='center', va='bottom', fontsize=13, color='#8a6d00')
    for c in cs:
        axb.axvline(c, color='gray', ls=':', lw=1.2, alpha=0.7, zorder=1)
print("crossover delta: n=0.622 ->", c6, " n=0.582 ->", c5)

axb.axhline(0, color='gray', ls='--', lw=1.0, alpha=0.6)
axb.set_xlabel(r'Anisotropy $\delta$', fontsize=20)
axb.set_ylabel(r'$N_{\zeta}^{\mathrm{Vertex}}$', fontsize=22)
axb.set_title(r'(b) $d$-wave crossover', fontsize=17)
axb.tick_params(axis='both', labelsize=14)
axb.legend(fontsize=13, frameon=False, loc='upper left')

plt.tight_layout()
plt.savefig("Fig6_pairing_crossover.png", dpi=600, bbox_inches='tight', facecolor='white')
print("saved Fig6_pairing_crossover.png")
