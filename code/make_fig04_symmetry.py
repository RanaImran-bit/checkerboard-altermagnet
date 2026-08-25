"""Figure 4: symmetry classification of the momentum-resolved spin polarization.

Two dimensionless ratios, built only from the symmetry of Dn(k) and therefore
insensitive to how the weight is spread among higher harmonics of the same
symmetry (a projection onto a single lattice harmonic is not):

    A_odd = sum_k |Dn(k) + Dn(R90 k)| / sum_k |Dn(k)|    0 => C4-ODD  => d-wave
    M_odd = sum_k |Dn(k) + Dn(M k)|   / sum_k |Dn(k)|    0 => d_x2-y2, 2 => d_xy

Panel (a) places every interacting parameter set at L = 8, 10, 12 in the
(M_odd, A_odd) plane, where the two d-wave channels sit at opposite corners.
Panels (b)-(d) show why: a representative Dn(k), its C4 sum, which cancels, and
its mirror sum, which does not.

Partner momenta are found by explicit lookup on folded coordinates. np.rot90 is
NOT the C4 map on an FFT-indexed mesh, and returns 0.42-0.63 for pure d-wave
harmonics that must give exactly 0.

  python make_fig04_symmetry.py
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

plt.rcParams.update({
    "font.family": "serif", "font.serif": ["DejaVu Serif"],
    "font.size": 8, "axes.labelsize": 9, "legend.fontsize": 7,
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
    "axes.linewidth": 0.7, "xtick.major.width": 0.7, "ytick.major.width": 0.7,
    "figure.dpi": 600, "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
})
CL = {8: "#0072B2", 10: "#D55E00", 12: "#009E73"}     # validated CVD-safe
MK = {8: "o", 10: "s", 12: "^"}
HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "..", "data")

cls = pd.read_csv(os.path.join(D, "classify_3L.csv"))
cls = cls[~cls.A_odd.isna()]

def snap(a):
    a = np.asarray(a, float)
    a = np.where(np.abs(np.abs(a) - np.pi) < 1e-4, -np.pi, a)
    return (a + np.pi) % (2 * np.pi) - np.pi

def grid_and_partners(arr):
    """Return Dn on a regular grid plus its C4 and mirror partner grids."""
    kx, ky, v = snap(arr[:, 0]), snap(arr[:, 1]), arr[:, 2].astype(float)
    seen = {}
    for a, b, val in zip(kx, ky, v):
        seen[(round(a, 6), round(b, 6))] = val
    ks = sorted({k[0] for k in seen})
    L = len(ks); idx = {k: i for i, k in enumerate(ks)}
    G = np.full((L, L), np.nan)
    for (a, b), val in seen.items():
        G[idx[a], idx[b]] = val
    look = {k: i for i, k in enumerate(seen)}
    def partner(op):
        P = np.full((L, L), np.nan)
        for (a, b), val in seen.items():
            pa, pb = op(a, b)
            key = (round(float(snap(np.array([pa]))[0]), 6),
                   round(float(snap(np.array([pb]))[0]), 6))
            if key in seen:
                P[idx[a], idx[b]] = seen[key]
        return P
    return ks, G, partner(lambda a, b: (-b, a)), partner(lambda a, b: (b, a))

maps = np.load(os.path.join(D, "fortran_L12_maps.npz"), allow_pickle=True)
meta = maps["meta"]
pick = int(np.argmin(np.abs(meta[:, 2] - 4.0) + np.abs(meta[:, 3] - 0.4) + 100*(meta[:, 2] == 0)))
ks, G, G4, GM = grid_and_partners(np.asarray(maps["maps"][pick], float))
U_s, d_s = meta[pick, 2], meta[pick, 3]

fig = plt.figure(figsize=(6.75, 2.55))
gs = GridSpec(1, 4, width_ratios=[1.5, 1, 1, 1], wspace=0.42,
              top=0.80, bottom=0.20, figure=fig)

ax = fig.add_subplot(gs[0, 0])
for L in sorted(cls.L.unique()):
    s = cls[cls.L == L]
    ax.plot(s.M_odd, s.A_odd, MK[int(L)], color=CL[int(L)], mfc="none", mew=0.9,
            ms=3.6, ls="none", label=rf"$L={int(L)}$")
ax.plot(2, 0, "*", color="0.25", ms=11, zorder=5)
ax.annotate(r"$d_{xy}$", xy=(2, 0), xytext=(2.10, 0.028), fontsize=8.5,
            color="0.2", ha="left")
ax.plot(0, 0, "*", color="0.6", ms=9, zorder=5)
ax.annotate(r"$d_{x^2-y^2}$", xy=(0, 0), xytext=(0.02, 0.075), fontsize=8.5,
            color="0.45", ha="left")
ax.set_xlabel(r"$M_{\mathrm{odd}}$"); ax.set_ylabel(r"$A_{\mathrm{odd}}$")
ax.set_xlim(-0.15, 2.48); ax.set_ylim(-0.035, 0.50)
ax.tick_params(direction="in", top=True, right=True)
ax.legend(frameon=False, loc="upper left", ncol=1, handlelength=1.0,
          columnspacing=0.8, borderaxespad=0.3, labelspacing=0.15)
ax.set_title("(a)", loc="left", fontsize=9, pad=6)

ext = [-np.pi, np.pi, -np.pi, np.pi]
vmax = np.nanmax(np.abs(G))
for j, (M, tt) in enumerate([(G, r"$\Delta n(\mathbf{k})$"),
                             (G + G4, r"$\Delta n(\mathbf{k})+\Delta n(R_{90}\mathbf{k})$"),
                             (G + GM, r"$\Delta n(\mathbf{k})+\Delta n(M\mathbf{k})$")]):
    a2 = fig.add_subplot(gs[0, j + 1])
    im = a2.imshow(M.T, origin="lower", extent=ext, cmap="RdBu_r",
                   vmin=-vmax, vmax=vmax, interpolation="nearest")
    a2.set_title(f"({'bcd'[j]})", loc="left", fontsize=9, pad=12)
    a2.text(0.5, 1.045, tt, transform=a2.transAxes, ha="center", fontsize=6.8)
    a2.set_xticks([-np.pi, 0, np.pi]); a2.set_yticks([-np.pi, 0, np.pi])
    a2.set_xticklabels([r"$-\pi$", "0", r"$\pi$"])
    a2.set_yticklabels([r"$-\pi$", "0", r"$\pi$"] if j == 0 else [])
    a2.set_xlabel(r"$k_x$")
    if j == 0: a2.set_ylabel(r"$k_y$", labelpad=1)
    a2.tick_params(direction="in", top=True, right=True, colors="0.3")
    if j == 2:
        cb = fig.colorbar(im, ax=a2, fraction=0.046, pad=0.04)
        cb.ax.tick_params(labelsize=6.5, width=0.6)
        cb.outline.set_linewidth(0.6)

fig.text(0.985, 0.005, rf"$L=12$, $U={U_s:g}$, $\delta={d_s:g}$",
         ha="right", fontsize=6.8, color="0.35")
out = os.path.join(HERE, "..", "manuscript", "figures", "fig04_symmetry")
fig.savefig(out + ".pdf"); fig.savefig(out + ".png")
print(f"wrote {out}.pdf / .png   (panel b-d cell: U={U_s:g}, delta={d_s:g})")
print(f"A_odd {cls.A_odd.min():.4f}-{cls.A_odd.max():.4f}   "
      f"M_odd {cls.M_odd.min():.4f}-{cls.M_odd.max():.4f}   n={len(cls)}")
print(f"C4 sum / |Dn| = {np.nansum(np.abs(G+G4))/np.nansum(np.abs(G)):.4f}   "
      f"mirror sum / |Dn| = {np.nansum(np.abs(G+GM))/np.nansum(np.abs(G)):.4f}")
