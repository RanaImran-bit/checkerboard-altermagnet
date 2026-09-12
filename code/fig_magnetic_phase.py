"""Magnetic ordering-vector phase diagram in the (n, delta) plane.

Built to mirror Fig. 4 of PRB 113, 134443, using the SAME Q coordinate:
  Q = 1        Neel      (pi, pi)
  Q = 0        stripe    (pi, 0) or (0, pi)
  0 < Q < 1    spiral    (pi, q) or (q, pi)

with delta on the vertical axis in place of t'. The peak of S_zz(q) is located
per seed on the full q grid, folded to the irreducible range, and Q is the
smaller of the two folded components.

This exists because magnetic_master.csv evaluates S at a hard-coded (pi,pi).
That is the right point only near half filling for delta = 0.2 to 0.5; elsewhere
the ordering vector has moved and the stored moment undercounts by up to 43%.
"""
import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import gaussian_filter

mpl.rcParams.update({
    "font.family": "serif", "font.serif": ["DejaVu Serif"],
    "mathtext.fontset": "dejavuserif", "axes.linewidth": 1.2,
    "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True,
    "ytick.right": True, "xtick.major.size": 6, "ytick.major.size": 6,
    "axes.labelsize": 19, "xtick.labelsize": 15, "ytick.labelsize": 15,
    "figure.dpi": 130})

D = "/Users/liujiaxin/Desktop/checkerboard-altermagnet/data"
d = pd.read_csv(f"{D}/peak_master.csv"); d = d[d.L == 12]
fx = np.where(d.qx > 1, 2 - d.qx, d.qx)
fy = np.where(d.qy > 1, 2 - d.qy, d.qy)
d["Q"] = np.minimum(fx, fy)                       # their coordinate exactly
g = d.groupby(["U", "n", "delta"]).Q.mean().reset_index()

US = [2.0, 4.0, 6.0, 8.0]
NS = np.sort(g.n.unique()); DS = np.sort(g.delta.unique())
nf = np.linspace(NS.min(), NS.max(), 400); df = np.linspace(DS.min(), DS.max(), 400)
NG, DG = np.meshgrid(nf, df); pts = np.stack([DG.ravel(), NG.ravel()], -1)

fig, ax = plt.subplots(1, 4, figsize=(20, 5.2), sharey=True)
for k, U in enumerate(US):
    a = ax[k]
    M = g[np.isclose(g.U, U)].pivot(index="delta", columns="n", values="Q").values
    F = gaussian_filter(
        RegularGridInterpolator((DS, NS), M)(pts).reshape(DG.shape), 7.0)
    im = a.pcolormesh(nf, df, F, cmap="viridis", shading="gouraud",
                      vmin=0, vmax=1, rasterized=True)
    a.set_title(rf"$U = {U:g}$", fontsize=22)
    a.set_xlabel(r"Filling  $n$")
    a.set_xticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    a.set_xlim(NS.min(), NS.max()); a.set_ylim(DS.min(), DS.max())
ax[0].set_ylabel(r"Anisotropy  $\delta$")
ax[0].text(0.60, 0.30, "spiral\n" r"$(\pi,q)$", fontsize=15, color="white",
           ha="center", va="center")
ax[3].text(0.62, 0.62, "stripe\n" r"$(\pi,0)$", fontsize=15, color="white",
           ha="center", va="center")
cb = fig.colorbar(im, ax=ax, pad=0.012, fraction=0.016)
cb.set_label(r"peak coordinate  $Q/\pi$", fontsize=17)
cb.set_ticks([0, 0.25, 0.5, 0.75, 1.0])
cb.ax.set_yticklabels(["0\nstripe", "0.25", "0.5\nspiral", "0.75", "1\nNéel"],
                      fontsize=13)
fig.savefig(f"{D}/../figs/fig_magnetic_phase_allU.png", dpi=600,
            bbox_inches="tight", facecolor="white")

print("majority ordering type per (n, delta), L=12, U=4 -- fraction of 6 seeds")
sub = d[d.U == 4]
print(f"{'delta':>6} | " + " ".join(f"{x:>8.3f}" for x in NS))
for dl in sorted(DS, reverse=True):
    row = []
    for nv in NS:
        s = sub[(sub.delta == dl) & np.isclose(sub.n, nv)]
        vc = s.kind.value_counts()
        row.append(f"{vc.index[0][:6]}" if len(vc) else "-")
    print(f"{dl:6.1f} | " + " ".join(f"{x:>8}" for x in row))
print("\ncensus over the whole L=12 cube:")
print(d.kind.value_counts().to_string())
