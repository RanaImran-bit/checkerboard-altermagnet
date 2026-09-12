"""Psi_dxy against delta at every available filling, with the wrong-symmetry control.

Psi_dxy = (1/N) sum_q sin(qx) sin(qy) S_zz(q). In real space this is the
difference between the "/" and "\\" diagonal spin correlations, which is what the
anisotropy makes inequivalent, so symmetry forces it to zero at delta = 0.

Psi_dx2y2 = (1/N) sum_q [cos(qx) - cos(qy)] S_zz(q) is the CONTROL. It lives on
the axis bonds, and delta never distinguishes x from y, so it must stay zero at
every delta. A bug, a normalisation slip or U-dependent noise would appear in
BOTH harmonics; the control staying flat is what excludes those.
"""
import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt

mpl.rcParams.update({
    "font.family": "serif", "font.serif": ["DejaVu Serif"],
    "mathtext.fontset": "dejavuserif", "axes.linewidth": 1.1,
    "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True,
    "ytick.right": True, "legend.frameon": False,
    "axes.labelsize": 14, "xtick.labelsize": 11, "ytick.labelsize": 11,
    "figure.dpi": 130})

D = "/Users/liujiaxin/Desktop/checkerboard-altermagnet/data"
MK = {2.: "o", 4.: "s", 6.: "^", 8.: "D"}
FC = {2.: "k", 4.: "white", 6.: "k", 8.: "white"}
LS = {2.: "-", 4.: "--", 6.: "-.", 8.: ":"}
US = [2., 4., 6., 8.]
KW = dict(color="k", mec="k", ms=4.5, lw=1.2, capsize=2.2)

d = pd.read_csv(f"{D}/magnetic_master.csv"); d = d[(d.L == 12) & (d.U > 0)]
NS = sorted(d.n.unique())

fig, ax = plt.subplots(2, len(NS), figsize=(3.5 * len(NS), 7.6),
                       sharex=True, squeeze=False)
for j, nv in enumerate(NS):
    h = d[np.isclose(d.n, nv)]
    g = h.groupby(["U", "delta"])[["Psi_dxy", "Psi_dx2y2"]].agg(["mean", "sem"])
    for row, col in enumerate(["Psi_dxy", "Psi_dx2y2"]):
        a = ax[row, j]
        for U in US:
            s = g.loc[U]
            a.errorbar(s.index, np.abs(s[(col, "mean")]) * 1e3,
                       yerr=s[(col, "sem")] * 1e3, marker=MK[U], mfc=FC[U],
                       ls=LS[U], label=f"$U={U:g}$", **KW)
        a.axhline(0, color="k", lw=1.0, ls="--")
        if row == 1: a.set_xlabel(r"$\delta$")
        if row == 0: a.set_title(rf"$n = {nv:.3f}$", fontsize=14)
    ax[0, j].set_ylim(-0.6, 3.4)          # one scale for every panel, both rows
    ax[1, j].set_ylim(-0.6, 3.4)
ax[0, 0].set_ylabel(r"$|\Psi_{d_{xy}}|\times 10^{3}$" "\n" "(a) AM order parameter")
ax[1, 0].set_ylabel(r"$|\Psi_{d_{x^2-y^2}}|\times 10^{3}$" "\n"
                    "(b) control: wrong symmetry")
ax[0, 0].legend(fontsize=9, loc="upper left")
fig.tight_layout()
fig.savefig(f"{D}/../figs/fig_psi_am_fillings.png", dpi=600,
            bbox_inches="tight", facecolor="white")

print("|Psi_dxy| x10^3 at delta=0.7, and its significance, per filling")
print(f"{'n':>7}" + "".join(f"{'U='+str(int(u)):>16}" for u in US))
for nv in NS:
    row = f"{nv:7.3f}"
    for U in US:
        s = d[(d.U == U) & np.isclose(d.n, nv) & (d.delta == 0.7)].Psi_dxy
        m, e = abs(s.mean()) * 1e3, s.sem() * 1e3
        row += f"{m:9.2f} ({m/e if e > 0 else 0:4.1f}s)"
    print(row)
print("\nsame at delta=0 -- symmetry demands zero:")
print(f"{'n':>7}" + "".join(f"{'U='+str(int(u)):>16}" for u in US))
for nv in NS:
    row = f"{nv:7.3f}"
    for U in US:
        s = d[(d.U == U) & np.isclose(d.n, nv) & (d.delta == 0.0)].Psi_dxy
        m, e = abs(s.mean()) * 1e3, s.sem() * 1e3
        row += f"{m:9.2f} ({m/e if e > 0 else 0:4.1f}s)"
    print(row)
