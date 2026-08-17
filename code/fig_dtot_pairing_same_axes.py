"""Delta_tot and the pairing channels on the SAME axes, plus the phase diagram.

Both quantities come from the same Fortran run folders, so they can be plotted
in the same (U, delta) and (n, delta) planes with no cross-code caveat.

Top row is the (U, delta) plane at L=14 half filling -- the direct analogue of
PRL Fig. 1, whose panel (a) is Delta_tot and panel (b) is the d-wave vertex.
Bottom row is the (n, delta) plane at L=10, U=4.

The phase-diagram panels use DISCRETE colour: one flat colour per cell, no
blending. An earlier version blended the leading two channels, which produced
violet and teal patches that match no channel in the legend and read as extra
phases.
"""
import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch
from scipy.interpolate import griddata

mpl.rcParams.update({
    "font.family": "serif", "font.serif": ["DejaVu Serif"],
    "mathtext.fontset": "dejavuserif", "axes.linewidth": 1.2,
    "xtick.direction": "out", "ytick.direction": "out",
    "xtick.top": False, "ytick.right": False,
    "axes.labelsize": 15, "xtick.labelsize": 12, "ytick.labelsize": 12,
    "figure.dpi": 130})

D = "/Users/liujiaxin/Desktop/checkerboard-altermagnet/data"
dt = pd.read_csv(f"{D}/fortran_dnk.csv");  dt["dtot_N"] = dt.delta_tot / dt.L**2
pr = pd.read_csv(f"{D}/fortran_pairing_all.csv")
d  = dt.merge(pr, on=["L", "n", "U", "delta"])
K   = ["son", "sext", "d", "dxy"]
LAB = ["on-site $s$", "extended $s$", r"$d_{x^2-y^2}$", r"$d_{xy}$"]
COL = ["#8C8C8C", "#1DB954", "#0047FF", "#FF1500"]

ROWS = [(d[(d.L == 14) & np.isclose(d.n, 1.0)], "delta", "U",
         r"anisotropy $\delta$", r"$U/t$", r"$L=14$, half filling"),
        (d[(d.L == 10) & (d.U == 4.0)],        "n",     "delta",
         r"filling $n$", r"anisotropy $\delta$", r"$L=10$, $U=4$")]
FIELDS = [("dtot_N", r"$\Delta_{tot}$"), ("d", r"$d_{x^2-y^2}$ vertex"),
          ("dxy", r"$d_{xy}$ vertex")]

fig, ax = plt.subplots(2, 4, figsize=(21.5, 9.4))
for r, (T, xk, yk, xl, yl, ttl) in enumerate(ROWS):
    gx, gy = np.meshgrid(np.linspace(T[xk].min(), T[xk].max(), 260),
                         np.linspace(T[yk].min(), T[yk].max(), 260))
    for c, (fld, flab) in enumerate(FIELDS):
        a = ax[r, c]
        gi = griddata((T[xk].values, T[yk].values), T[fld].values, (gx, gy),
                      method="linear")
        im = a.contourf(gx, gy, gi, levels=100, cmap="jet", extend="both")
        cb = fig.colorbar(im, ax=a, pad=0.02)
        cb.set_label(flab, fontsize=13); cb.ax.tick_params(labelsize=10)
        a.plot(T[xk], T[yk], "o", ms=4, mfc="none", mec="k", mew=0.9)
        a.set_title(f"({chr(97+4*r+c)}) {flab}", fontsize=14)
        a.set_xlabel(xl); a.set_ylabel(yl)
    # phase diagram: DISCRETE, one colour per measured cell
    a = ax[r, 3]
    win = T[K].values.argmax(axis=1)
    a.scatter(T[xk], T[yk], c=[COL[i] for i in win], s=320, marker="s",
              edgecolors="k", linewidths=0.6)
    a.set_title(f"({chr(97+4*r+3)}) leading channel", fontsize=14)
    a.set_xlabel(xl); a.set_ylabel(yl)
    a.margins(0.09)
    print(f"{ttl}: " + ", ".join(
        f"{LAB[i].replace('$','')} {int((win==i).sum())}" for i in range(4)))
fig.legend(handles=[Patch(fc=c, ec="k", lw=0.6, label=l) for c, l in zip(COL, LAB)],
           fontsize=13, loc="lower center", ncol=4, bbox_to_anchor=(0.5, -0.035),
           frameon=False)
for r, (_, _, _, _, _, ttl) in enumerate(ROWS):
    ax[r, 0].text(-0.30, 0.5, ttl, transform=ax[r, 0].transAxes, fontsize=15,
                  rotation=90, va="center", ha="center")
fig.tight_layout()
fig.savefig(f"{D}/../figs/fig_dtot_pairing_same_axes.png", dpi=600,
            bbox_inches="tight", facecolor="white")
