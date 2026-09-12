"""L=12 half filling, the complete 7 x 7 (U, delta) Fortran grid -- all 49 cells.

Every cell is measured, so these are true heatmaps: pcolormesh over the data
itself, no griddata, no convex hull, no white wedge where runs were missing.
That was the whole reason for running the grid rather than replotting the L=14
set, which held 19 of 35 cells and stopped at delta = 0.4.

The result the figure is built around: the altermagnetic splitting and the
pairing crossover move together. As U grows, the ridge of Delta_tot shifts to
larger delta, and so does the delta at which d_xy overtakes d_x2-y2. Panel (c)
puts the two boundaries on one axis.

Delta_tot = sum_k |n_up(k) - n_dn(k)| / L^2, PRL Eq. (2), absolute value inside
the sum so it stays finite while the net magnetisation vanishes. Pairing values
are the Fortran equal-time vertex at k = (0,0).

CAVEAT for the caption: the U = 0 row is identically zero, the correct
non-interacting control, but the onset in U at fixed delta tracks the UHF trial's
own Stoner threshold. Trends above onset are QMC; the onset position is inherited
from the trial.
"""
import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt

mpl.rcParams.update({
    "font.family": "serif", "font.serif": ["DejaVu Serif"],
    "mathtext.fontset": "dejavuserif", "axes.linewidth": 1.1,
    "xtick.direction": "out", "ytick.direction": "out",
    "legend.frameon": False, "axes.labelsize": 15,
    "xtick.labelsize": 12, "ytick.labelsize": 12, "figure.dpi": 130})

D = "/Users/liujiaxin/Desktop/checkerboard-altermagnet"
d = pd.read_csv(f"{D}/data/fortran_L12_half_49.csv")
d["dif"] = d.dxy - d.d
US, DS = np.sort(d.U.unique()), np.sort(d.delta.unique())
A = d.pivot(index="U", columns="delta", values="dtot_N").reindex(index=US, columns=DS).values
B = d.pivot(index="U", columns="delta", values="dif").reindex(index=US, columns=DS).values

fig, ax = plt.subplots(1, 3, figsize=(16.4, 4.7))

im = ax[0].pcolormesh(DS, US, A, cmap="jet", shading="nearest")
cb = fig.colorbar(im, ax=ax[0], pad=0.02)
cb.set_label(r"$\Delta_{\mathrm{tot}}$", fontsize=14)
ax[0].set_xlabel(r"anisotropy  $\delta$"); ax[0].set_ylabel(r"$U/t$")
ax[0].set_title(r"(a)  altermagnetic splitting", loc="left", fontsize=14)

v = np.nanmax(np.abs(B))
im2 = ax[1].pcolormesh(DS, US, B, cmap="coolwarm", vmin=-v, vmax=v, shading="nearest")
ax[1].contour(DS, US, B, levels=[0], colors="k", linewidths=2.2)
cb2 = fig.colorbar(im2, ax=ax[1], pad=0.02)
cb2.set_label(r"$d_{xy} - d_{x^2-y^2}$", fontsize=14)
ax[1].set_xlabel(r"anisotropy  $\delta$"); ax[1].set_ylabel(r"$U/t$")
ax[1].set_title(r"(b)  pairing channel,  red $= d_{xy}$", loc="left", fontsize=14)

# (c) the two boundaries together. The Delta_tot ridge is the delta of the
# maximum along each U row; the pairing boundary is the first sign change of
# d_xy - d_x2-y2, linearly interpolated between the two bracketing deltas.
ridge, cross, uu = [], [], []
for i, U in enumerate(US):
    if U == 0:
        continue
    uu.append(U)
    ridge.append(DS[np.nanargmax(A[i])])
    y, xc = B[i], None
    for j in range(len(y) - 1):
        if y[j] < 0 <= y[j + 1]:
            xc = DS[j] + (DS[j + 1] - DS[j]) * (-y[j]) / (y[j + 1] - y[j]); break
    cross.append(xc)
ax[2].plot(ridge, uu, "o-", ms=8, lw=2, color="tab:red",
           label=r"$\delta$ of max $\Delta_{\mathrm{tot}}$")
ax[2].plot(cross, uu, "s--", ms=8, lw=2, color="tab:blue",
           label=r"$d_{xy}$ overtakes $d_{x^2-y^2}$")
ax[2].set_xlabel(r"anisotropy  $\delta$"); ax[2].set_ylabel(r"$U/t$")
ax[2].set_xlim(0.05, 0.75); ax[2].legend(fontsize=11, loc="upper left")
ax[2].set_title("(c)  both boundaries move with $U$", loc="left", fontsize=14)

fig.tight_layout()
fig.savefig(f"{D}/figs/fig_L12_half49.png", dpi=600, bbox_inches="tight",
            facecolor="white")

print(f"{len(d)} cells, U={list(US)}, delta={list(DS)}")
print("\n  U   ridge(Delta_tot)   crossover(pairing)")
for U, r, c in zip(uu, ridge, cross):
    print(f"  {U:<4g} {r:^16g} {c if c is None else round(c,3):^18}")
