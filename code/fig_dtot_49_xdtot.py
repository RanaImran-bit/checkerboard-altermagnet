"""Delta_tot on the x-axis, with the third parameter carried by colour.

  (a) x = Delta_tot, y = U/t,   colour = anisotropy delta
  (b) x = Delta_tot, y = delta, colour = U/t

WHY THESE ARE POINTS AND NOT A FILLED MAP. Putting the measured quantity on x
makes the control parameters the dependent ones, and Delta_tot is not monotonic
in either. At fixed U it rises with delta, peaks on the ridge, then falls, so a
single Delta_tot value corresponds to TWO different delta (checked: every U from
2 to 5 has at least one turning point, and delta = 0.1, 0.2, 0.3 turn in U as
well). A filled colormap over this plane would have to choose one branch and
interpolate colour through regions where the function is double valued, i.e.
invent data. The 49 measured points are drawn exactly where they were measured.

Colour comes from a continuous colormap read off a continuous colorbar, so the
family trend appears as a gradient rather than as discrete legend swatches.

Delta_tot = sum_k |n_up(k) - n_dn(k)| / L^2, PRL Eq. (2). U = 0 is excluded: it
is identically zero at every delta and would pile 7 points on the origin.
"""
import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt

mpl.rcParams.update({
    "font.family": "serif", "font.serif": ["DejaVu Serif"],
    "mathtext.fontset": "dejavuserif", "axes.linewidth": 1.1,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
    "axes.labelsize": 16, "xtick.labelsize": 13, "ytick.labelsize": 13,
    "figure.dpi": 130})

D = "/Users/liujiaxin/Desktop/checkerboard-altermagnet"
d = pd.read_csv(f"{D}/data/fortran_L12_all49.csv")
d = d[np.isclose(d.n, 1.0) & (d.U >= 2)]

fig, ax = plt.subplots(1, 2, figsize=(13.6, 5.0))

s0 = ax[0].scatter(d.dtot_N, d.U, c=d.delta, cmap="plasma", s=330,
                   edgecolor="0.2", linewidth=0.7, zorder=3)
ax[0].set_xlabel(r"$\Delta_{\mathrm{tot}}$"); ax[0].set_ylabel(r"$U/t$")
ax[0].set_title("(a)", loc="left", fontsize=15)
cb0 = fig.colorbar(s0, ax=ax[0], pad=0.02)
cb0.set_label(r"anisotropy  $\delta$", fontsize=15)

s1 = ax[1].scatter(d.dtot_N, d.delta, c=d.U, cmap="viridis", s=330,
                   edgecolor="0.2", linewidth=0.7, zorder=3)
ax[1].set_xlabel(r"$\Delta_{\mathrm{tot}}$")
ax[1].set_ylabel(r"anisotropy  $\delta$")
ax[1].set_title("(b)", loc="left", fontsize=15)
cb1 = fig.colorbar(s1, ax=ax[1], pad=0.02)
cb1.set_label(r"$U/t$", fontsize=15)

for a in ax:
    a.grid(alpha=0.18, lw=0.7)
    a.set_axisbelow(True)

fig.tight_layout()
fig.savefig(f"{D}/figs/fig_dtot_49_xdtot.png", dpi=600, bbox_inches="tight",
            facecolor="white")

print(f"{len(d)} measured points (U >= 2)")
print(f"Delta_tot range {d.dtot_N.min():.4f} to {d.dtot_N.max():.4f}")
