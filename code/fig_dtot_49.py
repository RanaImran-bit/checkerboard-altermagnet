"""Delta_tot across the complete L=12 half-filling grid, 7 U x 7 delta = 49 cells.

Delta_tot is on the x-axis in both panels, since it is the measured quantity and
U and delta are the control parameters:

  (a) delta against Delta_tot, one curve per U      -- what changes as U changes
  (b) U     against Delta_tot, one curve per delta  -- what changes as delta changes

Curve colour is drawn from a continuous colormap and read off a continuous
colorbar rather than a legend of discrete swatches, so the trend across the
family is visible as a gradient instead of seven unrelated boxes.

Delta_tot = sum_k |n_up(k) - n_dn(k)| / L^2, PRL Eq. (2). The absolute value is
inside the sum, so it stays finite even though the net magnetisation vanishes.
U = 0 is included as the control and sits exactly at zero for every delta.
"""
import numpy as np, pandas as pd, matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

mpl.rcParams.update({
    "font.family": "serif", "font.serif": ["DejaVu Serif"],
    "mathtext.fontset": "dejavuserif", "axes.linewidth": 1.1,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
    "axes.labelsize": 15, "xtick.labelsize": 12, "ytick.labelsize": 12,
    "figure.dpi": 130})

D = "/Users/liujiaxin/Desktop/checkerboard-altermagnet"
d = pd.read_csv(f"{D}/data/fortran_L12_all49.csv")
d = d[np.isclose(d.n, 1.0)]
US, DS = np.sort(d.U.unique()), np.sort(d.delta.unique())

fig, ax = plt.subplots(1, 2, figsize=(13.2, 5.0))

# ---- (a) one curve per U, coloured continuously by U ----
nu = Normalize(vmin=US.min(), vmax=US.max())
cmu = plt.cm.viridis
for U in US:
    s = d[d.U == U].sort_values("delta")
    ax[0].plot(s.dtot_N, s.delta, "-o", ms=5.5, lw=2.0, color=cmu(nu(U)),
               mec="0.25", mew=0.5)
ax[0].set_xlabel(r"$\Delta_{\mathrm{tot}}$")
ax[0].set_ylabel(r"anisotropy  $\delta$")
ax[0].set_title(r"(a)  curves of constant $U$", loc="left", fontsize=14)
cb0 = fig.colorbar(ScalarMappable(norm=nu, cmap=cmu), ax=ax[0], pad=0.02)
cb0.set_label(r"$U/t$", fontsize=14)

# ---- (b) one curve per delta, coloured continuously by delta ----
nd = Normalize(vmin=DS.min(), vmax=DS.max())
cmd = plt.cm.plasma
for dl in DS:
    s = d[d.delta == dl].sort_values("U")
    ax[1].plot(s.dtot_N, s.U, "-o", ms=5.5, lw=2.0, color=cmd(nd(dl)),
               mec="0.25", mew=0.5)
ax[1].set_xlabel(r"$\Delta_{\mathrm{tot}}$")
ax[1].set_ylabel(r"$U/t$")
ax[1].set_title(r"(b)  curves of constant $\delta$", loc="left", fontsize=14)
cb1 = fig.colorbar(ScalarMappable(norm=nd, cmap=cmd), ax=ax[1], pad=0.02)
cb1.set_label(r"anisotropy  $\delta$", fontsize=14)

fig.tight_layout()
fig.savefig(f"{D}/figs/fig_dtot_49.png", dpi=600, bbox_inches="tight",
            facecolor="white")

print("Delta_tot, complete 7x7 grid (rows U, cols delta):")
print(d.pivot_table(index="U", columns="delta", values="dtot_N").round(4).to_string())
print("\ndelta at which Delta_tot peaks, per U:")
for U in US[1:]:
    s = d[d.U == U]
    print(f"   U={U:g}: delta* = {s.loc[s.dtot_N.idxmax(), 'delta']:.1f}"
          f"   (max {s.dtot_N.max():.4f})")
