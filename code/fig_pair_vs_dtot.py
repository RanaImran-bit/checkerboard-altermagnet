"""All four pairing channels against the altermagnetic splitting, half filling.

Built on the complete L=12 Fortran grid: 7 U x 7 delta = 49 cells, n = 1.

Structure follows the reference figure, which plots pairing against a MAGNETIC
quantity on x. Here that quantity is Delta_tot rather than m. The literal
reading -- Delta_tot on x with U (or delta) on y -- cannot show the channels at
all: Delta_tot and U are properties of the cell, not of the channel, so all four
channels would sit on one identical point. So the channel value takes the y
axis, and the swept parameter is carried by colour:

  (a) colour = U/t      -- what happens as U changes
  (b) colour = delta    -- what happens as delta changes

Both panels show the same 42 cells and the same four channels. Channel is set by
marker shape, exactly as the reference sets U by marker shape.

Pairing values are the Fortran vertex contributions at k = (0,0), full - L*Unpair.
U = 0 is excluded: every channel is identically zero there and Delta_tot is too,
so it contributes a pile of points at the origin.
"""
import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.lines import Line2D

mpl.rcParams.update({
    "font.family": "serif", "font.serif": ["DejaVu Serif"],
    "mathtext.fontset": "dejavuserif", "axes.linewidth": 1.1,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
    "axes.labelsize": 16, "xtick.labelsize": 12, "ytick.labelsize": 12,
    "figure.dpi": 130})

CH = [("son", "on-site $s$", "v"), ("sext", "extended $s$", "o"),
      ("d", r"$d_{x^2-y^2}$", "s"), ("dxy", r"$d_{xy}$", "D")]

D = "/Users/liujiaxin/Desktop/checkerboard-altermagnet"
d = pd.read_csv(f"{D}/data/fortran_L12_all49.csv")
d = d[np.isclose(d.n, 1.0) & (d.U > 0)].sort_values(["U", "delta"])

fig, ax = plt.subplots(1, 2, figsize=(15.0, 5.6))
PANEL = [("U", r"$U/t$", "viridis"), ("delta", r"$\delta$", "plasma")]

for k, (ccol, clab, cmap) in enumerate(PANEL):
    a = ax[k]
    norm = plt.Normalize(d[ccol].min(), d[ccol].max())
    for c, lab, mk in CH:
        a.scatter(d.dtot_N, d[c], c=d[ccol], cmap=cmap, norm=norm, s=88,
                  marker=mk, edgecolors="k", linewidths=0.6, zorder=3)
    a.axhline(0, color="k", lw=1.1, ls="--")
    a.set_xlabel(r"$\Delta_{\mathrm{tot}}$")
    a.set_ylabel(r"pairing vertex at $\mathbf{k}=0$")
    a.set_title(f"({'ab'[k]})  colour = {clab}", loc="left", fontsize=15)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm); sm.set_array([])
    cb = fig.colorbar(sm, ax=a, pad=0.02); cb.set_label(clab, fontsize=15)
    a.legend(handles=[Line2D([], [], marker=mk, ls="none", mfc="0.85", mec="k",
                             ms=8, label=lab) for _, lab, mk in CH],
             fontsize=10, loc="best", framealpha=0.9)

fig.tight_layout()
fig.savefig(f"{D}/figs/fig_pair_vs_dtot.png", dpi=600, bbox_inches="tight",
            facecolor="white")

print(f"{len(d)} cells, half filling, U > 0")
print("\ncorrelation of each channel with Delta_tot (within fixed U, averaged):")
for c, lab, _ in CH:
    r = np.mean([np.corrcoef(g.dtot_N, g[c])[0, 1] for _, g in d.groupby("U")])
    print(f"   {lab:16s} r = {r:+.2f}")
print("\nsame, within fixed delta:")
for c, lab, _ in CH:
    r = np.mean([np.corrcoef(g.dtot_N, g[c])[0, 1] for _, g in d.groupby("delta")])
    print(f"   {lab:16s} r = {r:+.2f}")
