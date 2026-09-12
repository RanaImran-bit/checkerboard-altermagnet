"""Figure 5: the delta = 0 null control.

Delta_tot sums absolute values, so statistical noise contributes with a definite
sign and the estimator has a positive floor that does not vanish when the physics
says it should. At delta = 0 both checkerboard diagonals are equal on every site,
the two sublattices become equivalent, and the model reduces to the plain t-t'
square lattice whose Neel sublattices are related by a TRANSLATION. There
n_up(k) = n_dn(k) identically, so Delta_tot must vanish by symmetry and whatever is
measured is the floor at production settings.

The U = 0 column cannot serve this purpose: there the Green function is built
without auxiliary fields, so its zero is arithmetic rather than statistical.

Log y-axis, because the point of the figure is that the separation is two orders
of magnitude, not a few percent.

  python make_fig05_null.py
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# REVTeX single column is 3.375 in; keep fonts at ~8 pt so they survive reduction.
plt.rcParams.update({
    "font.family": "serif", "font.serif": ["DejaVu Serif"],
    "font.size": 8, "axes.labelsize": 9, "axes.titlesize": 9,
    "legend.fontsize": 7, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "axes.linewidth": 0.7, "xtick.major.width": 0.7, "ytick.major.width": 0.7,
    "lines.linewidth": 1.3, "lines.markersize": 4.5,
    "figure.dpi": 600, "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
})
# Okabe-Ito subset, validated: worst CVD pair separation 11.4 (target >= 8)
C = {0.1: "#0072B2", 0.3: "#D55E00", 0.5: "#009E73", 0.7: "#E69F00"}
MK = {0.1: "o", 0.3: "s", 0.5: "^", 0.7: "D"}

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
sig = pd.read_csv(os.path.join(D, "fortran_L12_all49.csv"))
flo = pd.read_csv(os.path.join(D, "delta0_control_L12.csv"))

fmean, fstd = flo.dtot_N.mean(), flo.dtot_N.std(ddof=1)

fig, ax = plt.subplots(figsize=(3.375, 2.7))

# the symmetry-forbidden floor, drawn as a band so the eye reads it as a region
ax.axhspan(max(fmean - fstd, 1e-6), fmean + fstd, color="0.75", alpha=0.55, lw=0, zorder=1)
ax.axhline(fmean, color="0.35", lw=0.9, ls="--", zorder=2)
ax.plot(flo.U, flo.dtot_N, "x", color="0.35", ms=4, mew=1.0, zorder=3,
        label=r"$\delta=0$ (symmetry forbidden)")

for dl in [0.1, 0.3, 0.5, 0.7]:
    s = sig[(sig.delta == dl) & (sig.U > 0)].sort_values("U")
    ax.plot(s.U, s.dtot_N, MK[dl] + "-", color=C[dl], mfc="white", mew=1.1,
            zorder=4, label=rf"$\delta={dl}$")

ax.set_yscale("log")
ax.set_xlabel(r"$U/|t_0|$")
ax.set_ylabel(r"$\Delta_{\mathrm{tot}}/N$")
ax.set_xlim(1.6, 5.4)
ax.set_ylim(2e-4, 6e-1)
ax.tick_params(direction="in", top=True, right=True, which="both")
ax.grid(axis="y", color="0.9", lw=0.5, zorder=0)
ax.set_axisbelow(True)

ax.annotate("noise floor", xy=(5.25, fmean), xytext=(4.25, 1.35e-3),
            fontsize=6.5, color="0.3",
            arrowprops=dict(arrowstyle="-", color="0.5", lw=0.6))
ax.legend(frameon=False, loc="upper left", ncol=2, handlelength=1.5,
          borderaxespad=0.2, labelspacing=0.2, columnspacing=1.0)

out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..",
                   "manuscript", "figures", "fig05_null_control")
fig.savefig(out + ".pdf")
fig.savefig(out + ".png")
print(f"wrote {out}.pdf / .png")
print(f"floor = {fmean:.6f} +/- {fstd:.6f}")
for dl in [0.1, 0.3, 0.5, 0.7]:
    s = sig[(sig.delta == dl) & (sig.U > 0)]
    print(f"  delta={dl}: signal {s.dtot_N.min():.5f}-{s.dtot_N.max():.5f}"
          f"  = {s.dtot_N.min()/fmean:5.1f}-{s.dtot_N.max()/fmean:5.1f} x floor")
