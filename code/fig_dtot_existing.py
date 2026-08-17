"""Delta_tot from the runs that EXIST, with no interpolation.

The (U, delta) grid holds 19 of 35 cells, so a contour heatmap fills only the
convex hull and leaves a white diamond that reads as missing physics rather than
missing runs. These panels use every measured point directly and invent nothing:

  (a) Delta_tot against delta, one line per U
  (b) Delta_tot against U, one line per delta
  (c) the grid itself, as discrete cells -- present cells coloured, absent blank

Delta_tot = sum_k |n_up(k) - n_dn(k)| / L^2, the PRL definition.
"""
import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt

mpl.rcParams.update({
    "font.family": "serif", "font.serif": ["DejaVu Serif"],
    "mathtext.fontset": "dejavuserif", "axes.linewidth": 1.2,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True, "legend.frameon": False,
    "axes.labelsize": 16, "xtick.labelsize": 13, "ytick.labelsize": 13,
    "figure.dpi": 130})

D = "/Users/liujiaxin/Desktop/checkerboard-altermagnet/data"
d = pd.read_csv(f"{D}/fortran_dnk.csv")
d["dtot_N"] = d.delta_tot / d.L ** 2
h = d[(d.L == 14) & np.isclose(d.n, 1.0)].sort_values(["U", "delta"])

US = sorted(h.U.unique()); DS = sorted(h.delta.unique())
CU = plt.cm.viridis(np.linspace(0.05, 0.9, len(US)))
CD = plt.cm.plasma(np.linspace(0.05, 0.85, len(DS)))
MK = ["o", "s", "^", "D", "v", "P", "X"]

fig, ax = plt.subplots(1, 1, figsize=(7.6, 5.6))
ax = [ax]

# (a) against delta
for i, U in enumerate(US):
    s = h[h.U == U]
    ax[0].plot(s.delta, s.dtot_N, marker=MK[i % len(MK)], ms=7, lw=1.7,
               color=CU[i], label=f"$U={U:g}$")
ax[0].set_xlabel(r"anisotropy  $\delta$"); ax[0].set_ylabel(r"$\Delta_{tot}$")
ax[0].legend(fontsize=10.5, ncol=2, loc="upper left")

ax[0].set_title(r"$\Delta_{tot} = \sum_{\mathbf{k}}|n_\uparrow - n_\downarrow|/L^2$"
                "\n" r"$L=14$, half filling", fontsize=16)
fig.tight_layout()
fig.savefig(f"{D}/../figs/fig_dtot_vs_delta.png", dpi=600, bbox_inches="tight",
            facecolor="white")

for U in US:
    t = h[h.U == U].sort_values("delta")
    print(f"  U={U:g}: " + "  ".join(f"d={r.delta:g}:{r.dtot_N:.4f}" for _, r in t.iterrows()))
