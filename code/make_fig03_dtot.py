"""Figure 3: emergence of the momentum-space spin polarization.

(a) Delta_tot/N over the full 7 x 7 (U, delta) grid at L = 12. The U = 0 row is
    identically zero: with no auxiliary fields the Green function is the free one
    and n_up(k) = n_dn(k) exactly, so every feature above it is interaction
    generated.
(b) Cuts at fixed U, showing how the polarization varies with anisotropy. This is
    the complement of Fig. 5, which cuts the same surface at fixed delta.

Sequential single-hue map for a positive-definite magnitude (never a rainbow;
never a diverging map for data that does not change sign).

  python make_fig03_dtot.py
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.family": "serif", "font.serif": ["DejaVu Serif"],
    "font.size": 8, "axes.labelsize": 9, "legend.fontsize": 7,
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
    "axes.linewidth": 0.7, "xtick.major.width": 0.7, "ytick.major.width": 0.7,
    "lines.linewidth": 1.3, "lines.markersize": 4.2,
    "figure.dpi": 600, "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
})
HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "..", "data")
d = pd.read_csv(os.path.join(D, "fortran_L12_all49.csv"))

Us = sorted(d.U.unique()); Ds = sorted(d.delta.unique())
M = np.full((len(Ds), len(Us)), np.nan)
for i, dl in enumerate(Ds):
    for j, U in enumerate(Us):
        s = d[(d.delta == dl) & (d.U == U)]["dtot_N"]
        if len(s): M[i, j] = s.mean()

fig, (a1, a2) = plt.subplots(1, 2, figsize=(6.75, 2.6),
                             gridspec_kw={"width_ratios": [1.15, 1], "wspace": 0.40})

im = a1.imshow(M, origin="lower", aspect="auto", cmap="Blues",
               interpolation="nearest", vmin=0.0, vmax=np.nanmax(M))
a1.set_xticks(range(len(Us))); a1.set_xticklabels([f"{u:g}" for u in Us])
a1.set_yticks(range(len(Ds))); a1.set_yticklabels([f"{x:g}" for x in Ds])
a1.set_xlabel(r"$U/|t_0|$"); a1.set_ylabel(r"$\delta$")
a1.set_title("(a)", loc="left", fontsize=9)
a1.tick_params(direction="out", length=2.5)
# mark the exactly-zero noninteracting column
for i in range(len(Ds)):
    a1.text(0, i, "0", ha="center", va="center", fontsize=6.2, color="0.45")
cb = fig.colorbar(im, ax=a1, fraction=0.046, pad=0.03)
cb.set_label(r"$\Delta_{\mathrm{tot}}/N$", fontsize=8)
cb.ax.tick_params(labelsize=6.8, width=0.6); cb.outline.set_linewidth(0.6)

CU = {2.0: "#0072B2", 3.5: "#D55E00", 4.5: "#009E73", 5.0: "#E69F00"}
MKU = {2.0: "o", 3.5: "s", 4.5: "^", 5.0: "D"}
for U in [2.0, 3.5, 4.5, 5.0]:
    s = d[(d.U == U)].sort_values("delta")
    a2.plot(s.delta, s.dtot_N, MKU[U] + "-", color=CU[U], mfc="white", mew=1.1,
            label=rf"$U={U:g}$")
s0 = d[d.U == 0].sort_values("delta")
a2.plot(s0.delta, s0.dtot_N, "x--", color="0.45", ms=4, mew=1.0, lw=0.9,
        label=r"$U=0$")
a2.set_xlabel(r"$\delta$"); a2.set_ylabel(r"$\Delta_{\mathrm{tot}}/N$", labelpad=2)
a2.set_title("(b)", loc="left", fontsize=9)
a2.set_xlim(0.03, 0.77); a2.set_ylim(-0.004, 0.105)
a2.tick_params(direction="in", top=True, right=True)
a2.grid(color="0.92", lw=0.5); a2.set_axisbelow(True)
a2.legend(frameon=False, ncol=2, handlelength=1.5, labelspacing=0.2,
          columnspacing=1.0, borderaxespad=0.3)

out = os.path.join(HERE, "..", "manuscript", "figures", "fig03_dtot_emergence")
fig.savefig(out + ".pdf"); fig.savefig(out + ".png")
print(f"wrote {out}.pdf / .png")
print(f"U=0 row: max |Delta_tot/N| = {np.nanmax(np.abs(M[:, 0])):.3e}  (exactly zero)")
print(f"interacting range: {np.nanmin(M[:,1:]):.4f} - {np.nanmax(M[:,1:]):.4f}")
