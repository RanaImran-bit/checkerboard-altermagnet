"""Pairing against the AFM structure factor.

TERMINOLOGY, corrected. S_AFM(q) is a TWO-BODY CORRELATION FUNCTION,
sum_ij exp(iq.(r_i - r_j)) <S^z_i S^z_j> / N. It is not a polarisation.

A polarisation would be the one-body <S^z_i>, and that is identically zero here:
every driver runs CPMC(L, L, nup, nup, ...), i.e. equal up and down electron
numbers, so the calculation sits in the S_z = 0 sector by construction.

sqrt(S/N) equals a sublattice magnetisation ONLY under long-range order. The
finite-size scaling fails that test: S(q*) grows by a factor 1.0 from L = 8 to
L = 12 where order would need 2.25, and m falls as 1/L. So the old "m" and
"Delta m" labels claimed something the data does not support, and this file
plots and names the correlation itself.

S is evaluated at the TRUE PEAK of S(q), not at a fixed (pi,pi): the ordering
vector moves with delta and filling.
"""
import numpy as np, pandas as pd, matplotlib.pyplot as plt

D = "/Users/liujiaxin/Desktop/checkerboard-altermagnet/data"
CH = [("chi_son", "on-site $s$"), ("chi_sext", "extended $s$"),
      ("chi_d", r"$d_{x^2-y^2}$"), ("chi_dxy", r"$d_{xy}$")]
MK = {2.: "o", 4.: "s", 6.: "^", 8.: "D"}

pk = pd.read_csv(f"{D}/peak_master.csv"); pk = pk[pk.L == 12]
p = pd.read_csv(f"{D}/pairing_master.csv"); p = p[p.L == 12]

gm = pk.groupby(["U", "n", "delta"]).S_peak.mean().reset_index()
u0 = gm[gm.U == 0][["n", "delta", "S_peak"]].rename(columns={"S_peak": "S0"})
gm = gm.merge(u0, on=["n", "delta"])
gm["dS"] = gm.S_peak - gm.S0          # interaction-induced part of the correlation

gp = p.groupby(["U", "n", "delta"])[[c for c, _ in CH]].mean().reset_index()
d = gm.merge(gp, on=["U", "n", "delta"])
d = d[d.U > 0]
NS = [x for x in sorted(d.n.unique()) if x > 0.6]


def make(xcol, xlab, fname):
    fig, ax = plt.subplots(len(NS), 4, figsize=(19.5, 4.15 * len(NS)),
                           facecolor="white", squeeze=False)
    sm = plt.cm.ScalarMappable(cmap="viridis", norm=plt.Normalize(0, 0.7))
    for i, nv in enumerate(NS):
        h = d[np.isclose(d.n, nv)]
        for j, (c, lab) in enumerate(CH):
            a = ax[i, j]
            for U, g in h.groupby("U"):
                g = g.sort_values("delta")
                a.plot(g[xcol], g[c] / 144, "-", color="0.75", lw=1.2, zorder=1)
                a.scatter(g[xcol], g[c] / 144, c=g.delta, cmap="viridis",
                          vmin=0, vmax=0.7, s=95, marker=MK[U], edgecolors="k",
                          linewidths=0.7, zorder=3, label=f"$U={U:g}$")
            a.axhline(0, color="k", lw=1.1, ls="--")
            a.tick_params(labelsize=11)
            if i == 0: a.set_title(lab, fontsize=17, pad=20)
            if i == len(NS) - 1: a.set_xlabel(xlab, fontsize=14)
            if j == 0: a.set_ylabel(rf"$n={nv:.3f}$" "\n" r"$\chi_\alpha/N$", fontsize=13)
            rr = [np.corrcoef(g[xcol], g[c])[0, 1] for _, g in h.groupby("U")]
            rv = float(np.mean(rr))
            a.text(1.0, 1.012, rf"$\bar r$ = {rv:+.2f}", transform=a.transAxes,
                   fontsize=12.5, fontweight="bold", ha="right", va="bottom",
                   color="darkred" if rv > 0 else "navy")
            hh, ll = a.get_legend_handles_labels()
            a.legend(hh[:4], ll[:4], fontsize=9, loc="best", framealpha=0.85)
    cb = fig.colorbar(sm, ax=ax, fraction=0.012, pad=0.010)
    cb.set_label(r"$\delta$", fontsize=16); cb.ax.tick_params(labelsize=12)
    fig.savefig(fname, dpi=600, facecolor="white", bbox_inches="tight")
    plt.close(fig)


make("S_peak", r"AFM structure factor  $S_{\rm AFM}(\mathbf{q}^*)$",
     f"{D}/../figs/fig_safm_allfillings.png")
make("dS", r"$\Delta S_{\rm AFM} = S_{\rm AFM}(U) - S_{\rm AFM}(U{=}0)$",
     f"{D}/../figs/fig_dsafm_allfillings.png")

print("within-U correlation with the AFM structure factor, averaged over U")
print(f"{'n':>7} {'x-axis':>10} | " + " ".join(f"{l.replace('$','')[:9]:>10}" for _, l in CH))
for nv in NS:
    for tag, xc in [("S_AFM", "S_peak"), ("dS_AFM", "dS")]:
        h = d[np.isclose(d.n, nv)]
        row = [np.mean([np.corrcoef(g[xc], g[c])[0, 1] for _, g in h.groupby("U")])
               for c, _ in CH]
        print(f"{nv:7.3f} {tag:>10} | " + " ".join(f"{x:10.2f}" for x in row))
    print()
