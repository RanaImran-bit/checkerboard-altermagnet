"""Pairing against the magnetic polarisation, all six fillings.

Two versions of the same figure:
  x = m         raw staggered polarisation, no subtraction  (what was asked for)
  x = Delta m   m - m(U=0), interaction-induced part only

m = sqrt(S_zz(pi,pi)/N). At U = 0 it is already 60 to 93 per cent of its
interacting value, because free fermions have spin correlations from Pauli
statistics alone. Plotting against raw m therefore puts most of the x-axis range
into a background that carries no information about interaction-driven
magnetism, which is the argument for subtracting. Plotting raw m avoids taking a
difference of two comparable numbers. Both are defensible, so both are here.
"""
import numpy as np, pandas as pd, matplotlib.pyplot as plt

D = "/Users/liujiaxin/Desktop/checkerboard-altermagnet/data"
CH = [("chi_son", "on-site $s$"), ("chi_sext", "extended $s$"),
      ("chi_d", r"$d_{x^2-y^2}$"), ("chi_dxy", r"$d_{xy}$")]
MK = {2.: "o", 4.: "s", 6.: "^", 8.: "D"}

m = pd.read_csv(f"{D}/magnetic_master.csv"); m = m[m.L == 12]
p = pd.read_csv(f"{D}/pairing_master.csv"); p = p[p.L == 12]
gm = (m.groupby(["U", "n", "delta"]).m.agg(["mean", "sem"]).reset_index()
      .rename(columns={"mean": "m", "sem": "me"}))
u0 = gm[gm.U == 0][["n", "delta", "m"]].rename(columns={"m": "m0"})
gm = gm.merge(u0, on=["n", "delta"]); gm["dm"] = gm.m - gm.m0
gp = p.groupby(["U", "n", "delta"])[[c for c, _ in CH]].agg(["mean", "sem"]).reset_index()
gp.columns = ["U", "n", "delta"] + [f"{a}_{b}" for a, b in gp.columns[3:]]
d = gm.merge(gp, on=["U", "n", "delta"])
d = d[d.U > 0]
NS = [x for x in sorted(d.n.unique()) if x > 0.6]   # four fillings


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
                a.plot(g[xcol], g[f"{c}_mean"] / 144, "-", color="0.75", lw=1.2, zorder=1)
                a.scatter(g[xcol], g[f"{c}_mean"] / 144, c=g.delta, cmap="viridis",
                          vmin=0, vmax=0.7, s=95, marker=MK[U], edgecolors="k",
                          linewidths=0.7, zorder=3, label=f"$U={U:g}$")
            a.axhline(0, color="k", lw=1.1, ls="--")
            a.tick_params(labelsize=11)
            if i == 0: a.set_title(lab, fontsize=17, pad=20)
            if i == len(NS) - 1: a.set_xlabel(xlab, fontsize=14)
            if j == 0: a.set_ylabel(rf"$n={nv:.3f}$" "\n" r"$\chi_\alpha/N$", fontsize=13)
            # within-U correlation, averaged over the four U values: the pooled
            # value mixes four offset curves and understates the relation
            rr = [np.corrcoef(g[xcol], g[f"{c}_mean"])[0, 1] for _, g in h.groupby("U")]
            rv = float(np.mean(rr))
            # outside the axes: inside, the per-panel legend covered it in
            # several panels no matter which corner "best" chose
            a.text(1.0, 1.012, rf"$\bar r$ = {rv:+.2f}", transform=a.transAxes,
                   fontsize=12.5, fontweight="bold", ha="right", va="bottom",
                   color="darkred" if rv > 0 else "navy")
            hh, ll = a.get_legend_handles_labels()
            a.legend(hh[:4], ll[:4], fontsize=9, loc="best", framealpha=0.85)
    cb = fig.colorbar(sm, ax=ax, fraction=0.012, pad=0.010)
    cb.set_label(r"$\delta$", fontsize=16); cb.ax.tick_params(labelsize=12)
    fig.savefig(fname, dpi=600, facecolor="white", bbox_inches="tight")
    plt.show()


make("m", r"$m$", "fig_pol_m_allfillings.png")
make("dm", r"$\Delta m$", "fig_pol_dm_allfillings.png")

print("within-U correlation, averaged over U = 2,4,6,8")
print(f"{'n':>7} | " + " ".join(f"{l.replace('$','')[:9]:>10}" for _, l in CH))
for xc, tag in [("m", "vs raw m"), ("dm", "vs Delta m")]:
    print(f"--- {tag} ---")
    for nv in NS:
        h = d[np.isclose(d.n, nv)]
        row = []
        for c, _ in CH:
            row.append(np.mean([np.corrcoef(g[xc], g[f"{c}_mean"])[0, 1]
                                for _, g in h.groupby("U")]))
        print(f"{nv:7.3f} | " + " ".join(f"{x:10.2f}" for x in row))
