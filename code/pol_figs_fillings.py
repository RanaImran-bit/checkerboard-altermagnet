"""Delta polarisation against pairing, at every filling.

Dm = m - m(U=0) with m = sqrt(S_zz(pi,pi)/N). The U=0 subtraction removes the
free-fermion Pauli background, which at L=12 is 60-90% of m and carries no
information about interaction-driven magnetism.
"""
import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt
mpl.rcParams.update({
    "font.family": "serif", "font.serif": ["DejaVu Serif"],
    "mathtext.fontset": "dejavuserif", "axes.linewidth": 1.0,
    "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True,
    "ytick.right": True, "legend.frameon": False, "figure.dpi": 130,
    "axes.labelsize": 12, "xtick.labelsize": 10, "ytick.labelsize": 10})

D = "/Users/liujiaxin/Desktop/checkerboard-altermagnet/data"
CH = [("chi_sext", r"extended $s$", "#1b7837", "s"),
      ("chi_d", r"$d_{x^2-y^2}$", "#2166ac", "o"),
      ("chi_dxy", r"$d_{xy}$", "#b2182b", "^")]
MK = {2.: "o", 4.: "s", 6.: "^", 8.: "D"}

m = pd.read_csv(f"{D}/magnetic_master.csv"); m = m[m.L == 12]
p = pd.read_csv(f"{D}/pairing_master.csv"); p = p[p.L == 12]
gm = m.groupby(["U", "n", "delta"]).m.mean().reset_index()
u0 = gm[gm.U == 0][["n", "delta", "m"]].rename(columns={"m": "m0"})
gm = gm.merge(u0, on=["n", "delta"]); gm["dm"] = gm.m - gm.m0
gp = p.groupby(["U", "n", "delta"])[[c for c, *_ in CH]].mean().reset_index()
d = gm.merge(gp, on=["U", "n", "delta"])
d = d[d.U > 0]
NS = sorted(d.n.unique())

fig, ax = plt.subplots(2, 3, figsize=(14.6, 8.2))
for j, nv in enumerate(NS):
    a = ax.flat[j]; g0 = d[d.n == nv]
    a.axhline(0, color="k", lw=1.0, ls="--")
    txt = []
    for c, lab, col, mk in CH:
        for U, g in g0.groupby("U"):
            g = g.sort_values("delta")
            a.plot(g.dm, g[c] / 144, "-", color=col, lw=1.0, alpha=0.45, zorder=1)
            a.scatter(g.dm, g[c] / 144, color=col, s=42, marker=MK[U],
                      edgecolors="k", linewidths=0.5, zorder=3,
                      label=lab if U == 2.0 else None)
        rr = [np.corrcoef(g.dm, g[c])[0, 1] for _, g in g0.groupby("U")]
        txt.append((lab, np.mean(rr), col))
    a.set_title(rf"$n = {nv:.3f}$", fontsize=13)
    a.xaxis.set_major_locator(mpl.ticker.MaxNLocator(5))
    a.ticklabel_format(axis="x", style="sci", scilimits=(-2, 3))
    lines = "\n".join(f"{lab}:  r = {rv:+.2f}" for lab, rv, _ in txt)
    a.text(0.035, 0.035, lines, transform=a.transAxes, fontsize=9.5, va="bottom",
           bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="0.7", alpha=0.92))
    if j >= 3: a.set_xlabel(r"delta polarisation  $\Delta m$")
    if j % 3 == 0: a.set_ylabel(r"pairing vertex  $\chi_\alpha/N$")
ax.flat[0].legend(fontsize=9.5, loc="upper right")
fig.suptitle(r"Pairing against the delta polarisation $\Delta m = m - m(U{=}0)$ "
             r"at every filling ($L=12$)" "\n"
             r"lines join constant $U$; $\bar r$ is the correlation averaged over "
             r"the four $U$ values", fontsize=14, y=1.02)
fig.tight_layout()
fig.savefig("fig_pol_fillings.png", dpi=300, bbox_inches="tight", facecolor="white")

print("within-U correlation with Dm, averaged over U=2,4,6,8")
print(f"{'n':>7}{'ext-s':>9}{'dx2-y2':>9}{'dxy':>9}")
for nv in NS:
    g0 = d[d.n == nv]; row = []
    for c, *_ in CH:
        row.append(np.mean([np.corrcoef(g.dm, g[c])[0, 1] for _, g in g0.groupby("U")]))
    print(f"{nv:7.3f}" + "".join(f"{x:9.2f}" for x in row))
