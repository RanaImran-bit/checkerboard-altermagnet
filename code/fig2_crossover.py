"""FIGURE 2 -- the pairing channel crossover driven by the plaquette anisotropy.

The paper's headline. At delta = 0 the checkerboard reproduces the isotropic
square-lattice result, extended-s and dx2-y2 leading. Turning on delta drives
dx2-y2 down monotonically and dxy up, and the two cross.

Two boundaries are NOT the same and the figure keeps them apart:

  (b) delta_d : where dxy overtakes dx2-y2. This is the symmetry statement --
      the checkerboard replaces the square lattice's d-wave with the other d.
  (c) delta*  : where dxy overtakes EVERY other channel. This is larger,
      because extended-s is the overall leader at small delta and dxy has to
      pass it too. Quoting delta_d as though it were delta* would overstate
      the result, so both are drawn.

chi is the tau-integrated pair susceptibility. The equal-time Fortran vertex
data gives the same crossover independently, which is the cross-check reported
in the supplement -- it is not a second version of this figure.

Error bars are the standard error over 6 independent seeds.
"""
import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt

mpl.rcParams.update({
    "font.family": "serif", "font.serif": ["DejaVu Serif"],
    "mathtext.fontset": "dejavuserif", "axes.linewidth": 1.1,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True, "legend.frameon": False,
    "axes.labelsize": 14, "xtick.labelsize": 11.5, "ytick.labelsize": 11.5,
    "figure.dpi": 130})

D = "/Users/liujiaxin/Desktop/checkerboard-altermagnet"
p = pd.read_csv(f"{D}/data/pairing_master.csv")
h = p[(p.L == 12) & np.isclose(p.n, 1.0)]
CH = ["chi_son", "chi_sext", "chi_d", "chi_dxy"]

m = h.groupby(["U", "delta"])[CH].mean().reset_index()
e = h.groupby(["U", "delta"])[CH].sem().reset_index()

# on-site s and extended s in neutral tones so the two d channels carry the eye;
# blue for dx2-y2 and red for dxy throughout the manuscript
STY = {"chi_son":  ("0.60", "^", r"on-site $s$"),
       "chi_sext": ("0.15", "o", r"extended $s$"),
       "chi_d":    ("tab:blue", "s", r"$d_{x^2-y^2}$"),
       "chi_dxy":  ("tab:red", "D", r"$d_{xy}$")}


def cross(x, y):
    """First sign change of y, linearly interpolated. None if it never crosses."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    for i in range(len(y) - 1):
        if y[i] < 0 <= y[i + 1]:
            return x[i] + (x[i + 1] - x[i]) * (-y[i]) / (y[i + 1] - y[i])
    return None


fig, ax = plt.subplots(1, 3, figsize=(15.6, 4.6))

# ---------------- (a) the four channels at U = 4 ----------------
s, se = m[m.U == 4].sort_values("delta"), e[e.U == 4].sort_values("delta")
for c in CH:
    col, mk, lab = STY[c]
    ax[0].errorbar(s.delta, s[c], yerr=se[c], marker=mk, ms=5.5, lw=1.7,
                   color=col, capsize=2.5, elinewidth=0.9, label=lab)
ax[0].axhline(0, color="k", lw=0.7, ls=":")
dd = cross(s.delta.values, (s.chi_dxy - s.chi_d).values)
ds = cross(s.delta.values, (s.chi_dxy - s[["chi_son", "chi_sext", "chi_d"]].max(axis=1)).values)
for xv, lab, ls in ((dd, r"$\delta_d$", "--"), (ds, r"$\delta^*$", "-.")):
    if xv is not None:
        ax[0].axvline(xv, color="0.45", lw=1.0, ls=ls)
        ax[0].text(xv, ax[0].get_ylim()[1] * 0.97, lab, ha="center", va="top",
                   fontsize=12, color="0.3")
ax[0].set_xlabel(r"anisotropy  $\delta$")
ax[0].set_ylabel(r"$\chi_\alpha$")
ax[0].legend(fontsize=10.5, loc="lower right")
ax[0].set_title(r"(a)  $U=4$", loc="left", fontsize=14)

# ---------------- (b) dxy against dx2-y2, one line per U ----------------
UU = [2.0, 4.0, 6.0, 8.0]
CU = plt.cm.viridis(np.linspace(0.05, 0.82, len(UU)))
dc = {}
for i, U in enumerate(UU):
    s = m[m.U == U].sort_values("delta"); se = e[e.U == U].sort_values("delta")
    y = s.chi_dxy - s.chi_d
    ye = np.hypot(se.chi_dxy, se.chi_d)          # independent seeds, add in quadrature
    ax[1].errorbar(s.delta, y, yerr=ye, marker="o", ms=5, lw=1.7, color=CU[i],
                   capsize=2.5, elinewidth=0.9, label=f"$U={U:g}$")
    dc[U] = cross(s.delta.values, y.values)
ax[1].axhline(0, color="k", lw=0.9)
for i, U in enumerate(UU):
    if dc[U] is not None:
        ax[1].plot(dc[U], 0, "v", ms=8, color=CU[i], mec="k", mew=0.6, zorder=5)
ax[1].set_xlabel(r"anisotropy  $\delta$")
ax[1].set_ylabel(r"$\chi_{d_{xy}} - \chi_{d_{x^2-y^2}}$")
ax[1].legend(fontsize=10.5, loc="lower right")
ax[1].set_title(r"(b)  $d_{xy}$ overtakes $d_{x^2-y^2}$", loc="left", fontsize=14)

# ---------------- (c) margin over ALL other channels, in (U, delta) ----------
g = m[m.U > 0].copy()
g["marg"] = g.chi_dxy - g[["chi_son", "chi_sext", "chi_d"]].max(axis=1)
piv = g.pivot_table(index="U", columns="delta", values="marg")
X, Y = np.meshgrid(piv.columns.values, piv.index.values)
v = np.abs(piv.values).max()
# U = 0 is excluded: chi vanishes identically there, so it is a trivial zero row
# that would sit on the diverging colour map's midpoint and read as a boundary.
im = ax[2].pcolormesh(X, Y, piv.values, cmap="coolwarm", vmin=-v, vmax=v,
                      shading="nearest")
ax[2].contour(X, Y, piv.values, levels=[0], colors="k", linewidths=1.8)
cb = fig.colorbar(im, ax=ax[2], pad=0.02)
cb.set_label(r"$\chi_{d_{xy}} - \max_{\alpha \neq d_{xy}} \chi_\alpha$", fontsize=12)
ax[2].set_xlabel(r"anisotropy  $\delta$"); ax[2].set_ylabel(r"$U/t$")
ax[2].set_title(r"(c)  black line: $\delta^*$", loc="left", fontsize=14)

fig.tight_layout()
fig.savefig(f"{D}/figs/fig2_crossover.png", dpi=600, bbox_inches="tight",
            facecolor="white")

print("delta_d  (dxy overtakes dx2-y2):")
for U in UU:
    print(f"   U={U:g}: {dc[U]:.3f}" if dc[U] is not None else f"   U={U:g}: none")
print("\ndelta*   (dxy overtakes every other channel):")
for U in UU:
    s = m[m.U == U].sort_values("delta")
    x = cross(s.delta.values,
              (s.chi_dxy - s[["chi_son", "chi_sext", "chi_d"]].max(axis=1)).values)
    print(f"   U={U:g}: {x:.3f}" if x is not None else f"   U={U:g}: none")
plt.show()
