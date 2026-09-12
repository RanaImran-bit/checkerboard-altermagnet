"""Publication figures for the real-space pairing vertex V(R).

Follows Huang, Lin and Gubernatis PRB 64, 205101: split the pair-field vertex by
the separation R between the two pairs, then compare a LOCAL piece against a
LONG-RANGE piece. The driver stores a SUM over all site pairs at each separation,
and the multiplicity runs from 144 to 1728, so everything here is divided by the
pair count to give a per-pair correlation v(R) = V(R) / npairs(R). Without that
division the distant shells look artificially large.
"""
import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.lines import Line2D

mpl.rcParams.update({
    "font.family": "serif", "font.serif": ["DejaVu Serif"], "mathtext.fontset": "dejavuserif",
    "axes.linewidth": 1.0, "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True, "xtick.major.size": 5,
    "ytick.major.size": 5, "xtick.minor.size": 2.8, "ytick.minor.size": 2.8,
    "axes.labelsize": 13, "xtick.labelsize": 11, "ytick.labelsize": 11,
    "legend.frameon": False, "figure.dpi": 130,
})

CH = [("sext", r"extended $s$", "#1b7837", "s"),
      ("d",    r"$d_{x^2-y^2}$", "#2166ac", "o"),
      ("dxy",  r"$d_{xy}$",      "#b2182b", "^"),
      ("son",  r"on-site $s$",   "#7a7a7a", "D")]
LR = 2.0        # long range means R > LR, as in Huang

L = 12; n = L * L
xs, ys = np.arange(n) // L, np.arange(n) % L
dx = np.abs(xs[:, None] - xs[None, :]); dx = np.minimum(dx, L - dx)
dy = np.abs(ys[:, None] - ys[None, :]); dy = np.minimum(dy, L - dy)
Rm = np.round(np.sqrt(dx**2 + dy**2), 6)
rv, rc = np.unique(Rm, return_counts=True)
mult = dict(zip(rv, rc))

d = pd.read_csv("pr_all.csv")
d["npair"] = d.R.round(6).map(mult)
d["v"] = d.V / d.npair


def prof(U, delta, nup=72):
    """v(R) with a seed error bar, one row per (channel, R)."""
    g = d[(d.U == U) & (d.delta == delta) & (d.nup == nup)]
    a = g.groupby(["channel", "R"]).v.agg(["mean", "sem"]).reset_index()
    return a


def local_long(U, delta, nup=72):
    """v(R=0) and the pair-count-weighted mean of v over R > LR."""
    g = d[(d.U == U) & (d.delta == delta) & (d.nup == nup)]
    out = {}
    for c in g.channel.unique():
        s = g[g.channel == c]
        loc = s[s.R == 0].groupby("seed").v.mean()
        f = s[s.R > LR].groupby("seed").apply(
            lambda z: z.V.sum() / z.npair.sum(), include_groups=False)
        out[c] = (loc.mean(), loc.sem(), f.mean(), f.sem())
    return out


# ============================ FIGURE 1: v(R) profiles ======================
DL = [0.0, 0.2, 0.4, 0.7]
fig, ax = plt.subplots(1, 4, figsize=(15.5, 4.1), sharey=True)
for j, dl in enumerate(DL):
    a = ax[j]; t = prof(4.0, dl)
    for c, lab, col, mk in CH:
        s = t[t.channel == c].sort_values("R")
        y, e, R = s["mean"].values, s["sem"].values, s.R.values
        pos = y > 0
        a.errorbar(R[pos], np.abs(y[pos]), yerr=e[pos], color=col, marker=mk, ms=5.5,
                   lw=1.3, ls="-", mfc=col, mec=col, capsize=0, label=lab, zorder=3)
        a.errorbar(R[~pos], np.abs(y[~pos]), yerr=e[~pos], color=col, marker=mk, ms=5.5,
                   lw=1.3, ls=":", mfc="white", mec=col, capsize=0, zorder=3)
    a.set_yscale("log"); a.set_xlabel(r"$R$")
    a.set_title(rf"$\delta = {dl:g}$", fontsize=13)
    a.set_xlim(-0.4, 9); a.set_ylim(1e-6, 3e-2)
    a.axvline(LR, color="0.8", lw=1.0, ls="--", zorder=0)
ax[0].set_ylabel(r"$|v(R)|$   per-pair vertex")
ax[0].legend(fontsize=10.5, loc="upper right", ncol=1)
ex = [Line2D([], [], color="k", marker="o", ls="-", mfc="k", label="attractive  $v>0$"),
      Line2D([], [], color="k", marker="o", ls=":", mfc="white", label="repulsive  $v<0$")]
ax[1].legend(handles=ex, fontsize=10.5, loc="upper right")
fig.suptitle(r"Real-space pairing vertex $v(R)$ at $U=4$, half filling, $L=12$"
             "\n" r"dashed line marks the local / long-range split at $R=2$",
             fontsize=13.5, y=1.06)
fig.tight_layout(); fig.savefig("fig_pr_profiles.png", dpi=300, bbox_inches="tight",
                                facecolor="white")

# ================= FIGURE 2: local vs long range against delta ============
DLS = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
fig, ax = plt.subplots(1, 2, figsize=(11.6, 4.3))
tab = {dl: local_long(4.0, dl) for dl in DLS}
for k, (idx, ttl) in enumerate([(0, "local  " + r"$v(R=0)$"),
                                (2, "long range  " + r"$\langle v\rangle_{R>2}$")]):
    a = ax[k]
    for c, lab, col, mk in CH:
        y = [tab[dl][c][idx] * 1e3 for dl in DLS]
        e = [tab[dl][c][idx + 1] * 1e3 for dl in DLS]
        a.errorbar(DLS, y, yerr=e, color=col, marker=mk, ms=6.5, lw=1.7,
                   capsize=3, label=lab)
    a.axhline(0, color="k", lw=1.0, ls="--")
    a.set_xlabel(r"anisotropy  $\delta$"); a.set_title(ttl, fontsize=13)
    a.set_ylabel(r"$v \times 10^{3}$")
ax[0].legend(fontsize=10.5, loc="lower left")
fig.suptitle(r"Local and long-range pairing vertex at $U=4$, half filling",
             fontsize=13.5, y=1.02)
fig.tight_layout(); fig.savefig("fig_pr_local_long.png", dpi=300, bbox_inches="tight",
                                facecolor="white")

# ============ FIGURE 3: long-range channel against delta, per U ===========
US = [2.0, 4.0, 6.0, 8.0]
fig, ax = plt.subplots(1, 4, figsize=(15.5, 4.0), sharey=True)
for j, U in enumerate(US):
    a = ax[j]; t = {dl: local_long(U, dl) for dl in DLS}
    for c, lab, col, mk in CH:
        y = [t[dl][c][2] * 1e3 for dl in DLS]
        e = [t[dl][c][3] * 1e3 for dl in DLS]
        a.errorbar(DLS, y, yerr=e, color=col, marker=mk, ms=6, lw=1.6,
                   capsize=3, label=lab)
    a.axhline(0, color="k", lw=1.0, ls="--")
    a.set_xlabel(r"$\delta$"); a.set_title(rf"$U = {U:g}$", fontsize=13)
ax[0].set_ylabel(r"long-range $\langle v\rangle_{R>2}\times 10^{3}$")
ax[0].legend(fontsize=10, loc="upper left")
fig.suptitle(r"Long-range pairing vertex against anisotropy, half filling, $L=12$",
             fontsize=13.5, y=1.04)
fig.tight_layout(); fig.savefig("fig_pr_longrange_U.png", dpi=300, bbox_inches="tight",
                                facecolor="white")

# ---- numbers quoted in the write-up ----
print("U=4 half filling, per-pair vertex x1000")
print(f"{'delta':>6}  {'LOCAL: sext':>12}{'d':>9}{'dxy':>9}   {'LONG: sext':>12}{'d':>9}{'dxy':>9}")
for dl in DLS:
    t = tab[dl]
    print(f"{dl:6.1f}  {t['sext'][0]*1e3:12.2f}{t['d'][0]*1e3:9.2f}{t['dxy'][0]*1e3:9.2f}"
          f"   {t['sext'][2]*1e3:12.3f}{t['d'][2]*1e3:9.3f}{t['dxy'][2]*1e3:9.3f}")
print()
print("long-range winner per (U, delta):")
for U in US:
    t = {dl: local_long(U, dl) for dl in DLS}
    w = [max(["sext", "d", "dxy"], key=lambda c: t[dl][c][2]) for dl in DLS]
    print(f"  U={U:g}: " + "  ".join(f"{dl:g}:{x}" for dl, x in zip(DLS, w)))
