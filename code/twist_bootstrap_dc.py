"""delta_c with a bootstrap error bar, from the twist-averaged fine grid.

Combines the coarse twist grid (delta = 0 to 0.7, step 0.1) with the fine grid
(U=2: 0.22-0.28, U=4: 0.12-0.18) so the zero crossing of the d_xy vertex is
bracketed at step 0.02 instead of 0.1.

Twist average: only theta = 0 and pi keep the hopping real, which the engine
requires. The two mixed twists are degenerate by the lattice's x <-> y symmetry,
so three runs carry weights 1 : 2 : 1.

    (+1,+1) periodic          weight 1
    (-1,+1) half twist        weight 2   (stands for (+1,-1) as well)
    (-1,-1) both antiperiodic weight 1

Bootstrap: resample the 6 seeds with replacement independently for each
(delta, leg), rebuild the weighted average, and locate the zero crossing of
chi_dxy by linear interpolation between the two bracketing deltas. The spread
over resamples is the error on delta_c.
"""
import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt

mpl.rcParams.update({
    "font.family": "serif", "font.serif": ["DejaVu Serif"],
    "mathtext.fontset": "dejavuserif", "axes.linewidth": 1.1,
    "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True,
    "ytick.right": True, "legend.frameon": False,
    "axes.labelsize": 14, "xtick.labelsize": 12, "ytick.labelsize": 12,
    "figure.dpi": 130})

D = "/Users/liujiaxin/Desktop/checkerboard-altermagnet/data"
SCR = ("/private/tmp/claude-501/-Users-liujiaxin-Desktop-Susceptibility-"
       "qmc-platform-master/851bd175-76db-4789-ac47-3b3fb836a9b4/scratchpad/tw3")
W = {(1, 1): 1.0, (-1, 1): 2.0, (-1, -1): 1.0}
NB = 4000
rng = np.random.default_rng(0)

fine = pd.read_csv(f"{SCR}/tw3_merged.csv")
coarse = pd.read_csv(f"{D}/twist_merged.csv")
COL = "chi_dxy_vertex"
keep = ["U", "delta", "apx", "apy", "seed", COL]
d = pd.concat([coarse[keep], fine[keep]], ignore_index=True).drop_duplicates()


def weighted_curve(sub):
    """Weighted twist average per delta, from a table already reduced to one
    value per (delta, leg)."""
    out = {}
    for dl, g in sub.groupby("delta"):
        num = den = 0.0
        for (ax, ay), gg in g.groupby(["apx", "apy"]):
            w = W[(ax, ay)]
            num += w * gg[COL].mean(); den += w
        out[dl] = num / den
    return pd.Series(out).sort_index()


def crossing(s):
    """First sign change of the curve, linearly interpolated."""
    x, y = s.index.values.astype(float), s.values
    i = np.where(np.diff(np.sign(y)))[0]
    if not len(i):
        return np.nan
    k = i[0]
    return x[k] - y[k] * (x[k + 1] - x[k]) / (y[k + 1] - y[k])


fig, ax = plt.subplots(1, 2, figsize=(12.4, 4.8))
print("delta_c from the twist-averaged d_xy vertex, L=12, half filling\n")
for j, U in enumerate([2.0, 4.0]):
    sub = d[d.U == U]
    curve = weighted_curve(sub)
    dc = crossing(curve)
    # bootstrap over seeds, independently per (delta, leg)
    boots = []
    groups = {k: g[COL].values for k, g in sub.groupby(["delta", "apx", "apy"])}
    for _ in range(NB):
        rows = []
        for (dl, ax_, ay_), v in groups.items():
            rows.append((dl, ax_, ay_, rng.choice(v, size=len(v)).mean()))
        b = pd.DataFrame(rows, columns=["delta", "apx", "apy", COL])
        boots.append(crossing(weighted_curve(b)))
    boots = np.array([x for x in boots if np.isfinite(x)])
    lo, hi = np.percentile(boots, [16, 84])
    print(f"U = {U:g}:  delta_c = {dc:.4f}  +{hi-dc:.4f} -{dc-lo:.4f}   "
          f"(bootstrap 68% CI [{lo:.4f}, {hi:.4f}], {len(boots)}/{NB} valid)")

    a = ax[j]
    err = sub.groupby("delta")[COL].sem()
    a.errorbar(curve.index, curve.values, yerr=err.reindex(curve.index).values,
               marker="o", ms=5, lw=1.4, color="k", capsize=2.5)
    a.axhline(0, color="k", lw=1.0, ls="--")
    a.axvspan(lo, hi, color="#b2182b", alpha=0.20)
    a.axvline(dc, color="#b2182b", lw=1.8)
    a.set_xlabel(r"anisotropy  $\delta$")
    a.set_ylabel(r"twist-averaged $\chi_{d_{xy}}$ vertex")
    a.set_title(rf"$U = {U:g}$:  $\delta_c = {dc:.3f}^{{+{hi-dc:.3f}}}_{{-{dc-lo:.3f}}}$",
                fontsize=14)
fig.tight_layout()
fig.savefig(f"{D}/../figs/fig_twist_deltac.png", dpi=600, bbox_inches="tight",
            facecolor="white")
