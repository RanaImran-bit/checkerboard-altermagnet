"""Dense-filling dominant-channel panel at U = 4.

19 fillings on a uniform 1/36 step against the 6 uneven fillings of the coarse
grid, so the boundary positions here are measured rather than interpolated across
0.11-wide gaps. Uniform spacing also means cell counts and painted area finally
agree, which they did not on the coarse grid.

Only delta <= 0.4 exists in the dense set, so the panel stops there.
Same top-two colour blend as the coarse figure.
"""
import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
from matplotlib.patches import Patch
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import gaussian_filter

mpl.rcParams.update({
    "font.family": "serif", "font.serif": ["DejaVu Serif"],
    "mathtext.fontset": "dejavuserif", "axes.linewidth": 1.2,
    "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True,
    "ytick.right": True, "xtick.major.size": 6, "ytick.major.size": 6,
    "axes.labelsize": 20, "xtick.labelsize": 16, "ytick.labelsize": 16,
    "figure.dpi": 130})

D = "/Users/liujiaxin/Desktop/checkerboard-altermagnet/data"
CH = [("chi_son", r"on-site $s$", "#8c8c8c"),
      ("chi_sext", r"extended $s$", "#1b7837"),
      ("chi_d", r"$d_{x^2-y^2}$", "#2166ac"),
      ("chi_dxy", r"$d_{xy}$", "#b2182b")]
COLS = [c for c, _, _ in CH]; RGB = np.array([to_rgb(c) for _, _, c in CH])
BAND, BLUR = 0.22, 3.0


def panel(a, g, ttl, blur=BLUR):
    NS = np.sort(g.n.unique()); DS = np.sort(g.delta.unique())
    nf = np.linspace(NS.min(), NS.max(), 460)
    df = np.linspace(DS.min(), DS.max(), 460)
    NG, DG = np.meshgrid(nf, df); pts = np.stack([DG.ravel(), NG.ravel()], -1)
    F = np.stack([gaussian_filter(RegularGridInterpolator((DS, NS),
        g.pivot(index="delta", columns="n", values=c).values)(pts).reshape(DG.shape),
        blur) for c in COLS])
    order = np.argsort(-F, axis=0); i1, i2 = order[0], order[1]
    yy, xx = np.indices(i1.shape)
    spread = np.maximum(F.max(0) - F.min(0), 1e-12)
    f = np.clip((F[i1, yy, xx] - F[i2, yy, xx]) / (BAND * spread), 0, 1)
    img = (0.5 + 0.5 * f)[..., None] * RGB[i1] + (0.5 * (1 - f))[..., None] * RGB[i2]
    a.imshow(img, origin="lower", aspect="auto",
             extent=[nf.min(), nf.max(), df.min(), df.max()],
             interpolation="bilinear")
    a.set_title(ttl, fontsize=21)
    a.set_xlabel(r"Filling  $n$"); a.set_xticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    return F, NS, DS


dn = pd.read_csv(f"{D}/dense_L12_U4_all.csv")
gd = dn.groupby(["n", "delta"])[COLS].mean().reset_index()
co = pd.read_csv(f"{D}/pairing_master.csv")
co = co[(co.L == 12) & (co.U == 4) & (co.delta <= 0.4)]
gc = co.groupby(["n", "delta"])[COLS].mean().reset_index()

fig, ax = plt.subplots(1, 2, figsize=(15.2, 5.6), sharey=True)
panel(ax[1], gd, r"dense: 19 fillings, uniform $1/36$ step")
panel(ax[0], gc, r"coarse: 6 fillings, uneven spacing")
ax[0].set_ylabel(r"Anisotropy  $\delta$")
fig.legend(handles=[Patch(fc=c, ec="k", lw=0.6, label=l) for _, l, c in CH],
           fontsize=17, loc="center left", ncol=1, bbox_to_anchor=(0.875, 0.5),
           frameon=False, handlelength=1.9, handleheight=1.25, labelspacing=0.85)
fig.suptitle(r"Dominant pairing channel at $U=4$, $L=12$   "
             r"(same data source, different filling resolution)", fontsize=20, y=1.0)
fig.tight_layout(rect=[0, 0, 0.865, 1])
fig.savefig("fig_dense_vs_coarse_U4.png", dpi=600, bbox_inches="tight",
            facecolor="white")

# does the extra resolution change the answer?
for tag, g in [("coarse", gc), ("dense", gd)]:
    w = np.argmax(g[COLS].values, axis=1)
    tot = len(g)
    print(f"{tag:>7} ({tot:3d} cells): " + "  ".join(
        f"{CH[i][1].replace('$','').replace('\\','')} {100*(w==i).mean():4.1f}%"
        for i in range(4) if (w == i).any()))
print()
print("d_xy region boundary at each delta (dense grid): lowest n where d_xy leads")
for dl in sorted(gd.delta.unique()):
    s = gd[gd.delta == dl].sort_values("n")
    w = np.argmax(s[COLS].values, axis=1)
    ns = s.n.values[w == 3]
    print(f"   delta={dl:.1f}: " + (f"n >= {ns.min():.3f}  ({len(ns)} of {len(s)} fillings)"
                                    if len(ns) else "d_xy never leads"))
