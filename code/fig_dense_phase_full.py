"""Dense-filling pairing phase diagram at full resolution, U = 4, L = 12.

19 fillings on a uniform 1/36 step, all eight anisotropies: 152 of 152 cells,
912 rows. Two things this fixes over the coarse version:

  - the coarse grid used 6 UNEVEN fillings (0.500, 0.556, then steps of 0.111),
    which over-weighted the low-filling corner. Cell counts and painted area
    disagreed there by up to a factor of four for d_xy. With uniform spacing
    they now agree.
  - the earlier dense run stopped at delta = 0.4, cutting off exactly where the
    d_xy region grows. It now runs to 0.7.

Colour is a blend of the LEADING TWO channels, weighted by their gap relative to
the local spread: pure where one leads clearly, mixed at boundaries, so the
transition width is itself the margin of dominance.
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
    "axes.labelsize": 19, "xtick.labelsize": 15, "ytick.labelsize": 15,
    "figure.dpi": 130})

D = "/Users/liujiaxin/Desktop/checkerboard-altermagnet/data"
# Vivid, high-saturation palette. The old ColorBrewer set was too dark: the
# top-two blend mutes whatever it mixes, so starting from muted colours left the
# whole plane looking grey. These stay distinct after blending, and the pairwise
# mixes (blue+green -> teal, blue+rose -> violet) are still readable.
CH = [("chi_son", "on-site $s$", "#8C8C8C"),      # grey, the control
      ("chi_sext", "extended $s$", "#1DB954"),    # teal
      ("chi_d", r"$d_{x^2-y^2}$", "#0047FF"),     # blue
      ("chi_dxy", r"$d_{xy}$", "#FF1500")]        # red
COLS = [c for c, _, _ in CH]
RGB = np.array([to_rgb(c) for _, _, c in CH])
BAND, BLUR = 0.14, 7.0   # narrower band -> more area at full saturation

d = pd.read_csv(f"{D}/dense_L12_U4_all.csv")
g = d.groupby(["n", "delta"])[COLS].mean().reset_index()
NS = np.sort(g.n.unique()); DS = np.sort(g.delta.unique())
print(f"grid: {len(NS)} fillings x {len(DS)} anisotropies = {len(NS)*len(DS)} cells, "
      f"{len(g)} present")

nf = np.linspace(NS.min(), NS.max(), 460); df = np.linspace(DS.min(), DS.max(), 460)
NG, DG = np.meshgrid(nf, df); pts = np.stack([DG.ravel(), NG.ravel()], -1)
F = np.stack([gaussian_filter(
    RegularGridInterpolator((DS, NS),
        g.pivot(index="delta", columns="n", values=c).values)(pts).reshape(DG.shape),
    BLUR) for c in COLS])

order = np.argsort(-F, axis=0); i1, i2 = order[0], order[1]
yy, xx = np.indices(i1.shape)
spread = np.maximum(F.max(0) - F.min(0), 1e-12)
frac = np.clip((F[i1, yy, xx] - F[i2, yy, xx]) / (BAND * spread), 0, 1)
img = (0.5 + 0.5*frac)[..., None]*RGB[i1] + (0.5*(1-frac))[..., None]*RGB[i2]

fig, ax = plt.subplots(figsize=(9.2, 6.6))
ax.imshow(img, origin="lower", aspect="auto",
          extent=[nf.min(), nf.max(), df.min(), df.max()], interpolation="bilinear")
ax.set_xlabel(r"Filling  $n$"); ax.set_ylabel(r"Anisotropy  $\delta$")
ax.set_xticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
ax.set_title(r"Dominant pairing channel, $U=4$, $L=12$" "\n"
             r"19 fillings $\times$ 8 anisotropies, uniform $\Delta n = 1/36$",
             fontsize=17, pad=12)
fig.legend(handles=[Patch(fc=c, ec="k", lw=0.6, label=l) for _, l, c in CH],
           fontsize=15, loc="center left", bbox_to_anchor=(0.90, 0.5),
           frameon=False, handlelength=1.8, handleheight=1.2, labelspacing=0.8)
fig.tight_layout(rect=[0, 0, 0.885, 1])
fig.savefig(f"{D}/../figs/fig_dense_phase_full.png", dpi=600,
            bbox_inches="tight", facecolor="white")

w = np.argmax(g[COLS].values, axis=1)
print("\ncell counts (now equal to area, since the spacing is uniform):")
for i in range(4):
    print(f"  {CH[i][1].replace('$',''):14s} {int((w==i).sum()):3d}/{len(g)}"
          f"  = {100*(w==i).mean():5.1f}%")
print("\nleading channel per (n, delta):")
piv = g.copy(); piv["win"] = [COLS[i].replace("chi_", "") for i in w]
t = piv.pivot(index="delta", columns="n", values="win")
print(t.iloc[::-1].to_string())
