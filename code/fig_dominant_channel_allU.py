"""Dominant pairing channel in the (n, delta) plane, one panel per U.

Uses the total susceptibility chi_alpha from pairing_master.csv. No R-resolved
decomposition here: the question is simply which channel leads at each point.

Smooth colour without pretending the quantity is continuous. Each channel owns an
RGB colour, and every point on the fine grid is painted with a softmax-weighted
BLEND of the four. Where one channel leads clearly its colour comes through pure;
near a boundary the two colours mix, so the transition is continuous instead of a
hard edge, and the width of the mixed band is itself the margin of dominance.

The softmax temperature is set to a fraction of the local spread across channels,
so the sharpness is the same at U = 2 and U = 8 even though chi shrinks by an
order of magnitude between them.
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

CSV = "/Users/liujiaxin/Desktop/checkerboard-altermagnet/data/pairing_master.csv"
L = 12
CH = [("chi_son", r"on-site $s$", "#8C8C8C"),
      ("chi_sext", r"extended $s$", "#1DB954"),
      ("chi_d", r"$d_{x^2-y^2}$", "#0047FF"),
      ("chi_dxy", r"$d_{xy}$", "#FF1500")]
COLS = [c for c, _, _ in CH]
RGB = np.array([to_rgb(c) for _, _, c in CH])
BAND = 0.14          # width of the blended boundary band, as a fraction of spread
BLUR = 5.0           # grid steps of smoothing after interpolation

d = pd.read_csv(CSV); d = d[d.L == L]
g = d.groupby(["U", "n", "delta"])[COLS].mean().reset_index()
US = [2.0, 4.0, 6.0, 8.0]
NS = np.sort(g.n.unique()); DS = np.sort(g.delta.unique())
nf = np.linspace(NS.min(), NS.max(), 420); df = np.linspace(DS.min(), DS.max(), 420)
NG, DG = np.meshgrid(nf, df); pts = np.stack([DG.ravel(), NG.ravel()], -1)


def field(U):
    """Each channel interpolated onto the fine grid, shape (4, ny, nx)."""
    s = g[np.isclose(g.U, U)]
    out = []
    for c in COLS:
        M = s.pivot(index="delta", columns="n", values=c).values
        F = RegularGridInterpolator((DS, NS), M)(pts).reshape(DG.shape)
        out.append(gaussian_filter(F, BLUR))
    return np.stack(out)


fig, ax = plt.subplots(1, 4, figsize=(22.2, 5.4), sharey=True)
for k, U in enumerate(US):
    a = ax[k]; F = field(U)
    # Blend only the LEADING TWO channels. A softmax over all four mixed three
    # colours at once and produced grey-brown everywhere, which hid the regions.
    order = np.argsort(-F, axis=0)
    i1, i2 = order[0], order[1]
    yy, xx = np.indices(i1.shape)
    top1, top2 = F[i1, yy, xx], F[i2, yy, xx]
    spread = np.maximum(F.max(0) - F.min(0), 1e-12)
    f = np.clip((top1 - top2) / (BAND * spread), 0, 1)   # 0 = tie, 1 = clear win
    c1, c2 = RGB[i1], RGB[i2]
    img = (0.5 + 0.5 * f)[..., None] * c1 + (0.5 * (1 - f))[..., None] * c2
    a.imshow(img, origin="lower", aspect="auto",
             extent=[nf.min(), nf.max(), df.min(), df.max()],
             interpolation="bilinear")
    a.set_title(rf"$U = {U:g}$", fontsize=24)
    a.set_xlabel(r"Filling  $n$")
    a.set_xticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
ax[0].set_ylabel(r"Anisotropy  $\delta$")
fig.legend(handles=[Patch(fc=c, ec="k", lw=0.6, label=l) for _, l, c in CH],
           fontsize=18, loc="center left", ncol=1, bbox_to_anchor=(0.905, 0.5),
           frameon=False, handlelength=1.9, handleheight=1.25, labelspacing=0.85)
fig.suptitle(rf"Dominant pairing channel, $L={L}$   "
             r"(colours blend where channels are close)", fontsize=22, y=1.02)
fig.tight_layout(rect=[0, 0, 0.895, 1])
fig.savefig("fig_dominant_channel_allU.png", dpi=600, bbox_inches="tight",
            facecolor="white")

print("winner counts per U (on the measured 6 x 8 grid):")
for U in US:
    s = g[np.isclose(g.U, U)]
    wv = np.argmax(s[COLS].values, axis=1)
    print(f"  U={U:g}: " + ", ".join(
        f"{CH[i][1].replace('$','')} {int((wv==i).sum())}" for i in range(4)))
