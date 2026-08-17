"""Which d-wave channel leads, one panel per time sector.

Only the two d-wave symmetries are compared here, so the map needs only two
colours: blue where d_x2-y2 leads, red where d_xy leads, blending through
violet where they are close. Same blend machinery as the four-channel phase
diagrams, so the figures sit together consistently.

Both sectors are divided so they carry the same units: chi alone has an extra
factor of imaginary time and is not comparable to C(0). The blend and the
boundary are unaffected by that scaling, since a positive rescaling cannot
change which channel is larger.
"""
import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
from matplotlib.patches import Patch
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import gaussian_filter

mpl.rcParams.update({
    "font.family": "serif", "font.serif": ["DejaVu Serif"],
    "mathtext.fontset": "dejavuserif", "axes.linewidth": 1.2,
    "xtick.direction": "out", "ytick.direction": "out",
    "xtick.top": False, "ytick.right": False,
    "axes.labelsize": 17, "xtick.labelsize": 13, "ytick.labelsize": 13,
    "figure.dpi": 130})

D = "/Users/liujiaxin/Desktop/checkerboard-altermagnet/data"
NSITES, TAU_MAX = 100, 16 * 0.05    # L = 10;  MUST match BP, DT in the driver
BAND, BLUR = 0.14, 4.0
CH = [("d", r"$d_{x^2-y^2}$", "#0047FF"),      # blue
      ("dxy", r"$d_{xy}$", "#FF1500")]         # red
RGB = np.array([to_rgb(x[2]) for x in CH])
SECT = [("eq", NSITES, r"(a) equal time   $\tau=0$"),
        ("chi", NSITES * TAU_MAX, r"(b) unequal time   $\int d\tau$")]

d = pd.read_csv(f"{D}/eqtime_L10_U4_all.csv")
NS = np.sort(d.n.unique()); DS = np.sort(d.delta.unique())
nf = np.linspace(NS.min(), NS.max(), 400); df = np.linspace(DS.min(), DS.max(), 400)
NG, DG = np.meshgrid(nf, df); pts = np.stack([DG.ravel(), NG.ravel()], -1)


def smooth_field(col, scale):
    g = d.groupby(["n", "delta"])[col].mean()
    assert len(g) == len(NS) * len(DS), "grid has holes"
    M = g.unstack("n").values          # (len(DS), len(NS)), the order the
    F = RegularGridInterpolator((DS, NS), M)(pts).reshape(DG.shape)  # interpolator wants
    return gaussian_filter(F, BLUR) / scale


fig, ax = plt.subplots(1, 2, figsize=(13.6, 5.8), sharey=True)
for k, (pre, sc, ttl) in enumerate(SECT):
    a = ax[k]
    F = np.stack([smooth_field(f"{pre}_{key}_vertex", sc) for key, *_ in CH])
    lead = F[1] - F[0]                                  # >0 where d_xy leads
    # colour saturates once the gap exceeds BAND times the local spread, and
    # blends toward violet inside it, so the transition width is the margin
    spread = np.maximum(np.abs(F[0] - F[1]).max(), 1e-12)
    frac = np.clip(np.abs(lead) / (BAND * spread), 0, 1)
    win = (lead > 0).astype(int)
    c1, c2 = RGB[win], RGB[1 - win]
    img = (0.5 + 0.5 * frac)[..., None] * c1 + (0.5 * (1 - frac))[..., None] * c2
    a.imshow(img, origin="lower", aspect="auto",
             extent=[nf.min(), nf.max(), df.min(), df.max()],
             interpolation="bilinear")
    a.set_title(ttl, fontsize=17)
    a.set_xlabel(r"Filling  $n$")
    a.set_xticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    if k == 0: a.set_ylabel(r"Anisotropy  $\delta$")
    print(f"{ttl:34s} d_xy leads over {100*(lead>0).mean():.0f}% of the plane")
fig.legend(handles=[Patch(fc=c, ec="k", lw=0.6, label=l) for _, l, c in CH],
           fontsize=15, loc="center left", bbox_to_anchor=(0.905, 0.5),
           frameon=False, handlelength=1.9, handleheight=1.25, labelspacing=0.9)
fig.tight_layout(rect=[0, 0, 0.90, 1])
fig.savefig(f"{D}/../figs/fig_dchannel_diff.png", dpi=600, bbox_inches="tight",
            facecolor="white")
