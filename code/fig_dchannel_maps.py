"""The two d-wave channels as continuous (n, delta) maps, in both time sectors.

The companion figures ask WHICH channel leads. This one asks HOW STRONG the
d-wave channels are, plotted as a continuous field the way the PRL shows
N^Vertex for a single channel.

Rows are the two time sectors, columns the two d-wave symmetries. Both rows are
divided so they carry the same units: chi alone has an extra factor of imaginary
time and is not comparable to C(0).

Each panel keeps its own colour scale. The four quantities differ by more than an
order of magnitude, so a shared scale would flatten three of them -- which also
means colours must not be compared across panels, only the numbers.
"""
import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import gaussian_filter

mpl.rcParams.update({
    "font.family": "serif", "font.serif": ["DejaVu Serif"],
    "mathtext.fontset": "dejavuserif", "axes.linewidth": 1.2,
    "xtick.direction": "out", "ytick.direction": "out",
    # ticks on the bottom and left only: with contourf the top/right
    # marks sit on the filled area and read as part of the data
    "xtick.top": False, "ytick.right": False,
    "axes.labelsize": 17, "xtick.labelsize": 13, "ytick.labelsize": 13,
    "figure.dpi": 130})

D = "/Users/liujiaxin/Desktop/checkerboard-altermagnet/data"
NSITES, TAU_MAX = 100, 16 * 0.05    # L = 10;  MUST match BP, DT in the driver
BLUR = 4.0
CHAN = [("d", r"$d_{x^2-y^2}$"), ("dxy", r"$d_{xy}$")]
SECT = [("eq", NSITES, "equal time"), ("chi", NSITES * TAU_MAX, "unequal time")]

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


fig, ax = plt.subplots(2, 2, figsize=(13.2, 9.6), sharex=True, sharey=True)
for r, (pre, sc, sect) in enumerate(SECT):
    for c, (key, sym) in enumerate(CHAN):
        a = ax[r, c]
        F = smooth_field(f"{pre}_{key}_vertex", sc)
        im = a.contourf(nf, df, F, levels=100, cmap="jet", extend="both")
        cb = fig.colorbar(im, ax=a, pad=0.02, shrink=0.92)
        cb.set_ticks(np.linspace(F.min(), F.max(), 5))
        cb.ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter("%.3f"))
        cb.ax.tick_params(labelsize=11)
        a.set_title(f"({chr(97 + 2*r + c)}) {sym}, {sect}", fontsize=16)
        if r == 1: a.set_xlabel(r"Filling  $n$")
        if c == 0: a.set_ylabel(r"Anisotropy  $\delta$")
        a.set_xticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
        print(f"({chr(97 + 2*r + c)}) {key:3s} {sect:13s} "
              f"{F.min():+9.4f} to {F.max():+9.4f}")
fig.tight_layout()
fig.savefig(f"{D}/../figs/fig_dchannel_maps.png", dpi=600, bbox_inches="tight",
            facecolor="white")
