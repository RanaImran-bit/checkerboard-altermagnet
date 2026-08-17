"""Equal-time vs unequal-time pairing, BOTH sectors from our own code.

The earlier version took equal time from the group's Fortran and unequal time
from ours. Their vertex is full - L*Unpair, ours is full - bubble, so the
magnitudes were not comparable and only the channel ORDERING could be read.

checkerboard_eqtime.py now saves both from the SAME driver at the SAME
parameters -- C_a(tau=0) alongside chi_a = int C_a(tau) dtau -- so that caveat
is gone and the two sectors can be compared directly, in the same units.

L = 10, U = 4, 11 fillings n = 0.50 to 1.00, delta = 0 to 0.4, 6 seeds.
"""
import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
from matplotlib.patches import Patch
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import gaussian_filter

mpl.rcParams.update({
    "font.family": "serif", "font.serif": ["DejaVu Serif"],
    "mathtext.fontset": "dejavuserif", "axes.linewidth": 1.1,
    "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True,
    "ytick.right": True, "legend.frameon": False,
    "axes.labelsize": 14, "xtick.labelsize": 11, "ytick.labelsize": 11,
    "figure.dpi": 130})

D = "/Users/liujiaxin/Desktop/checkerboard-altermagnet/data"
CH = [("son", "on-site $s$", "#8C8C8C", "D"),
      ("sext", "extended $s$", "#1DB954", "s"),
      ("d", r"$d_{x^2-y^2}$", "#0047FF", "o"),
      ("dxy", r"$d_{xy}$", "#FF1500", "^")]
KEYS = [c for c, *_ in CH]
RGB = np.array([to_rgb(x[2]) for x in CH])   # x = (key, label, colour, marker)
NSITES = 100
TAU_MAX = 16 * 0.05     # BP = 16 steps at dt = 0.05
BAND, BLUR = 0.14, 5.5

# UNITS. C_a(0) is a correlation. chi_a = int C_a(tau) dtau carries an extra
# factor of imaginary time, so C(0)/N and chi/N are NOT the same quantity and
# must not share an axis. Dividing chi by the integration window gives the
# TIME-AVERAGED correlation, <C_a>_tau, which has the same units as C_a(0) and
# is directly comparable to it. The phase diagrams below are unaffected either
# way, since argmax is invariant under a positive rescaling.

d = pd.read_csv(f"{D}/eqtime_L10_U4_all.csv")
DL = sorted(d.delta.unique())
EQ = d.groupby(["n", "delta"])[[f"eq_{c}_vertex" for c in KEYS]].mean()
UN = d.groupby(["n", "delta"])[[f"chi_{c}_vertex" for c in KEYS]].mean()
EQ.columns = KEYS; UN.columns = KEYS

# ---------- (1) channel curves, both sectors, same units ----------
fig, ax = plt.subplots(2, len(DL), figsize=(4.0 * len(DL), 8.0), sharex=True)
for j, dl in enumerate(DL):
    for row, (T, lab, sc) in enumerate([
            (EQ, r"equal time" "\n" r"$C_\alpha(0)/N$", NSITES),
            (UN, r"time averaged" "\n" r"$\langle C_\alpha\rangle_\tau/N$",
             NSITES * TAU_MAX)]):
        a = ax[row, j]
        s = T.xs(dl, level="delta").sort_index()
        for c, lb, col, mk in CH:
            a.plot(s.index, s[c] / sc, marker=mk, ms=5, lw=1.4,
                   color=col, label=lb)
        a.axhline(0, color="k", lw=1.0, ls="--")
        a.set_xlim(0.48, 1.02)
        if row == 0: a.set_title(rf"$\delta = {dl:g}$", fontsize=14)
        if row == 1: a.set_xlabel(r"filling  $n$")
        if j == 0: a.set_ylabel(lab)
ax[0, 0].legend(fontsize=9, loc="best")
fig.tight_layout()
fig.savefig(f"{D}/../figs/fig_eq_vs_uneq_onecode.png", dpi=600,
            bbox_inches="tight", facecolor="white")

# ---------- (2) dominant-channel phase diagram, both sectors ----------
NS = np.sort(d.n.unique()); DS = np.array(DL)
nf = np.linspace(NS.min(), NS.max(), 380); df = np.linspace(DS.min(), DS.max(), 380)
NG, DG = np.meshgrid(nf, df); pts = np.stack([DG.ravel(), NG.ravel()], -1)

fig, ax = plt.subplots(1, 2, figsize=(14.2, 5.8), sharey=True)
for k, (T, lab) in enumerate([(EQ, r"(a) equal time   $\tau=0$"),
                              (UN, r"(b) unequal time   $\int d\tau$")]):
    a = ax[k]
    F = np.stack([gaussian_filter(
        RegularGridInterpolator((DS, NS),
            T[c].unstack("n").values)(pts).reshape(DG.shape), BLUR) for c in KEYS])
    order = np.argsort(-F, axis=0); i1, i2 = order[0], order[1]
    yy, xx = np.indices(i1.shape)
    spread = np.maximum(F.max(0) - F.min(0), 1e-12)
    frac = np.clip((F[i1, yy, xx] - F[i2, yy, xx]) / (BAND * spread), 0, 1)
    img = (0.5 + 0.5*frac)[..., None]*RGB[i1] + (0.5*(1-frac))[..., None]*RGB[i2]
    a.imshow(img, origin="lower", aspect="auto",
             extent=[nf.min(), nf.max(), df.min(), df.max()], interpolation="bilinear")
    a.set_xlabel(r"Filling  $n$"); a.set_title(lab, fontsize=17)
    a.set_xticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
ax[0].set_ylabel(r"Anisotropy  $\delta$")
fig.legend(handles=[Patch(fc=c, ec="k", lw=0.6, label=l) for _, l, c, _ in CH],
           fontsize=14, loc="center left", bbox_to_anchor=(0.905, 0.5),
           frameon=False, handlelength=1.8, handleheight=1.2, labelspacing=0.8)
fig.tight_layout(rect=[0, 0, 0.90, 1])
fig.savefig(f"{D}/../figs/fig_eq_vs_uneq_phase_onecode.png", dpi=600,
            bbox_inches="tight", facecolor="white")

for T, lab in [(EQ, "EQUAL TIME"), (UN, "UNEQUAL TIME")]:
    w = T[KEYS].values.argmax(axis=1)
    print(f"{lab:13s}: " + ", ".join(
        f"{CH[i][1].replace('$',''):12s}{int((w==i).sum()):3d}/{len(T)}" for i in range(4)))
print("\nleading channel, EQUAL TIME (rows delta, cols n):")
print(pd.Series([KEYS[i] for i in EQ[KEYS].values.argmax(1)], index=EQ.index)
      .unstack("n").iloc[::-1].to_string())
print("\nleading channel, UNEQUAL TIME:")
print(pd.Series([KEYS[i] for i in UN[KEYS].values.argmax(1)], index=UN.index)
      .unstack("n").iloc[::-1].to_string())
