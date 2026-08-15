"""(n, delta) dominant-channel phase diagram in BOTH time sectors, L=10, U=4.

  left   EQUAL TIME    Fortran CPQMC, Vertex_*.dat at k=(0,0), tau=0
  right  UNEQUAL TIME  our Python CPQMC, chi = int C(tau) dtau

Same top-two colour blend as the main phase diagram: pure colour where one
channel leads clearly, blended where two compete, so boundaries merge smoothly.

The Fortran grid is irregular -- 27 runs, not a full 10 x 5 product -- so its
panel is interpolated from scattered points and the measured points are marked.
"""
import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
from matplotlib.patches import Patch
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter

mpl.rcParams.update({
    "font.family": "serif", "font.serif": ["DejaVu Serif"],
    "mathtext.fontset": "dejavuserif", "axes.linewidth": 1.2,
    "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True,
    "ytick.right": True, "axes.labelsize": 17, "xtick.labelsize": 14,
    "ytick.labelsize": 14, "figure.dpi": 130})

D = "/Users/liujiaxin/Desktop/checkerboard-altermagnet/data"
SCR = ("/private/tmp/claude-501/-Users-liujiaxin-Desktop-Susceptibility-"
       "qmc-platform-master/851bd175-76db-4789-ac47-3b3fb836a9b4/scratchpad")
MAP = {"sowave": "son", "swave": "sext", "dwave": "d", "dd12wave": "dxy"}
CH = [("son", "on-site $s$", "#8c8c8c"), ("sext", "extended $s$", "#1b7837"),
      ("d", r"$d_{x^2-y^2}$", "#2166ac"), ("dxy", r"$d_{xy}$", "#b2182b")]
KEYS = [c for c, _, _ in CH]
RGB = np.array([to_rgb(c) for _, _, c in CH])
BAND = 0.22

f = pd.read_csv(f"{SCR}/fortran_L10.csv"); f["ch"] = f.channel.map(MAP)
EQ = f.pivot_table(index=["n", "delta"], columns="ch", values="vertex").reset_index()

p = pd.read_csv(f"{D}/pairing_master.csv"); p = p[(p.L == 10) & (p.U == 4.0)]
UN = (p.groupby(["n", "delta"])[["chi_son", "chi_sext", "chi_d", "chi_dxy"]]
      .mean().rename(columns=lambda c: c.replace("chi_", "")).reset_index())
UN = UN[UN.delta <= 0.4]                    # match the Fortran delta range

nf = np.linspace(0.5, 1.0, 320); df = np.linspace(0.0, 0.4, 320)
NG, DG = np.meshgrid(nf, df)

fig, ax = plt.subplots(1, 2, figsize=(14.5, 5.8), sharey=True)
for k, (T, lab) in enumerate([(EQ, "(a) equal time   " r"$\tau=0$"),
                              (UN, "(b) unequal time   " r"$\int d\tau$")]):
    a = ax[k]
    pts = T[["n", "delta"]].values
    F = np.stack([gaussian_filter(
        griddata(pts, T[c].values, (NG, DG), method="linear",
                 fill_value=np.nan), 4.0) for c in KEYS])
    order = np.argsort(-F, axis=0)
    i1, i2 = order[0], order[1]
    yy, xx = np.indices(i1.shape)
    spread = np.maximum(F.max(0) - F.min(0), 1e-12)
    frac = np.clip((F[i1, yy, xx] - F[i2, yy, xx]) / (BAND * spread), 0, 1)
    img = (0.5 + 0.5 * frac)[..., None] * RGB[i1] + (0.5 * (1 - frac))[..., None] * RGB[i2]
    img = np.where(np.isnan(F[0])[..., None], 1.0, img)      # blank outside the hull
    a.imshow(img, origin="lower", aspect="auto",
             extent=[nf.min(), nf.max(), df.min(), df.max()],
             interpolation="bilinear")
    a.plot(T.n, T.delta, "o", ms=3.5, color="k", alpha=0.55, lw=0)
    a.set_xlabel(r"Filling  $n$"); a.set_title(lab, fontsize=18)
    a.set_xticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
ax[0].set_ylabel(r"Anisotropy  $\delta$")
fig.legend(handles=[Patch(fc=c, ec="k", lw=0.6, label=l) for _, l, c in CH],
           fontsize=14, loc="center left", bbox_to_anchor=(0.905, 0.5),
           frameon=False, handlelength=1.8, handleheight=1.2, labelspacing=0.8)
fig.tight_layout(rect=[0, 0, 0.90, 1])
fig.savefig(f"{D}/../figs/fig_eq_vs_uneq_phase.png", dpi=600,
            bbox_inches="tight", facecolor="white")

for T, lab in [(EQ, "EQUAL TIME"), (UN, "UNEQUAL TIME")]:
    w = T[KEYS].values.argmax(axis=1)
    print(f"{lab}: " + ", ".join(
        f"{CH[i][1].replace('$','')} {int((w == i).sum())}/{len(T)}" for i in range(4)))
