"""Pairing channel competition: EQUAL-TIME vs UNEQUAL-TIME.

Two independent data sets, both L = 10, U = 4, filling 0.5 to 1:

  EQUAL-TIME   the group's Fortran CPQMC, Vertex_*.dat at k = (0,0).
               tau = 0 only, no imaginary-time integration. 27 runs.
  UNEQUAL-TIME our Python CPQMC, chi_a = int C_a(tau) dtau. pairing_master.csv.

Channel correspondence, read off mc2duph.f90 (lines 1353-1355, 1386):
    sowave   -> on-site s        swave  -> extended s
    dwave    -> d_x2-y2          dd12wave -> d_xy

CAUTION on absolute magnitudes. The two codes use different conventions: the
Fortran vertex is full - L*Unpair (implied factor 10.03 +/- 0.15 at L = 10),
ours is full - bubble. So the two sectors are NOT on a common scale and only
the ORDERING of channels, and its evolution with n and delta, is compared here.
Each panel is normalised to its own largest channel for that reason.
"""
import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.patches import Patch

mpl.rcParams.update({
    "font.family": "serif", "font.serif": ["DejaVu Serif"],
    "mathtext.fontset": "dejavuserif", "axes.linewidth": 1.1,
    "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True,
    "ytick.right": True, "legend.frameon": False,
    "axes.labelsize": 14, "xtick.labelsize": 11, "ytick.labelsize": 11,
    "figure.dpi": 130})

D = "/Users/liujiaxin/Desktop/checkerboard-altermagnet/data"
SCR = ("/private/tmp/claude-501/-Users-liujiaxin-Desktop-Susceptibility-"
       "qmc-platform-master/851bd175-76db-4789-ac47-3b3fb836a9b4/scratchpad")
MAP = {"sowave": "son", "swave": "sext", "dwave": "d", "dd12wave": "dxy"}
CH = [("son", "on-site $s$", "#8c8c8c", "D"),
      ("sext", "extended $s$", "#1b7837", "s"),
      ("d", r"$d_{x^2-y^2}$", "#2166ac", "o"),
      ("dxy", r"$d_{xy}$", "#b2182b", "^")]
DL = [0.0, 0.1, 0.2, 0.3, 0.4]
NSITES = 100          # L = 10

# ---- equal time (Fortran) ----
f = pd.read_csv(f"{SCR}/fortran_L10.csv")
f["ch"] = f.channel.map(MAP)
EQ = f.pivot_table(index=["n", "delta"], columns="ch", values="vertex")

# ---- unequal time (ours) ----
p = pd.read_csv(f"{D}/pairing_master.csv")
p = p[(p.L == 10) & (p.U == 4.0)]
UN = (p.groupby(["n", "delta"])[["chi_son", "chi_sext", "chi_d", "chi_dxy"]]
      .mean().rename(columns=lambda c: c.replace("chi_", "")))

fig, ax = plt.subplots(2, len(DL), figsize=(4.0 * len(DL), 8.0), sharex=True)
for j, dl in enumerate(DL):
    for row, (T, lab) in enumerate([(EQ, "equal time"), (UN, "unequal time")]):
        a = ax[row, j]
        s = T.xs(dl, level="delta").sort_index()
        # Our panel keeps its physical units, chi/N. The Fortran panel is
        # normalised to its own maximum because the two codes use different
        # vertex conventions (theirs is full - L*Unpair, ours full - bubble),
        # so their magnitudes are not on a common scale. Only the ORDERING of
        # channels may be compared between the rows.
        scale = np.abs(s.values).max() if row == 0 else NSITES
        for c, lb, col, mk in CH:
            if c not in s: continue
            a.plot(s.index, s[c] / scale, marker=mk, ms=5, lw=1.4,
                   color=col, label=lb)
        a.axhline(0, color="k", lw=1.0, ls="--")
        a.set_xlim(0.48, 1.02)
        if row == 0: a.set_title(rf"$\delta = {dl:g}$", fontsize=14)
        if row == 1: a.set_xlabel(r"filling  $n$")
        if j == 0:
            a.set_ylabel("equal time (Fortran)\nvertex / max"  if row == 0 else
                         r"unequal time (ours)" "\n" r"vertex  $\chi_\alpha/N$")
ax[0, 0].legend(fontsize=9, loc="best")
fig.tight_layout()
fig.savefig(f"{D}/../figs/fig_eq_vs_uneq_channels.png", dpi=600,
            bbox_inches="tight", facecolor="white")

# ---- which channel leads, printed ----
print("LEADING CHANNEL (largest vertex), L=10, U=4\n")
for T, lab in [(EQ, "EQUAL TIME (Fortran)"), (UN, "UNEQUAL TIME (ours)")]:
    print(f"=== {lab} ===")
    ns = sorted(T.index.get_level_values("n").unique())
    print(f"{'delta':>6} | " + " ".join(f"{x:>7.3f}" for x in ns))
    for dl in DL:
        row = f"{dl:6.1f} | "
        for nv in ns:
            try:
                r = T.loc[(nv, dl)]
                row += f"{r.idxmax():>7} "
            except KeyError:
                row += f"{'-':>7} "
        print(row)
    print()
