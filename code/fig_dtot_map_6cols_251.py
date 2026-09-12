# --- Delta_tot(delta, U) heat map, six panels, one delta column marked in each. ---
# Paste-and-run on 251.
#
# Same background in every panel: the polarisation Delta_tot/N over the (delta, U)
# plane at L=12, the largest lattice with a complete 7x7 grid. What changes panel to
# panel is which delta column is picked out, with a marker at every measured U and
# the value written beside it. So the panels are six readings of one map.
#
# VMAX is the 0.0916 the other heat maps use, so colours are directly comparable
# with them. LBG switches the background lattice size (8, 10 or 12).
import os, glob, numpy as np, pandas as pd, matplotlib.pyplot as plt, matplotlib
import matplotlib.patheffects as pe
from scipy.interpolate import RegularGridInterpolator

DELTAS = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
NCOL   = 3
LBG    = 12                 # lattice size supplying the background map
VMAX   = 0.0916
CUTC   = "#FF00FF"          # magenta: the one strong hue jet never produces

def find(name):
    for c in [name, f"../data/{name}", f"data/{name}",
              os.path.expanduser(f"~/results/figdata/{name}"),
              os.path.expanduser(f"~/results/{name}"),
              os.path.expanduser(f"~/checkerboard-altermagnet/data/{name}")]:
        if os.path.exists(c): return c
    h = glob.glob(os.path.expanduser(f"~/results/**/{name}"), recursive=True)
    if h: return h[0]
    raise SystemExit(f"cannot find {name}. Try:  find ~ -name {name}")

plt.rcParams.update({"font.family":"serif","mathtext.fontset":"dejavuserif",
    "font.size":14,"axes.labelsize":18,"xtick.labelsize":14,"ytick.labelsize":14,
    "axes.linewidth":1.4})

f = pd.read_csv(find("fortran_L8_L10_L12_dedup.csv"))
f = f[np.isclose(f.n, 1.0) & (f.U >= 2) & (f.L == LBG)]

gx, gy = np.meshgrid(np.linspace(.1, .7, 400), np.linspace(2., 5., 400))

def grid(col):
    """Bilinear on the measured rectangle. NOT scipy.griddata: that triangulates,
    and on a rectangular lattice the triangulation is neither unique nor stable."""
    p = f.pivot_table(index="U", columns="delta", values=col)
    assert not p.isna().any().any(), "grid has holes"
    fn = RegularGridInterpolator((p.index.values, p.columns.values), p.values,
                                 method="linear", bounds_error=False, fill_value=None)
    return fn(np.stack([gy.ravel(), gx.ravel()], -1)).reshape(gx.shape)

BG   = grid("dtot_N")
LEV  = np.linspace(0, VMAX, 60)
CASE = [pe.withStroke(linewidth=3.2, foreground="k")]
Us   = np.sort(f.U.unique())

NROW = int(np.ceil(len(DELTAS)/NCOL))
fig, axg = plt.subplots(NROW, NCOL, figsize=(4.9*NCOL, 4.6*NROW),
                        sharex=True, sharey=True,
                        gridspec_kw={"wspace":.07, "hspace":.20})
ax = axg.ravel()
for a in ax[len(DELTAS):]: a.set_visible(False)

for i, (a, d0) in enumerate(zip(ax, DELTAS)):
    im = a.contourf(gx, gy, BG, levels=LEV, cmap="jet", extend="both")
    ln = a.axvline(d0, color=CUTC, lw=2.4, zorder=4); ln.set_path_effects(CASE)
    col = f[np.isclose(f.delta, d0)].sort_values("U")
    a.plot([d0]*len(col), col.U, "o", ms=8, mfc=CUTC, mec="k", mew=1.3, zorder=5)
    for _, r in col.iterrows():           # write the value beside each point
        side = 12 if d0 < 0.45 else -12
        a.annotate(f"{r.dtot_N:.3f}", xy=(d0, r.U), xytext=(side, 0),
                   textcoords="offset points", fontsize=11.5, color="w",
                   ha="left" if side > 0 else "right", va="center", zorder=6,
                   path_effects=[pe.withStroke(linewidth=2.6, foreground="k")])
    a.set_title(rf"$\delta = {d0:g}$", fontsize=19, pad=9)
    a.set_xlim(.1, .7); a.set_ylim(2., 5.)
    a.set_xticks([.1,.2,.3,.4,.5,.6,.7]); a.set_yticks([2,3,4,5])
    a.tick_params(direction="out", top=False, right=False, length=6, width=1.4)
    if i % NCOL == 0: a.set_ylabel(r"$U/t$", labelpad=8)
    if i >= len(DELTAS) - NCOL: a.set_xlabel(r"anisotropy  $\delta$", labelpad=7)

cb = fig.colorbar(im, ax=list(ax), pad=.013, fraction=.019)
cb.set_label(r"$\Delta_{\rm tot}/N$", fontsize=19, labelpad=9)
cb.ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%.3f"))
cb.ax.tick_params(labelsize=13)
fig.suptitle(rf"half filling, $L={LBG}$.  One map, six columns read off it.",
             fontsize=17, y=0.975)
plt.show()

print(f"Delta_tot/N at L={LBG}, by column")
print(f"{'U':>5}" + "".join(f"{d:>9.1f}" for d in DELTAS))
for U in Us:
    row = [float(f[(f.U == U) & np.isclose(f.delta, d)].dtot_N.iloc[0]) for d in DELTAS]
    star = max(range(len(row)), key=lambda j: row[j])
    print(f"{U:>5.1f}" + "".join(f"{v:>9.4f}" + ("*" if j == star else " ")
                                 for j, v in enumerate(row)))
print("  * marks the largest value in each row: the ridge, moving right as U rises.")
