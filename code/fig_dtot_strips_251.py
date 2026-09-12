# --- Narrow strip around each delta: the three lattice sizes side by side. ---
# Paste-and-run on 251.
#
# One panel per delta. The x axis is zoomed to a small window around it, and the
# three lattice sizes are placed at slight horizontal offsets so they do not
# overlap. Each marker is FILLED with its own Delta_tot/N on the same jet scale
# as the background, so the comparison is by colour:
#
#     background = the L=12 map, interpolated
#     markers    = L = 8 (circle), 10 (triangle), 12 (square)
#
# A marker that blends into the background agrees with L=12 at that (delta, U).
# A marker that stands out against it does not. Marker SHAPE carries the lattice
# size because fill colour is already carrying Delta_tot.
import os, glob, numpy as np, pandas as pd, matplotlib.pyplot as plt, matplotlib
import matplotlib.patheffects as pe
from scipy.interpolate import RegularGridInterpolator

DELTAS = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
NCOL   = 3
LBG    = 12                       # lattice supplying the background
VMAX   = 0.0916
HALF   = 0.058                    # half width of each strip, in delta
OFF    = {8: -0.030, 10: 0.0, 12: +0.030}
MK     = {8: "o", 10: "^", 12: "s"}

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
    "font.size":14,"axes.labelsize":18,"legend.fontsize":13,
    "xtick.labelsize":13,"ytick.labelsize":14,"axes.linewidth":1.4})

f = pd.read_csv(find("fortran_L8_L10_L12_dedup.csv"))
f = f[np.isclose(f.n, 1.0) & (f.U >= 2)]
bg = f[f.L == LBG]

def strip(d0):
    """Background for one strip: the L=LBG value AT delta0, held constant across
    the strip and interpolated only in U.

    It must not vary horizontally. If it did, a marker at the L=8 offset would be
    sitting on the L=12 value at a DIFFERENT delta, and 'blends in' would stop
    meaning 'agrees'. Constant in x makes the comparison exact at every offset."""
    gx, gy = np.meshgrid(np.linspace(d0-HALF, d0+HALF, 60),
                         np.linspace(2., 5., 300))
    col = bg[np.isclose(bg.delta, d0)].sort_values("U")
    prof = np.interp(gy[:, 0], col.U.values, col.dtot_N.values)
    return gx, gy, np.repeat(prof[:, None], gx.shape[1], axis=1)

NORM = matplotlib.colors.Normalize(0, VMAX)
LEV  = np.linspace(0, VMAX, 60)
NROW = int(np.ceil(len(DELTAS)/NCOL))
fig, axg = plt.subplots(NROW, NCOL, figsize=(4.6*NCOL, 4.6*NROW),
                        sharey=True, gridspec_kw={"wspace":.13, "hspace":.24})
ax = axg.ravel()
for a in ax[len(DELTAS):]: a.set_visible(False)

for i, (a, d0) in enumerate(zip(ax, DELTAS)):
    gx, gy, Z = strip(d0)
    im = a.contourf(gx, gy, Z, levels=LEV, cmap="jet", extend="both")
    for L in (8, 10, 12):
        g = f[(f.L == L) & np.isclose(f.delta, d0)].sort_values("U")
        if not len(g): continue
        a.scatter([d0 + OFF[L]]*len(g), g.U, c=g.dtot_N, cmap="jet", norm=NORM,
                  marker=MK[L], s=190, edgecolors="k", linewidths=1.5, zorder=5)
    a.set_title(rf"$\delta = {d0:g}$", fontsize=19, pad=9)
    a.set_xlim(d0-HALF, d0+HALF); a.set_ylim(2., 5.)
    a.set_xticks([d0+OFF[L] for L in (8, 10, 12)])
    a.set_xticklabels([f"$L={L}$" for L in (8, 10, 12)])
    a.set_yticks([2, 3, 4, 5])
    a.tick_params(direction="out", top=False, right=False, length=6, width=1.4)
    if i % NCOL == 0: a.set_ylabel(r"$U/t$", labelpad=8)

H = [plt.Line2D([], [], ls="none", marker=MK[L], ms=11, mfc="0.85", mec="k",
                mew=1.5, label=f"$L={L}$") for L in (8, 10, 12)]
fig.legend(handles=H, loc="lower center", ncol=3, frameon=False,
           bbox_to_anchor=(.5, -.02), fontsize=15, handletextpad=.5,
           columnspacing=2.4)
cb = fig.colorbar(im, ax=list(ax), pad=.014, fraction=.019)
cb.set_label(r"$\Delta_{\rm tot}/N$", fontsize=19, labelpad=9)
cb.ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%.3f"))
cb.ax.tick_params(labelsize=13)
fig.suptitle(rf"half filling.  Background is the $L={LBG}$ map; a marker that blends"
             rf" into it agrees with $L={LBG}$.", fontsize=16, y=0.975)
plt.show()

print(f"Delta_tot/N, spread across lattice size at each (delta, U)")
print(f"{'delta':>7}{'U':>6}{'L=8':>10}{'L=10':>10}{'L=12':>10}{'max/min':>10}")
for d0 in DELTAS:
    for U in np.sort(f.U.unique()):
        v = [float(f[(f.L == L) & np.isclose(f.delta, d0) & (f.U == U)].dtot_N.iloc[0])
             for L in (8, 10, 12)]
        sp = max(v)/min(v) if min(v) > 0 else np.nan
        flag = "  <-- sizes disagree" if sp > 2 else ""
        print(f"{d0:>7.1f}{U:>6.1f}" + "".join(f"{x:>10.4f}" for x in v)
              + f"{sp:>10.1f}" + flag)
