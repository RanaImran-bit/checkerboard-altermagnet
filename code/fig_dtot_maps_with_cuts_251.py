# --- Delta_tot heat map per lattice size, with the delta cuts marked. ---
# Paste-and-run on 251.
#
# Companion to the line-cut figure. Each vertical line marks one column of the map,
# and that column is exactly what the corresponding curve in the line-cut figure
# plots. The dots sit on the measured U values, so they are the same points.
#
# One meaning per axis here: x is delta, y is U, colour is Delta_tot/N. Nothing is
# displaced or rescaled, so the figure can be read directly.
#
# Read it across panels: does the bright spot at low U survive as L grows?
import os, glob, numpy as np, pandas as pd, matplotlib.pyplot as plt, matplotlib
import matplotlib.patheffects as pe
from scipy.interpolate import RegularGridInterpolator

CUTS  = [0.1, 0.2, 0.3]
SIZES = [8, 10, 12]
CUTC  = "#FF00FF"       # magenta: the one strong hue jet never produces

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
    "font.size":14,"axes.labelsize":19,"legend.fontsize":14,
    "xtick.labelsize":15,"ytick.labelsize":15,"axes.linewidth":1.4})

f = pd.read_csv(find("fortran_L8_L10_L12_dedup.csv"))
f = f[np.isclose(f.n, 1.0) & (f.U >= 2)]

gx, gy = np.meshgrid(np.linspace(.1, .7, 400), np.linspace(2., 5., 400))
LEVBG  = np.linspace(0, f.dtot_N.max(), 60)     # one scale across all panels
CASE   = [pe.withStroke(linewidth=4.0, foreground="k")]

def grid(s, col):
    """Bilinear on the measured rectangle. NOT scipy.griddata: that triangulates,
    and on a rectangular lattice the triangulation is neither unique nor stable."""
    p = s.pivot_table(index="U", columns="delta", values=col)
    assert not p.isna().any().any(), "grid has holes"
    fn = RegularGridInterpolator((p.index.values, p.columns.values), p.values,
                                 method="linear", bounds_error=False, fill_value=None)
    return fn(np.stack([gy.ravel(), gx.ravel()], -1)).reshape(gx.shape)

fig, ax = plt.subplots(1, len(SIZES), figsize=(5.4*len(SIZES), 5.4), sharex=True,
                       sharey=True, gridspec_kw={"wspace":.07})
for a, L in zip(ax, SIZES):
    s = f[f.L == L]
    im = a.contourf(gx, gy, grid(s, "dtot_N"), levels=LEVBG, cmap="jet", extend="both")
    Us = np.sort(s.U.unique())
    for d0 in CUTS:
        ln = a.axvline(d0, color=CUTC, lw=2.2, zorder=4)
        ln.set_path_effects(CASE)
        mk, = a.plot([d0]*len(Us), Us, "o", ms=6, mfc=CUTC, mec="k", mew=1.1, zorder=5)
        a.annotate(f"{d0:g}", xy=(d0, 5.0), xytext=(0, 6), textcoords="offset points",
                   ha="center", va="bottom", fontsize=14, color=CUTC,
                   path_effects=[pe.withStroke(linewidth=2.6, foreground="k")])
    a.set_title(rf"$L={L}$", fontsize=20, pad=22)
    a.set_xlabel(r"anisotropy  $\delta$", labelpad=7)
    a.set_xlim(.1, .7); a.set_ylim(2., 5.)
    a.set_xticks([.1,.2,.3,.4,.5,.6,.7]); a.set_yticks([2,3,4,5])
    a.tick_params(direction="out", top=False, right=False, length=6, width=1.4)
ax[0].set_ylabel(r"$U/t$", labelpad=8)

H = [plt.Line2D([],[],color=CUTC,lw=2.2,marker="o",ms=6,mec="k",mew=1.1)]
ax[0].legend(H, [r"cuts shown in the line figure"], loc="upper right",
             frameon=True, framealpha=.9, edgecolor="0.75",
             handlelength=2.0, borderpad=.45)
cb = fig.colorbar(im, ax=ax.tolist(), pad=.012, fraction=.017)
cb.set_label(r"$\Delta_{\rm tot}/N$", fontsize=20, labelpad=9)
cb.ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%.3f"))
cb.ax.tick_params(labelsize=14)
fig.suptitle("half filling.  Same colour scale in all three panels.", fontsize=16, y=1.03)
plt.show()
