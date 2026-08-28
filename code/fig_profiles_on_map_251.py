# --- Delta_tot heat map per lattice size, with the delta-cut profiles drawn on it. ---
# Paste-and-run on 251. One panel per L.
#
# Background: Delta_tot/N over (U, delta), the same heat map as before.
# Overlaid: the SAME curves as the line-cut figure. At each anchor delta_0 the
# profile Delta_tot(U) is drawn as a horizontal excursion from a dotted baseline at
# delta_0, so the curve sticks out to the right in proportion to Delta_tot. That
# puts the curve and the colour it came from in the same frame: where the profile
# bulges, the map should be red.
#
# SCALE sets how far a full-strength Delta_tot pushes the curve, in delta units.
# With SCALE=1.0 the largest Delta_tot (~0.09) displaces the curve by 0.09 in delta,
# just under the 0.1 spacing between anchors, so the three profiles never collide.
import os, glob, numpy as np, pandas as pd, matplotlib.pyplot as plt, matplotlib
import matplotlib.patheffects as pe
from scipy.interpolate import RegularGridInterpolator

# Baseline colour: magenta is the one strong hue jet never produces, so it stays
# visible from the deep blue through to the dark red. White vanished on the yellow.
BASE    = "#FF00FF"
ANCHORS = [0.1, 0.2, 0.3]
SIZES   = [8, 10, 12]
SCALE   = 1.0

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
    "font.size":14,"axes.labelsize":19,"legend.fontsize":13,
    "xtick.labelsize":15,"ytick.labelsize":15,"axes.linewidth":1.4})

f = pd.read_csv(find("fortran_L8_L10_L12_dedup.csv"))
f = f[np.isclose(f.n, 1.0) & (f.U >= 2)]

gx, gy = np.meshgrid(np.linspace(.1, .7, 400), np.linspace(2., 5., 400))
LEVBG  = np.linspace(0, f.dtot_N.max(), 60)
HALO   = [pe.withStroke(linewidth=3.0, foreground="w")]

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
    for d0 in ANCHORS:
        g = s[np.isclose(s.delta, d0)].sort_values("U")
        bl = a.axvline(d0, color=BASE, lw=2.0, ls=(0, (2.5, 2.0)), zorder=3)
        bl.set_path_effects([pe.withStroke(linewidth=3.6, foreground="k")])
        x = d0 + SCALE * g.dtot_N.values
        ln, = a.plot(x, g.U.values, marker="o", ms=6, lw=2.4, color="k", zorder=5)
        ln.set_path_effects(HALO)
        for xi, yi in zip(x, g.U.values):          # tie each point to its baseline
            a.plot([d0, xi], [yi, yi], color=BASE, lw=1.1, alpha=.75, zorder=4)
    a.set_title(rf"$L={L}$", fontsize=20, pad=10)
    a.set_xlabel(r"anisotropy  $\delta$", labelpad=7)
    a.set_xlim(.1, .7); a.set_ylim(2., 5.)
    a.set_xticks([.1,.2,.3,.4,.5,.6,.7]); a.set_yticks([2,3,4,5])
    a.tick_params(direction="out", top=False, right=False, length=6, width=1.4)
ax[0].set_ylabel(r"$U/t$", labelpad=8)

H = [plt.Line2D([],[],color=BASE,ls=(0,(2.5,2.0)),lw=2.0),
     plt.Line2D([],[],color="k",ls="-",lw=2.4,marker="o",ms=6)]
ax[0].legend(H, [r"baseline at $\delta_0$", r"$\Delta_{\rm tot}(U)$ profile"],
             loc="upper right", frameon=True, framealpha=.9, edgecolor="0.75",
             handlelength=1.8, labelspacing=.3, borderpad=.45)
cb = fig.colorbar(im, ax=ax.tolist(), pad=.012, fraction=.017)
cb.set_label(r"$\Delta_{\rm tot}/N$", fontsize=20, labelpad=9)
cb.ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%.3f"))
cb.ax.tick_params(labelsize=14)
fig.suptitle(r"half filling. Each profile is $\Delta_{\rm tot}(U)$ at its own $\delta_0$,"
             r" drawn rightward from the dotted baseline.", fontsize=16, y=1.00)
plt.show()
