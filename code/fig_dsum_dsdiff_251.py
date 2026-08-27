# --- heat maps of d_xy + d_x2-y2 and d_xy - d_x2-y2. Paste-and-run on 251. ---
# No contours to squint at: the channel combination IS the colour.
#   top row     SUM,  d_xy + d_x2-y2  = total d-wave pairing strength, sign-blind
#   bottom row  DIFF, d_xy - d_x2-y2  = which channel wins, and by how much
# The sum uses jet since it is positive everywhere. The difference changes sign, so
# it uses a diverging map pinned at zero: red = d_xy wins, blue = d_x2-y2 wins,
# white = the crossover. Left column equal time, right column time integrated.
import os, glob, numpy as np, pandas as pd, matplotlib.pyplot as plt, matplotlib
from scipy.interpolate import RegularGridInterpolator

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
    "font.size":15,"axes.labelsize":19,"xtick.labelsize":15,"ytick.labelsize":15,
    "axes.linewidth":1.4})

f = pd.read_csv(find("fortran_L8_L10_L12_dedup.csv"))
f = f[np.isclose(f.n,1.0) & (f.L==12) & (f.U>=2)]
v = pd.read_csv(find("chi_L12_vertex_summary.csv")); v = v[v.U>=2]
UMAX = min(f.U.max(), v.U.max()); f = f[f.U <= UMAX]
gx, gy = np.meshgrid(np.linspace(.1,.7,400), np.linspace(2.,UMAX,400))

def grid(s, col):
    """Bilinear on the measured rectangle. NOT griddata: that triangulates, and on
    a rectangular lattice the triangulation is neither unique nor stable."""
    p = s.pivot_table(index="U", columns="delta", values=col)
    assert not p.isna().any().any(), "grid has holes"
    fn = RegularGridInterpolator((p.index.values, p.columns.values), p.values,
                                 method="linear", bounds_error=False, fill_value=None)
    return fn(np.stack([gy.ravel(), gx.ravel()], -1)).reshape(gx.shape)

def bar(fig, im, a, lab):
    """Colorbar with few, round ticks. The default put four decimals on every tick
    and the rotated label landed on the next column."""
    cb = fig.colorbar(im, ax=a, pad=.025, fraction=.05)
    cb.locator = matplotlib.ticker.MaxNLocator(6, prune=None); cb.update_ticks()
    cb.ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%.2f"))
    cb.ax.tick_params(labelsize=13)
    cb.set_label(lab, fontsize=17, labelpad=6)
    return cb

COLS = [(f, "equal time"), (v, "time integrated")]
fig, ax = plt.subplots(2, 2, figsize=(12.6, 9.4), sharex=True, sharey=True,
                       gridspec_kw={"wspace":.42, "hspace":.22})
for j, (s, lab) in enumerate(COLS):
    t = s.assign(_sum=s.dxy + s.d, _dif=s.dxy - s.d)

    a = ax[0, j]                                    # SUM
    Z = grid(t, "_sum")
    im = a.contourf(gx, gy, Z, levels=np.linspace(Z.min(), Z.max(), 80),
                    cmap="jet", extend="both")
    bar(fig, im, a, r"$d_{xy}+d_{x^2-y^2}$")
    a.set_title(rf"{lab}:  total $d$-wave strength", fontsize=17, pad=10)

    a = ax[1, j]                                    # DIFFERENCE
    Z = grid(t, "_dif"); m = np.abs(Z).max()
    im = a.contourf(gx, gy, Z, levels=np.linspace(-m, m, 80), cmap="RdBu_r",
                    vmin=-m, vmax=m, extend="both")
    a.contour(gx, gy, Z, levels=[0.], colors="k", linewidths=3.0)
    bar(fig, im, a, r"$d_{xy}-d_{x^2-y^2}$")
    a.set_title(rf"{lab}:  which channel wins", fontsize=17, pad=10)

for a in ax.ravel():
    a.set_xticks([.2,.3,.4,.5,.6]); a.set_xlim(.1,.7); a.set_ylim(2., UMAX)
    a.tick_params(direction="out", top=False, right=False, length=6, width=1.4)
for a in ax[1]: a.set_xlabel(r"anisotropy  $\delta$", labelpad=8)
for a in ax[:,0]: a.set_ylabel(r"$U/t$", labelpad=8)
fig.suptitle(r"$L=12$, half filling.   black line: $d_{xy}=d_{x^2-y^2}$",
             fontsize=18, y=0.965)
plt.show()
