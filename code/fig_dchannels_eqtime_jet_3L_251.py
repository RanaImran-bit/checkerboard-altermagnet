# --- equal time only, both d channels over Delta_tot, three sizes, U up to 5. ---
# Jet background. Contours black, told apart by line style, because jet spans blue
# through red and any coloured contour disappears somewhere on it. The crossover is
# a white line with a black casing so it reads over every jet value, and contour
# labels carry a white halo for the same reason.
# The measured grid is delta in steps of 0.1 and U in {2, 3, 3.5, 4, 4.5, 5}.
# Contour position between those points is interpolation, so say so in the caption.
import os, glob, numpy as np, pandas as pd, matplotlib.pyplot as plt, matplotlib
import matplotlib.patheffects as pe
from scipy.interpolate import RegularGridInterpolator

BGMAP = "jet"
STYLE = "black"                     # "black" or "colour"
CD, CX = ("k", "k") if STYLE == "black" else ("#0072B2", "#D55E00")
LSEL   = [8, 10, 12]
UMAX   = 5.0
INTERP = "linear"                   # see the note on grid() before changing this

def find(name):
    """Locate a data file without caring which directory you launched from."""
    for c in [name, f"../data/{name}", f"data/{name}",
              os.path.expanduser(f"~/results/figdata/{name}"),
              os.path.expanduser(f"~/results/{name}"),
              os.path.expanduser(f"~/checkerboard-altermagnet/data/{name}")]:
        if os.path.exists(c): return c
    hits = glob.glob(os.path.expanduser(f"~/results/**/{name}"), recursive=True)
    if hits: return hits[0]
    raise SystemExit(f"cannot find {name}. Try:  find ~ -name {name}")

plt.rcParams.update({"font.family":"serif","mathtext.fontset":"dejavuserif",
    "font.size":16,"axes.labelsize":18,"legend.fontsize":17,
    "xtick.labelsize":16,"ytick.labelsize":16,"axes.linewidth":1.4})

f = pd.read_csv(find("fortran_L8_L10_L12_dedup.csv"))
f = f[np.isclose(f.n,1.0) & (f.U>=2) & (f.U<=UMAX) & f.L.isin(LSEL)]
LS = sorted(f.L.unique())

gx, gy = np.meshgrid(np.linspace(.1,.7,400), np.linspace(2.,UMAX,400))
LEVBG  = np.linspace(0, f.dtot_N.max(), 60)     # one Delta_tot scale for all panels
LD, LX = [.08,.16,.24], [.12,.24,.36]
HALO   = [pe.withStroke(linewidth=2.6, foreground="w")]

def halo(cs):
    """White outline on a contour set, across matplotlib versions.

    Up to 3.7 a ContourSet held .collections; from 3.8 it is itself a Collection
    and 3.10 removed .collections outright."""
    if hasattr(cs, "set_path_effects"):     # 3.8+, touching .collections warns
        cs.set_path_effects(HALO)
    else:
        for c in cs.collections: c.set_path_effects(HALO)

def grid(s, col):
    """Interpolate a column onto the display mesh.

    Do NOT use scipy.griddata here. It treats the points as scattered and
    triangulates them, but this data is a rectangular lattice, where every cell
    can be split along either diagonal and both are valid Delaunay triangulations.
    Qhull's choice is not stable: adding the U=5 row flipped 22 of the 48
    triangles lying below U=4.5, moving drawn contours by up to 20% of the panel
    range at U=3, far from any new data. The figure was then not a function of the
    data alone.

    RegularGridInterpolator works on the rectangle directly. With method="linear"
    (bilinear) each cell depends only on its own four corners, so the result is
    unique, local, reproducible, and bounded by the measurements. Verified: the
    shared window is identical to machine precision with and without the U=5 row,
    and the range matches the measured range to 1e-4.

    method="cubic" is smoother but is neither: it moves the U=3 contours by 7.9%
    when the U=5 row is added, and undershoots to -0.05 where the smallest
    measurement is +0.02. Use it for slides, never for the paper."""
    p = s.pivot_table(index="U", columns="delta", values=col)
    assert not p.isna().any().any(), "grid has holes, pivot_table left NaN"
    fn = RegularGridInterpolator((p.index.values, p.columns.values), p.values,
                                 method=INTERP, bounds_error=False, fill_value=None)
    return fn(np.stack([gy.ravel(), gx.ravel()], -1)).reshape(gx.shape)

def place(cs, ax, taken):               # labels off the frame and apart from each other
    x0,x1 = ax.get_xlim(); y0,y1 = ax.get_ylim(); out = []
    for segs in cs.allsegs:
        best, score = None, -9.
        for sg in segs:
            if len(sg) < 2: continue
            u = (sg[:,0]-x0)/(x1-x0); w = (sg[:,1]-y0)/(y1-y0)
            sc = np.minimum.reduce([u,1-u,w,1-w])
            for tu,tw in taken: sc = np.minimum(sc, .55*np.hypot(u-tu, w-tw))
            k = int(np.argmax(sc))
            if sc[k] > score: score, best = sc[k], (tuple(sg[k]), (u[k],w[k]))
        if best is not None and score > .045: out.append(best[0]); taken.append(best[1])
    return out

fig, ax = plt.subplots(1, len(LS), figsize=(6.2*len(LS),5.6), sharex=True, sharey=True,
                       squeeze=False, gridspec_kw={"wspace":.10})
ax = ax[0]
for a, L in zip(ax, LS):
    s = f[f.L==L]
    im = a.contourf(gx, gy, grid(s, "dtot_N"), levels=LEVBG, cmap=BGMAP, extend="both")
    cd = a.contour(gx, gy, grid(s, "d"),   levels=LD, colors=CD, linewidths=2.2)
    cx = a.contour(gx, gy, grid(s, "dxy"), levels=LX, colors=CX, linewidths=2.2,
                   linestyles="dashed")
    halo(cd); halo(cx)
    taken = []
    t1 = a.clabel(cd, manual=place(cd,a,taken), inline=True, inline_spacing=12, fontsize=15, fmt="%.2f")
    t2 = a.clabel(cx, manual=place(cx,a,taken), inline=True, inline_spacing=12, fontsize=15, fmt="%.2f")
    for t in list(t1) + list(t2): t.set_path_effects(HALO)
    t = s.assign(diff=s.dxy - s.d)
    ceq = grid(t, "diff")                             # crossover: white, black cased
    a.contour(gx, gy, ceq, levels=[0.], colors="k", linewidths=6.0)
    a.contour(gx, gy, ceq, levels=[0.], colors="w", linewidths=3.0)
    a.set_title(rf"equal time ($L={L}$, half filling)", fontsize=17, pad=12)
    a.set_xlabel(r"anisotropy  $\delta$", labelpad=8)
    a.set_xticks([.2,.3,.4,.5,.6]); a.set_ylim(2., UMAX); a.set_xlim(.1,.7)
    a.tick_params(direction="out", top=False, right=False, length=6, width=1.4)

ax[0].set_ylabel(r"$U/t$", labelpad=8)
H = [plt.Line2D([],[],color=CD,ls="-",lw=2.2), plt.Line2D([],[],color=CX,ls="--",lw=2.2),
     plt.Line2D([],[],color="w",ls="-",lw=3.0,path_effects=[pe.withStroke(linewidth=6.0,foreground="k")])]
fig.legend(H, [r"$d_{x^2-y^2}$", r"$d_{xy}$", r"$d_{xy}=d_{x^2-y^2}$"], loc="lower center",
           ncol=3, frameon=False, bbox_to_anchor=(.5,-.19), handlelength=2.4, columnspacing=2.6)
cb = fig.colorbar(im, ax=ax.tolist(), pad=.014, fraction=.020)
cb.set_label(r"$\Delta_{\rm tot}/N$", fontsize=20, labelpad=10)
cb.ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%.3f"))
cb.ax.tick_params(labelsize=15)
plt.show()
