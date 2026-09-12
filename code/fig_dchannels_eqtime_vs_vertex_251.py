# --- equal-time vs time-integrated vertex, both d channels, over Delta_tot. L=12. ---
# Paste-and-run on 251. Needs, anywhere under ~/results or the cwd:
#     fortran_L8_L10_L12_dedup.csv   equal-time pair correlations (Fortran)
#     chi_L12_vertex_summary.csv     vertex, beta=32, N_w=500, tau window 0.8, 6 seeds
import os, glob, numpy as np, pandas as pd, matplotlib.pyplot as plt, matplotlib
from matplotlib.colors import LinearSegmentedColormap
from scipy.interpolate import RegularGridInterpolator

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
    "font.size":16,"axes.labelsize":22,"legend.fontsize":19,
    "xtick.labelsize":18,"ytick.labelsize":18,"axes.linewidth":1.4})
CD, CX = "#0072B2", "#D55E00"
GREY = LinearSegmentedColormap.from_list("g", plt.get_cmap("Greys")(np.linspace(0,.72,256)))

f = pd.read_csv(find("fortran_L8_L10_L12_dedup.csv"))
f = f[np.isclose(f.n,1.0) & (f.L==12) & (f.U>=2)]
v = pd.read_csv(find("chi_L12_vertex_summary.csv")); v = v[v.U>=2]
UMAX = min(f.U.max(), v.U.max())        # equal time reaches U=5, the vertex 4.5
f = f[f.U <= UMAX]; bg = f[["U","delta","dtot_N"]]

# per panel: (data, d levels, dxy levels, title). Levels differ, the two
# quantities are ~20x apart in magnitude, so each panel carries its own numbers.
PAN = [(f, [.08,.16,.24], [.12,.24,.36], "equal time"),
       (v, [1.,2.,3.],    [1.,3.,5.],    "vertex, time integrated")]

gx, gy = np.meshgrid(np.linspace(.1,.7,400), np.linspace(2.,UMAX,400))
LEVBG  = np.linspace(0, bg.dtot_N.max(), 60)

def grid(s, col):
    """Interpolate a column onto the display mesh.

    NOT scipy.griddata. That treats the points as scattered and triangulates, but
    this data is a rectangular lattice where every cell can split along either
    diagonal, both valid Delaunay. Qhull's tie-break is unstable: dropping one row
    of the scan moved contours elsewhere in the panel by up to 64% of its range.
    The figure was then not a function of the data alone.

    RegularGridInterpolator with method="linear" is bilinear per cell: unique,
    local, reproducible, bounded by the measurements."""
    p = s.pivot_table(index="U", columns="delta", values=col)
    assert not p.isna().any().any(), "grid has holes, pivot_table left NaN"
    fn = RegularGridInterpolator((p.index.values, p.columns.values), p.values,
                                 method="linear", bounds_error=False, fill_value=None)
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

BGZ = grid(bg, "dtot_N")
fig, ax = plt.subplots(1, 2, figsize=(12.4,5.6), sharex=True, sharey=True,
                       gridspec_kw={"wspace":.13})
for a, (s, ld, lx, ttl) in zip(ax, PAN):
    im = a.contourf(gx, gy, BGZ, levels=LEVBG, cmap=GREY, extend="both")
    cd = a.contour(gx, gy, grid(s, "d"),   levels=ld, colors=CD, linewidths=2.2)
    cx = a.contour(gx, gy, grid(s, "dxy"), levels=lx, colors=CX, linewidths=2.2,
                   linestyles="dashed")
    taken = []; fmt = "%.2f" if max(ld) < 1 else "%.1f"
    a.clabel(cd, manual=place(cd,a,taken), inline=True, inline_spacing=12, fontsize=15, fmt=fmt)
    a.clabel(cx, manual=place(cx,a,taken), inline=True, inline_spacing=12, fontsize=15, fmt=fmt)
    a.contour(gx, gy, grid(s.assign(_diff=s.dxy-s.d), "_diff"), levels=[0.], colors="k", linewidths=3.4)
    a.plot(s.delta, s.U, ".", ms=3.4, color="0.35", zorder=5)
    a.set_title(ttl, fontsize=25, pad=12)
    a.set_xlabel(r"anisotropy  $\delta$", labelpad=8)
    a.set_xticks([.2,.3,.4,.5,.6]); a.set_ylim(2., UMAX); a.set_xlim(.1,.7)
    a.tick_params(direction="out", top=False, right=False, length=6, width=1.4)

ax[0].set_ylabel(r"$U/t$", labelpad=8)
H = [plt.Line2D([],[],color=CD,ls="-",lw=2.2), plt.Line2D([],[],color=CX,ls="--",lw=2.2),
     plt.Line2D([],[],color="k",ls="-",lw=3.4)]
fig.legend(H, [r"$d_{x^2-y^2}$", r"$d_{xy}$", r"$d_{xy}=d_{x^2-y^2}$"], loc="lower center",
           ncol=3, frameon=False, bbox_to_anchor=(.5,-.21), handlelength=2.4, columnspacing=2.6)
cb = fig.colorbar(im, ax=ax.tolist(), pad=.014, fraction=.020)
cb.set_label(r"$\Delta_{\rm tot}/N$", fontsize=24, labelpad=12)
cb.ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%.3f"))
cb.ax.tick_params(labelsize=17)
fig.suptitle(r"$L=12$, half filling", fontsize=25, y=1.06)
plt.show()
