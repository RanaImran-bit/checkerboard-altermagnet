# --- both d channels, EQUAL TIME only, over Delta_tot. Paste-and-run on 251. ---
# LSEL = [8,10,12] for three panels, [12] for the single panel.
# Background Blues truncated at 0.55, contours navy + vermillion: the original palette.
import os, glob, numpy as np, pandas as pd, matplotlib.pyplot as plt, matplotlib
from matplotlib.colors import LinearSegmentedColormap
from scipy.interpolate import RegularGridInterpolator

LSEL          = [8, 10, 12]        # <- [12] for the single panel
BGMAP, BGTOP  = "Blues", 0.55      # ramp, and where to truncate it
CD, CX        = "#08306B", "#D55E00"

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
    "font.size":16,"axes.labelsize":22,"legend.fontsize":19,
    "xtick.labelsize":18,"ytick.labelsize":18,"axes.linewidth":1.4})

d  = pd.read_csv(find("fortran_L8_L10_L12_dedup.csv"))
d  = d[np.isclose(d.n,1.0) & (d.U>=2) & d.L.isin(LSEL)]
LS = sorted(d.L.unique())
CMAP = LinearSegmentedColormap.from_list("bg", plt.get_cmap(BGMAP)(np.linspace(0,BGTOP,256)))
gx, gy = np.meshgrid(np.linspace(d.delta.min(),d.delta.max(),400),
                     np.linspace(d.U.min(),d.U.max(),400))
LEVBG  = np.linspace(0, d.dtot_N.max(), 60)
LD, LX = [.08,.16,.24], [.12,.24,.36]

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

def peak(x, y):
    i = int(np.argmax(y))
    if i in (0, len(y)-1): return x[i]
    den = y[i-1] - 2*y[i] + y[i+1]
    return x[i] if den == 0 else x[i] + .5*(y[i-1]-y[i+1])/den*(x[i+1]-x[i])

def place(cs, ax, taken):
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
                       squeeze=False, gridspec_kw={"wspace":.13})
ax = ax[0]
for a, L in zip(ax, LS):
    s = d[d.L==L]
    im = a.contourf(gx, gy, grid(s, "dtot_N"), levels=LEVBG, cmap=CMAP, extend="both")
    cd = a.contour(gx, gy, grid(s, "d"),   levels=LD, colors=CD, linewidths=2.2)
    cx = a.contour(gx, gy, grid(s, "dxy"), levels=LX, colors=CX, linewidths=2.2,
                   linestyles="dashed")
    taken = []
    a.clabel(cd, manual=place(cd,a,taken), inline=True, inline_spacing=12, fontsize=15, fmt="%.2f")
    a.clabel(cx, manual=place(cx,a,taken), inline=True, inline_spacing=12, fontsize=15, fmt="%.2f")
    a.contour(gx, gy, grid(s.assign(_diff=s.dxy-s.d), "_diff"), levels=[0.], colors="k", linewidths=3.4)
    a.plot(s.delta, s.U, ".", ms=3.4, color="0.25", zorder=5)
    Us = np.sort(s.U.unique())
    rd = [peak(g.delta.values, g.dtot_N.values)
          for g in (s[s.U==U].sort_values("delta") for U in Us)]
    a.plot(rd, Us, "o", ms=9.5, mfc="none", mec="w", mew=3.4, zorder=6)
    a.plot(rd, Us, "o", ms=9.5, mfc="none", mec="k", mew=1.8, zorder=7)
    a.set_title(rf"$L={L}$", fontsize=25, pad=12)
    a.set_xlabel(r"anisotropy  $\delta$", labelpad=8)
    a.set_xticks([.2,.3,.4,.5,.6])
    a.tick_params(direction="out", top=False, right=False, length=6, width=1.4)

ax[0].set_ylabel(r"$U/t$", labelpad=8)
H = [plt.Line2D([],[],color=CD,ls="-",lw=2.2), plt.Line2D([],[],color=CX,ls="--",lw=2.2),
     plt.Line2D([],[],color="k",ls="-",lw=3.4),
     plt.Line2D([],[],color="k",ls="none",marker="o",ms=9.5,mfc="none",mew=1.8)]
fig.legend(H, [r"$d_{x^2-y^2}$", r"$d_{xy}$", r"$d_{xy}=d_{x^2-y^2}$",
               r"ridge of $\Delta_{\rm tot}$"], loc="lower center", ncol=4, frameon=False,
           bbox_to_anchor=(.5,-.21), handlelength=2.4, columnspacing=2.6)
cb = fig.colorbar(im, ax=ax.tolist(), pad=.014, fraction=0.060/len(LS))
cb.set_label(r"$\Delta_{\rm tot}/N$", fontsize=24, labelpad=12)
cb.ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%.3f"))
cb.locator = matplotlib.ticker.MaxNLocator(6); cb.update_ticks()
cb.ax.tick_params(labelsize=17)
plt.show()
