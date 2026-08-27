# --- both d channels over Delta_tot: equal time (top) vs time integrated (bottom),
#     L = 8, 10, 12. Paste-and-run on 251. ---
# Background: Delta_tot/N, one jet scale shared by all six panels, so the growth of
# the magnetic signal with L is visible rather than normalised away. Delta_tot is a
# static quantity, so the background is the same field in both rows; only the
# contours change.
#
# Contours black, told apart by line style, since jet spans blue through red and any
# coloured contour vanishes somewhere on it. Crossover is white with a black casing.
# Contour LEVELS are fixed across each row, so the sparser contours at L=8 in the
# bottom row are real: the vertex there is genuinely about 2.3x weaker than at L=12.
#
# Sources: fortran_L8_L10_L12_dedup.csv for equal time and Delta_tot; the raw
# eqtime_L*_U*.csv runs for the vertex, summarised here on the fly.
import os, glob, numpy as np, pandas as pd, matplotlib.pyplot as plt, matplotlib
import matplotlib.patheffects as pe
from scipy.interpolate import RegularGridInterpolator

CH   = ["son", "sext", "d", "dxy"]
DIRS = ["~/collected/eqt_b32", "~/results/chi_L12", "~/eqt_b32", "~/chi_L12",
        "../data/chi_fss_b32", "../data/chi_L12", "chi_fss_b32", "chi_L12"]

def find(name):
    for c in [name, f"../data/{name}", f"data/{name}",
              os.path.expanduser(f"~/results/figdata/{name}"),
              os.path.expanduser(f"~/results/{name}"),
              os.path.expanduser(f"~/checkerboard-altermagnet/data/{name}")]:
        if os.path.exists(c): return c
    h = glob.glob(os.path.expanduser(f"~/results/**/{name}"), recursive=True)
    if h: return h[0]
    raise SystemExit(f"cannot find {name}. Try:  find ~ -name {name}")

# ---- equal time and the Delta_tot background ----
F = pd.read_csv(find("fortran_L8_L10_L12_dedup.csv"))
F = F[np.isclose(F.n, 1.0) & (F.U >= 2)]

# ---- vertex, summarised from the raw runs ----
files = []
for d in DIRS:
    files += glob.glob(os.path.join(os.path.expanduser(d), "eqtime_L*_U*.csv"))
if not files:
    raise SystemExit("no eqtime_L*_U*.csv found. Try:  find ~ -name 'eqtime_L*_U*.csv'")
raw = pd.concat([pd.read_csv(f) for f in sorted(set(files))], ignore_index=True)
raw = raw.drop_duplicates(subset=["L", "U", "delta", "seed"])
for c, want in (("beta_proj", 32.0), ("tau_max", 0.8), ("nw", 500)):
    assert sorted(raw[c].unique()) == [want], f"{c} differs across files, do not merge"
V = (raw.groupby(["L", "U", "delta"])[[f"chi_{c}_vertex" for c in CH]].mean()
        .rename(columns=lambda c: c.replace("chi_", "").replace("_vertex", ""))
        .reset_index())
V = V[V.U >= 2]

LS   = sorted(set(F.L.unique()) & set(V.L.unique()))
UMAX = min(F.U.max(), V.U.max())
F, V = F[F.U <= UMAX], V[V.U <= UMAX]
print(f"sizes {LS}, U up to {UMAX}")

plt.rcParams.update({"font.family":"serif","mathtext.fontset":"dejavuserif",
    "font.size":15,"axes.labelsize":20,"xtick.labelsize":16,"ytick.labelsize":16,
    "axes.linewidth":1.4})
gx, gy = np.meshgrid(np.linspace(.1,.7,400), np.linspace(2.,UMAX,400))
LEVBG  = np.linspace(0, F.dtot_N.max(), 60)     # one Delta_tot scale for all panels
HALO   = [pe.withStroke(linewidth=2.6, foreground="w")]

# (data, d levels, dxy levels, row label)
ROWS = [(F, [.08, .16, .24], [.12, .24, .36], "equal time"),
        (V, [1., 2., 3., 4.], [1., 2., 4.],   "time integrated")]

def grid(s, col):
    """Bilinear on the measured rectangle. NOT scipy.griddata: that triangulates,
    and on a rectangular lattice the triangulation is neither unique nor stable."""
    p = s.pivot_table(index="U", columns="delta", values=col)
    assert not p.isna().any().any(), "grid has holes"
    fn = RegularGridInterpolator((p.index.values, p.columns.values), p.values,
                                 method="linear", bounds_error=False, fill_value=None)
    return fn(np.stack([gy.ravel(), gx.ravel()], -1)).reshape(gx.shape)

def halo(cs):
    """White outline on a contour set, across matplotlib versions."""
    if hasattr(cs, "set_path_effects"): cs.set_path_effects(HALO)
    else:
        for c in cs.collections: c.set_path_effects(HALO)

def place(cs, ax, taken):
    """Label points scored by distance from the frame and from labels already set."""
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
        if best is not None and score > .05: out.append(best[0]); taken.append(best[1])
    return out

fig, ax = plt.subplots(2, len(LS), figsize=(5.7*len(LS), 10.4), sharex=True,
                       sharey=True, gridspec_kw={"wspace":.09, "hspace":.16})
for i, (T, ld, lx, rlab) in enumerate(ROWS):
    fmt = "%.2f" if max(ld) < 1 else "%.0f"
    for j, L in enumerate(LS):
        a = ax[i, j]; s = T[T.L == L]; bgs = F[F.L == L]
        im = a.contourf(gx, gy, grid(bgs, "dtot_N"), levels=LEVBG, cmap="jet",
                        extend="both")
        cd = a.contour(gx, gy, grid(s, "d"),   levels=ld, colors="k", linewidths=2.2)
        cx = a.contour(gx, gy, grid(s, "dxy"), levels=lx, colors="k", linewidths=2.2,
                       linestyles="dashed")
        halo(cd); halo(cx)
        taken = []
        t1 = a.clabel(cd, manual=place(cd,a,taken), inline=True, inline_spacing=12,
                      fontsize=14, fmt=fmt)
        t2 = a.clabel(cx, manual=place(cx,a,taken), inline=True, inline_spacing=12,
                      fontsize=14, fmt=fmt)
        for t in list(t1) + list(t2): t.set_path_effects(HALO)
        ceq = grid(s.assign(_v=s.dxy - s.d), "_v")
        a.contour(gx, gy, ceq, levels=[0.], colors="k", linewidths=6.0)
        a.contour(gx, gy, ceq, levels=[0.], colors="w", linewidths=3.0)
        a.set_title(rf"{rlab},  $L={int(L)}$", fontsize=18, pad=10)
        a.set_xticks([.2,.3,.4,.5,.6]); a.set_xlim(.1,.7); a.set_ylim(2., UMAX)
        a.tick_params(direction="out", top=False, right=False, length=6, width=1.4)
    ax[i, 0].set_ylabel(r"$U/t$", labelpad=8)
for a in ax[1]: a.set_xlabel(r"anisotropy  $\delta$", labelpad=8)

H = [plt.Line2D([],[],color="k",ls="-",lw=2.2), plt.Line2D([],[],color="k",ls="--",lw=2.2),
     plt.Line2D([],[],color="w",ls="-",lw=3.0,
                path_effects=[pe.withStroke(linewidth=6.0, foreground="k")])]
fig.legend(H, [r"$d_{x^2-y^2}$", r"$d_{xy}$", r"$d_{xy}=d_{x^2-y^2}$"],
           loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(.5, -.035),
           handlelength=2.4, columnspacing=2.6, fontsize=19)
cb = fig.colorbar(im, ax=ax.ravel().tolist(), pad=.014, fraction=.017)
cb.set_label(r"$\Delta_{\rm tot}/N$", fontsize=21, labelpad=10)
cb.ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%.3f"))
cb.ax.tick_params(labelsize=15)
plt.show()
