# --- Delta_tot heat map with the pairing lines on top, one panel per lattice size. ---
# Paste-and-run on 251.
#
# This is the "do the pairing lines follow the polarisation?" figure. Background is
# the magnetic order parameter Delta_tot/N in the (U, delta) plane; on top are the
# two d-channel contours and the white line where they are equal. If pairing tracks
# magnetism, the contours should lean the same way as the coloured structure.
#
# L = 8, 10, 12  full grid, delta 0.1-0.7 and U 2-5
# L = 14         Delta_tot from fortran_dnk.csv, pairing from fortran_pairing_all.csv.
#                Its scan is a cross, so only the sub-box delta 0.1-0.3, U 3-5 is a
#                complete rectangle. The panel keeps the same axes as the others and
#                simply draws nothing outside that box, which makes the coverage
#                difference visible rather than hidden.
#                Two U=2 cells (delta 0.1 and 0.3) are running now; when they land
#                the box becomes U 2-5 and this panel can be regenerated as-is.
# L = 16         two cells only (U=4, delta 0.2 and 0.3). Not enough for a panel.
import os, glob, numpy as np, pandas as pd, matplotlib.pyplot as plt, matplotlib
import matplotlib.patheffects as pe
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
    "font.size":14,"axes.labelsize":18,"legend.fontsize":14,
    "xtick.labelsize":14,"ytick.labelsize":14,"axes.linewidth":1.4})

F = pd.read_csv(find("fortran_L8_L10_L12_dedup.csv"))
F = F[np.isclose(F.n, 1.0) & (F.U >= 2)]
M14 = pd.read_csv(find("fortran_dnk.csv"))
M14 = M14[np.isclose(M14.n, 1.0) & (M14.L == 14)][["U","delta","dtot_N"]]
P14 = pd.read_csv(find("fortran_pairing_all.csv"))
P14 = P14[np.isclose(P14.n, 1.0) & (P14.L == 14)][["U","delta","d","dxy"]]
D14 = M14.merge(P14, on=["U","delta"])

def cells(L):
    """(dataframe, delta range, U range) for one lattice size."""
    if L <= 12:
        s = F[F.L == L]
        return s, (0.1, 0.7), (2.0, 5.0)
    s = D14[(D14.delta.between(0.1, 0.3)) & (D14.U >= 2)]
    piv = s.pivot_table(index="U", columns="delta", values="d")
    piv = piv.dropna(axis=0, how="any")          # keep only complete U rows
    s = s[s.U.isin(piv.index)]
    return s, (float(s.delta.min()), float(s.delta.max())), \
              (float(s.U.min()), float(s.U.max()))

LS   = [8, 10, 12, 14]
XLIM, YLIM = (0.1, 0.7), (2.0, 5.0)
LEVBG = np.linspace(0, F.dtot_N.max(), 60)
LD, LX = [.08, .16, .24], [.12, .24, .36]
HALO = [pe.withStroke(linewidth=2.6, foreground="w")]

def grid(s, col, gx, gy):
    """Bilinear on the measured rectangle. NOT scipy.griddata: that triangulates,
    and on a rectangular lattice the triangulation is neither unique nor stable."""
    p = s.pivot_table(index="U", columns="delta", values=col)
    assert not p.isna().any().any(), f"{col}: grid has holes"
    fn = RegularGridInterpolator((p.index.values, p.columns.values), p.values,
                                 method="linear", bounds_error=False, fill_value=None)
    return fn(np.stack([gy.ravel(), gx.ravel()], -1)).reshape(gx.shape)

def halo(cs):
    if hasattr(cs, "set_path_effects"): cs.set_path_effects(HALO)
    else:
        for c in cs.collections: c.set_path_effects(HALO)

fig, ax = plt.subplots(1, len(LS), figsize=(4.9*len(LS), 5.0), sharex=True,
                       sharey=True, gridspec_kw={"wspace":.08})
for a, L in zip(ax, LS):
    s, (x0, x1), (y0, y1) = cells(L)
    gx, gy = np.meshgrid(np.linspace(x0, x1, 400), np.linspace(y0, y1, 400))
    im = a.contourf(gx, gy, grid(s, "dtot_N", gx, gy), levels=LEVBG,
                    cmap="jet", extend="both")
    cd = a.contour(gx, gy, grid(s, "d", gx, gy),   levels=LD, colors="k", linewidths=2.0)
    cx = a.contour(gx, gy, grid(s, "dxy", gx, gy), levels=LX, colors="k", linewidths=2.0,
                   linestyles="dashed")
    halo(cd); halo(cx)
    ceq = grid(s.assign(_v=s.dxy - s.d), "_v", gx, gy)
    a.contour(gx, gy, ceq, levels=[0.], colors="k", linewidths=5.5)
    a.contour(gx, gy, ceq, levels=[0.], colors="w", linewidths=2.8)
    if (x0, x1) != XLIM or (y0, y1) != YLIM:
        a.add_patch(plt.Rectangle((x0, y0), x1-x0, y1-y0, fill=False,
                                  ec="k", lw=1.2, ls=":", zorder=6))
        a.set_title(rf"$L={L}$   (partial scan)", fontsize=17, pad=10)
    else:
        a.set_title(rf"$L={L}$", fontsize=18, pad=10)
    a.set_xlabel(r"anisotropy  $\delta$", labelpad=7)
    a.set_xlim(*XLIM); a.set_ylim(*YLIM)
    a.set_xticks([.2,.3,.4,.5,.6]); a.set_yticks([2,3,4,5])
    a.tick_params(direction="out", top=False, right=False, length=6, width=1.4)
ax[0].set_ylabel(r"$U/t$", labelpad=8)

H = [plt.Line2D([],[],color="k",ls="-",lw=2.0), plt.Line2D([],[],color="k",ls="--",lw=2.0),
     plt.Line2D([],[],color="w",ls="-",lw=2.8,
                path_effects=[pe.withStroke(linewidth=5.5, foreground="k")])]
fig.legend(H, [r"$d_{x^2-y^2}$", r"$d_{xy}$", r"$d_{xy}=d_{x^2-y^2}$"],
           loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(.5, -.17),
           handlelength=2.2, columnspacing=2.4, fontsize=17)
cb = fig.colorbar(im, ax=ax.tolist(), pad=.012, fraction=.016)
cb.set_label(r"$\Delta_{\rm tot}/N$", fontsize=19, labelpad=9)
cb.ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%.3f"))
cb.ax.tick_params(labelsize=13)
fig.suptitle("equal time, half filling.  background = magnetic polarisation,"
             " lines = pairing", fontsize=17, y=1.02)
plt.show()
