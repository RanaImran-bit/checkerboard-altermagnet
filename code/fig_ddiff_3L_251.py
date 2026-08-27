# --- which d channel wins, time-integrated vertex, three sizes. Paste-and-run on 251. ---
# Red: d_xy is the more enhanced channel. Blue: d_x2-y2 is. Black line: they are equal.
# One shared symmetric colour scale across the three panels, so the panels show both
# where the boundary sits AND that d_xy strengthens with system size (L=8 is paler
# because its signal is genuinely ~2.3x weaker, not because of the scaling).
#
# Data: chi_vertex_summary_3L.csv, built by build_chi_summary_3L.py from the b32 runs
# collected off 251. All three sizes share beta=32, tau window 0.8, N_w=500, 6 seeds,
# and a complete 7x7 (U, delta) grid.
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
    "font.size":15,"axes.labelsize":20,"xtick.labelsize":16,"ytick.labelsize":16,
    "axes.linewidth":1.4})

S = pd.read_csv(find("chi_vertex_summary_3L.csv")); S = S[S.U >= 2]
LS = sorted(S.L.unique()); UMAX = S.U.max()
gx, gy = np.meshgrid(np.linspace(.1,.7,400), np.linspace(2.,UMAX,400))

def grid(s, col):
    """Bilinear on the measured rectangle. NOT griddata: that triangulates, and on
    a rectangular lattice the triangulation is neither unique nor stable."""
    p = s.pivot_table(index="U", columns="delta", values=col)
    assert not p.isna().any().any(), "grid has holes"
    fn = RegularGridInterpolator((p.index.values, p.columns.values), p.values,
                                 method="linear", bounds_error=False, fill_value=None)
    return fn(np.stack([gy.ravel(), gx.ravel()], -1)).reshape(gx.shape)

def nice_step(span, target=6):
    raw = span / target; mag = 10.0 ** np.floor(np.log10(raw))
    for m in (1, 2, 2.5, 5, 10):
        if raw <= m * mag: return m * mag
    return 10 * mag

Z = {L: grid(S[S.L == L].assign(_v=lambda t: t.dxy - t.d), "_v") for L in LS}
M = max(np.abs(z).max() for z in Z.values())            # one scale for all panels
step = nice_step(2 * M)
ticks = np.arange(np.ceil(-M/step)*step, M + step/2, step)
ticks = ticks[(ticks >= -M) & (ticks <= M)]

fig, ax = plt.subplots(1, len(LS), figsize=(5.6*len(LS), 5.2), sharex=True,
                       sharey=True, squeeze=False, gridspec_kw={"wspace":.09})
ax = ax[0]
for a, L in zip(ax, LS):
    im = a.contourf(gx, gy, Z[L], levels=np.linspace(-M, M, 81), cmap="RdBu_r",
                    vmin=-M, vmax=M)
    a.contour(gx, gy, Z[L], levels=[0.], colors="k", linewidths=3.2)
    a.set_title(rf"$L={int(L)}$", fontsize=21, pad=10)
    a.set_xlabel(r"anisotropy  $\delta$", labelpad=8)
    a.set_xticks([.2,.3,.4,.5,.6]); a.set_xlim(.1,.7); a.set_ylim(2., UMAX)
    a.tick_params(direction="out", top=False, right=False, length=6, width=1.4)
ax[0].set_ylabel(r"$U/t$", labelpad=8)
cb = fig.colorbar(im, ax=ax.tolist(), pad=.014, fraction=.019, ticks=ticks)
cb.set_label(r"$d_{xy}-d_{x^2-y^2}$", fontsize=20, labelpad=10)
cb.ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%.0f"))
cb.ax.tick_params(labelsize=15)
fig.suptitle(r"time-integrated vertex, half filling."
             r"   black line: $d_{xy}=d_{x^2-y^2}$", fontsize=18, y=1.00)
plt.show()
