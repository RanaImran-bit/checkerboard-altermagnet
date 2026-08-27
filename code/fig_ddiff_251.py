# --- heat map of d_xy - d_x2-y2, which channel wins. L=12. ---
# Red: d_xy is the more enhanced channel. Blue: d_x2-y2 is. White and the black
# line: the two are equal. The colour scale is pinned symmetric about zero so white
# really means zero rather than the middle of the data.
#
# This is the map that answers whether anisotropy selects d_xy, which the sum
# cannot: adding the channels throws away exactly the competition between them.
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

def nice_step(span, target=6):
    """Round step giving roughly `target` intervals across `span`."""
    raw = span / target
    mag = 10.0 ** np.floor(np.log10(raw))
    for m in (1, 2, 2.5, 5, 10):
        if raw <= m * mag: return m * mag
    return 10 * mag

def limits(Z, symmetric=False, target=6):
    """Colour limits plus evenly spaced round tick values that lie inside them.

    The ticks were uneven before because extend="both" drew triangular ends while
    vmin and vmax were the raw data min and max, so the first and last ticks fell
    outside the bar and were squashed against its ends. Here the bar still spans
    exactly the data, wasting no colour, and the ticks are whole multiples of a
    round step taken from inside that span. Since contourf levels are evenly
    spaced in value, tick position is proportional to value and the spacing comes
    out even."""
    lo, hi = (-np.abs(Z).max(), np.abs(Z).max()) if symmetric else (Z.min(), Z.max())
    step = nice_step(hi - lo, target)
    ticks = np.arange(np.ceil(lo / step) * step, hi + step / 2, step)
    return lo, hi, step, ticks[(ticks >= lo) & (ticks <= hi)]

def panel(fig, a, Z, cmap, lab, symmetric=False, zeroline=False):
    lo, hi, step, ticks = limits(Z, symmetric)
    im = a.contourf(gx, gy, Z, levels=np.linspace(lo, hi, 81), cmap=cmap,
                    vmin=lo, vmax=hi)
    if zeroline:
        a.contour(gx, gy, Z, levels=[0.], colors="k", linewidths=3.0)
    cb = fig.colorbar(im, ax=a, pad=.025, fraction=.05, ticks=ticks)
    dec = max(0, int(np.ceil(-np.log10(step))) + 1)
    cb.ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter(f"%.{dec}f"))
    cb.ax.tick_params(labelsize=14)
    cb.set_label(lab, fontsize=19, labelpad=8)
    return im

COLS = [(f, "equal time"), (v, "time integrated")]
fig, ax = plt.subplots(1, 2, figsize=(13.0, 5.2), sharex=True, sharey=True,
                       gridspec_kw={"wspace":.34})
for a, (s, lab) in zip(ax, COLS):
    panel(fig, a, grid(s.assign(_v=s.dxy - s.d), "_v"), "RdBu_r",
          r"$d_{xy}-d_{x^2-y^2}$", symmetric=True, zeroline=True)
    a.set_title(f"{lab}  ($L=12$, half filling)", fontsize=18, pad=10)
    a.set_xlabel(r"anisotropy  $\delta$", labelpad=8)
    a.set_xticks([.2,.3,.4,.5,.6]); a.set_xlim(.1,.7); a.set_ylim(2., UMAX)
    a.tick_params(direction="out", top=False, right=False, length=6, width=1.4)
ax[0].set_ylabel(r"$U/t$", labelpad=8)
fig.suptitle(r"black line: $d_{xy}=d_{x^2-y^2}$", fontsize=17, y=1.02)
plt.show()
