"""Both d channels as contours over the (U, delta) plane, with Delta_tot beneath.

One panel per lattice size. The background is the magnetic order parameter
Delta_tot/N on a single grey scale shared by all three panels. Over it sit two
contour sets and their crossover:

    solid blue        d_x2-y2 pairing vertex
    dashed vermillion d_xy pairing vertex
    heavy black       where the two are equal
    open circles      the ridge of Delta_tot, its maximum in delta at each U

Both channels stay attractive across the whole grid. d_x2-y2 is positive in
126/126 cells; what delta does is let d_xy overtake it, not turn it repulsive.

The crossover and the ridge track each other closely (pooled r = +0.94) and the
offset between them shrinks with size, +0.096, +0.052, +0.021 at L = 8, 10, 12.
See check_crossover_ridge.py.

The background is grey rather than the jet of the earlier four-panel figures.
Two coloured contour sets have to be read on top of it, and any strongly hued
background steals from them. Grey is capped at 0.72 so the dark end never
approaches the black crossover line.

Fields are interpolated from the measured grid (delta in steps of 0.1, U in
steps of 0.5) with a cubic scheme, linear-filled at the hull. Grey dots mark the
points actually measured. Contour position between dots is interpolation, not
data.
"""
import os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.interpolate import griddata
from check_crossover_ridge import peak
from gridinterp import reggrid

plt.rcParams.update({"font.family":"serif","font.serif":["DejaVu Serif"],
    "mathtext.fontset":"dejavuserif","font.size":16,"axes.labelsize":22,
    "legend.fontsize":19,"xtick.labelsize":18,"ytick.labelsize":18,
    "axes.linewidth":1.4,"figure.dpi":300,"savefig.bbox":"tight",
    "savefig.pad_inches":0.05})

HERE = os.path.dirname(os.path.abspath(__file__))
D    = os.path.join(HERE, "..", "data")
d = pd.read_csv(os.path.join(D, "fortran_L8_L10_L12_dedup.csv"))
d = d[np.isclose(d.n, 1.0) & (d.U >= 2)]
LS = sorted(d.L.unique())
CD, CX = "#0072B2", "#D55E00"                      # blue, vermillion

GREY = LinearSegmentedColormap.from_list(          # white -> mid grey, never black
    "grey72", plt.get_cmap("Greys")(np.linspace(0.0, 0.72, 256)))

gx, gy = np.meshgrid(np.linspace(d.delta.min(), d.delta.max(), 400),
                     np.linspace(d.U.min(), d.U.max(), 400))
LEVBG = np.linspace(0, d.dtot_N.max(), 60)
LD = [0.08, 0.16, 0.24]                             # d_x2-y2 levels
LX = [0.12, 0.24, 0.36]                             # d_xy levels

def place(cs, ax, taken):
    """One label point per contour level, chosen deterministically.

    Walk every vertex of every segment and score it by how far it sits from the
    panel edge and from labels already placed, both in axes fractions. Take the
    best. This keeps labels off the frame, where they were being clipped, and
    apart from each other, where blue and orange run parallel near the crossover.
    """
    x0, x1 = ax.get_xlim(); y0, y1 = ax.get_ylim()
    out = []
    for segs in cs.allsegs:
        best, score = None, -9.0
        for sg in segs:
            if len(sg) < 2: continue
            u = (sg[:,0]-x0)/(x1-x0); v = (sg[:,1]-y0)/(y1-y0)
            sc = np.minimum.reduce([u, 1-u, v, 1-v])
            for tu, tv in taken:
                sc = np.minimum(sc, 0.55*np.hypot(u-tu, v-tv))
            k = int(np.argmax(sc))
            if sc[k] > score:
                score, best = sc[k], (tuple(sg[k]), (u[k], v[k]))
        if best is not None and score > 0.045:
            out.append(best[0]); taken.append(best[1])
    return out


def grid(s, col):
    """Cubic where it is defined, linear elsewhere, so the hull stays filled."""
    return reggrid(s, "delta", "U", col, gx, gy)

fig, ax = plt.subplots(1, len(LS), figsize=(6.2*len(LS), 5.6),
                       sharex=True, sharey=True, gridspec_kw={"wspace":0.13})
print(f"{'L':>4} {'d_x2-y2 range':>26} {'d_xy range':>26}   sign of d_x2-y2")
for j, L in enumerate(LS):
    a = ax[j]; s = d[d.L == L]
    BG = grid(s, "dtot_N")
    assert not np.isnan(BG).any(), f"L={L}: hull hole"
    im = a.contourf(gx, gy, BG, levels=LEVBG, cmap=GREY, extend="both")

    cd = a.contour(gx, gy, grid(s, "d"),   levels=LD, colors=CD,
                   linewidths=2.2, linestyles="solid")
    cx = a.contour(gx, gy, grid(s, "dxy"), levels=LX, colors=CX,
                   linewidths=2.2, linestyles="dashed")
    taken = []
    a.clabel(cd, manual=place(cd, a, taken), inline=True, inline_spacing=12,
             fontsize=15, fmt="%.2f")
    a.clabel(cx, manual=place(cx, a, taken), inline=True, inline_spacing=12,
             fontsize=15, fmt="%.2f")
    a.contour(gx, gy, grid(s.assign(_diff=s.dxy - s.d), "_diff"), levels=[0.0],
              colors="k", linewidths=3.4)                 # no inline label
    a.plot(s.delta, s.U, ".", ms=3.4, color="0.35", zorder=5)

    Us = np.sort(s.U.unique())                       # ridge of Delta_tot
    rd = [peak(g.delta.values, g.dtot_N.values)
          for g in (s[s.U == U].sort_values("delta") for U in Us)]
    a.plot(rd, Us, "o", ms=9.5, mfc="none", mec="w", mew=3.4, zorder=6)
    a.plot(rd, Us, "o", ms=9.5, mfc="none", mec="k", mew=1.8, zorder=7)

    a.set_title(rf"$L={L}$", fontsize=25, pad=12)
    a.set_xlabel(r"anisotropy  $\delta$", labelpad=8)
    a.set_xticks([0.2, 0.3, 0.4, 0.5, 0.6])          # endpoints dropped, they collide
    a.tick_params(direction="out", top=False, right=False, length=6, width=1.4)

ax[0].set_ylabel(r"$U/t$", labelpad=8)

H = [plt.Line2D([], [], color=CD, ls="-",  lw=2.2),
     plt.Line2D([], [], color=CX, ls="--", lw=2.2),
     plt.Line2D([], [], color="k", ls="-", lw=3.4),
     plt.Line2D([], [], color="k", ls="none", marker="o", ms=9.5,
                mfc="none", mew=1.8)]
Lb = [r"$d_{x^2-y^2}$", r"$d_{xy}$", r"$d_{xy}=d_{x^2-y^2}$",
      r"ridge of $\Delta_{\rm tot}$"]
fig.legend(H, Lb, loc="lower center", ncol=4, frameon=False,
           bbox_to_anchor=(0.5, -0.21), handlelength=2.4, columnspacing=2.6)

cb = fig.colorbar(im, ax=ax.tolist(), pad=0.014, fraction=0.020)
cb.set_label(r"$\Delta_{\rm tot}/N$", fontsize=24, labelpad=12)
cb.ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%.3f"))
cb.ax.tick_params(labelsize=17)

for L in LS:
    s = d[d.L == L]
    print(f"{L:>4} {s.d.min():>12.4f} .. {s.d.max():<12.4f}"
          f" {s.dxy.min():>12.4f} .. {s.dxy.max():<12.4f}"
          f"   {'all +' if (s.d>0).all() else 'mixed'}")

out = os.path.join(HERE, "..", "manuscript", "figures", "fig_dchannels_over_dtot")
fig.savefig(out + ".pdf"); fig.savefig(out + ".png")
print(f"\nwrote {out}.png / .pdf")
print(f"background Delta_tot/N: 0 .. {d.dtot_N.max():.4f}, one scale for all panels")
