"""Equal-time vs time-integrated vertex, both d channels, over Delta_tot. L=12.

Left panel  EQUAL TIME. The Fortran driver computes only equal-time pair
            correlations (see checkerboard_eqtime.py), so this is what
            fortran_L8_L10_L12_dedup.csv carries in its d and dxy columns.

Right panel THE VERTEX, time integrated. Full correlator minus its uncorrelated
            Wick piece, at beta = 32, N_w = 500, tau window 0.8, six seeds.
            chi_L12_vertex_summary.csv, built by analyze_chi_L12.py.

Both panels are L=12 and share the (U, delta) grid, the Delta_tot/N background and
its colour scale. Contour levels differ between panels because the two quantities
differ by roughly a factor of twenty in magnitude, so each panel is labelled with
its own values.

Only L=12 has the vertex on a full (U, delta) grid. L=8 and L=10 were run at U=4
alone, so the three-panel layout of the equal-time figure has no vertex counterpart.

The equal-time grid reaches U=5 and the vertex stops at U=4.5, so both panels are
trimmed to U <= 4.5 and compare like with like.

U=0 is excluded from the map. The vertex is measured there and vanishes, 5.6e-13
and 3.5e-13 in the two channels, but there is no data between U=0 and U=2 and
interpolating across that gap would invent structure.
"""
import os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.interpolate import griddata
from gridinterp import reggrid

plt.rcParams.update({"font.family":"serif","font.serif":["DejaVu Serif"],
    "mathtext.fontset":"dejavuserif","font.size":16,"axes.labelsize":22,
    "legend.fontsize":19,"xtick.labelsize":18,"ytick.labelsize":18,
    "axes.linewidth":1.4,"figure.dpi":300,"savefig.bbox":"tight",
    "savefig.pad_inches":0.05})

HERE = os.path.dirname(os.path.abspath(__file__))
D    = os.path.join(HERE, "..", "data")
CD, CX = "#0072B2", "#D55E00"
GREY = LinearSegmentedColormap.from_list(
    "grey72", plt.get_cmap("Greys")(np.linspace(0.0, 0.72, 256)))

f = pd.read_csv(os.path.join(D, "fortran_L8_L10_L12_dedup.csv"))
f = f[np.isclose(f.n, 1.0) & (f.L == 12) & (f.U >= 2)]
v = pd.read_csv(os.path.join(D, "chi_L12_vertex_summary.csv"))
v = v[v.U >= 2]
UMAX = min(f.U.max(), v.U.max())          # equal time reaches U=5, the vertex 4.5
f = f[f.U <= UMAX]
bg = f[["U", "delta", "dtot_N"]]                  # Delta_tot is static either way

# each panel: (data, d levels, dxy levels, title)
PAN = [(f, [0.08, 0.16, 0.24], [0.12, 0.24, 0.36], "equal time"),
       (v, [1.0,  2.0,  3.0 ], [1.0,  3.0,  5.0 ], "vertex, time integrated")]

gx, gy = np.meshgrid(np.linspace(0.1, 0.7, 400), np.linspace(2.0, UMAX, 400))
LEVBG = np.linspace(0, bg.dtot_N.max(), 60)


def grid(s, col, col_x="delta", col_y="U"):
    return reggrid(s, col_x, col_y, col, gx, gy)


def place(cs, ax, taken):
    """Label points scored by distance from the frame and from labels already set."""
    x0, x1 = ax.get_xlim(); y0, y1 = ax.get_ylim(); out = []
    for segs in cs.allsegs:
        best, score = None, -9.0
        for sg in segs:
            if len(sg) < 2: continue
            u = (sg[:,0]-x0)/(x1-x0); w = (sg[:,1]-y0)/(y1-y0)
            sc = np.minimum.reduce([u, 1-u, w, 1-w])
            for tu, tw in taken:
                sc = np.minimum(sc, 0.55*np.hypot(u-tu, w-tw))
            k = int(np.argmax(sc))
            if sc[k] > score: score, best = sc[k], (tuple(sg[k]), (u[k], w[k]))
        if best is not None and score > 0.045:
            out.append(best[0]); taken.append(best[1])
    return out


BGZ = grid(bg, "dtot_N")
assert not np.isnan(BGZ).any(), "hull hole in the background"

fig, ax = plt.subplots(1, 2, figsize=(6.2*2, 5.6), sharex=True, sharey=True,
                       gridspec_kw={"wspace":0.13})
for a, (s, ld, lx, ttl) in zip(ax, PAN):
    im = a.contourf(gx, gy, BGZ, levels=LEVBG, cmap=GREY, extend="both")
    cd = a.contour(gx, gy, grid(s, "d"),   levels=ld, colors=CD, linewidths=2.2)
    cx = a.contour(gx, gy, grid(s, "dxy"), levels=lx, colors=CX, linewidths=2.2,
                   linestyles="dashed")
    taken = []
    fmt = "%.2f" if max(ld) < 1 else "%.1f"
    a.clabel(cd, manual=place(cd, a, taken), inline=True, inline_spacing=12,
             fontsize=15, fmt=fmt)
    a.clabel(cx, manual=place(cx, a, taken), inline=True, inline_spacing=12,
             fontsize=15, fmt=fmt)
    a.contour(gx, gy, grid(s.assign(_diff=s.dxy - s.d), "_diff"), levels=[0.0],
              colors="k", linewidths=3.4)
    a.plot(s.delta, s.U, ".", ms=3.4, color="0.35", zorder=5)
    a.set_title(ttl, fontsize=25, pad=12)
    a.set_xlabel(r"anisotropy  $\delta$", labelpad=8)
    a.set_xticks([0.2, 0.3, 0.4, 0.5, 0.6])
    a.set_ylim(2.0, UMAX); a.set_xlim(0.1, 0.7)
    a.tick_params(direction="out", top=False, right=False, length=6, width=1.4)

ax[0].set_ylabel(r"$U/t$", labelpad=8)
H = [plt.Line2D([], [], color=CD, ls="-",  lw=2.2),
     plt.Line2D([], [], color=CX, ls="--", lw=2.2),
     plt.Line2D([], [], color="k", ls="-", lw=3.4)]
fig.legend(H, [r"$d_{x^2-y^2}$", r"$d_{xy}$", r"$d_{xy}=d_{x^2-y^2}$"],
           loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, -0.21),
           handlelength=2.4, columnspacing=2.6)
cb = fig.colorbar(im, ax=ax.tolist(), pad=0.014, fraction=0.020)
cb.set_label(r"$\Delta_{\rm tot}/N$", fontsize=24, labelpad=12)
cb.ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%.3f"))
cb.ax.tick_params(labelsize=17)
fig.suptitle(r"$L=12$, half filling", fontsize=25, y=1.06)

out = os.path.join(HERE, "..", "manuscript", "figures", "fig_dchannels_eqtime_vs_vertex")
fig.savefig(out + ".pdf"); fig.savefig(out + ".png")
print(f"wrote {out}.png / .pdf")
for s, _, _, ttl in PAN:
    print(f"  {ttl:>24}: d {s.d.min():+.3f}..{s.d.max():+.3f}   "
          f"dxy {s.dxy.min():+.3f}..{s.dxy.max():+.3f}")
