"""Both d channels, EQUAL TIME only, over Delta_tot. Original blue palette.

Set LSEL below: [8, 10, 12] for the three-size figure, [12] for the single panel
that matches the left half of fig_dchannels_eqtime_vs_vertex.

Background is Blues, contours are dark navy and vermillion, the palette used in
the first version of this figure. The Blues ramp is truncated at 0.55 rather than
running to full saturation, so the navy contours keep something to read against.
At full saturation the dark end of the background and the navy lines were the
same value and the contours disappeared into it.

The equal-time numbers are the interaction-induced part: exactly zero at U=0 in
both channels, the same convention as the time-integrated vertex.
"""
import os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.interpolate import griddata
from check_crossover_ridge import peak
from gridinterp import reggrid

LSEL = [8, 10, 12]                 # <- [12] for the single panel
BGMAP, BGTOP = "Blues", 0.55       # background ramp and where to truncate it
CD, CX = "#08306B", "#D55E00"      # navy, vermillion

plt.rcParams.update({"font.family":"serif","font.serif":["DejaVu Serif"],
    "mathtext.fontset":"dejavuserif","font.size":16,"axes.labelsize":22,
    "legend.fontsize":19,"xtick.labelsize":18,"ytick.labelsize":18,
    "axes.linewidth":1.4,"figure.dpi":300,"savefig.bbox":"tight",
    "savefig.pad_inches":0.05})

HERE = os.path.dirname(os.path.abspath(__file__))
D    = os.path.join(HERE, "..", "data")
d = pd.read_csv(os.path.join(D, "fortran_L8_L10_L12_dedup.csv"))
d = d[np.isclose(d.n, 1.0) & (d.U >= 2) & d.L.isin(LSEL)]
LS = sorted(d.L.unique())

CMAP = LinearSegmentedColormap.from_list(
    "bg", plt.get_cmap(BGMAP)(np.linspace(0.0, BGTOP, 256)))
gx, gy = np.meshgrid(np.linspace(d.delta.min(), d.delta.max(), 400),
                     np.linspace(d.U.min(), d.U.max(), 400))
LEVBG = np.linspace(0, d.dtot_N.max(), 60)
LD, LX = [0.08, 0.16, 0.24], [0.12, 0.24, 0.36]


def grid(s, col):
    return reggrid(s, "delta", "U", col, gx, gy)


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


fig, ax = plt.subplots(1, len(LS), figsize=(6.2*len(LS), 5.6), sharex=True,
                       sharey=True, squeeze=False, gridspec_kw={"wspace":0.13})
ax = ax[0]
for a, L in zip(ax, LS):
    s = d[d.L == L]
    im = a.contourf(gx, gy, grid(s, "dtot_N"), levels=LEVBG,
                    cmap=CMAP, extend="both")
    cd = a.contour(gx, gy, grid(s, "d"),   levels=LD, colors=CD, linewidths=2.2)
    cx = a.contour(gx, gy, grid(s, "dxy"), levels=LX, colors=CX, linewidths=2.2,
                   linestyles="dashed")
    taken = []
    a.clabel(cd, manual=place(cd, a, taken), inline=True, inline_spacing=12,
             fontsize=15, fmt="%.2f")
    a.clabel(cx, manual=place(cx, a, taken), inline=True, inline_spacing=12,
             fontsize=15, fmt="%.2f")
    a.contour(gx, gy, grid(s.assign(_diff=s.dxy - s.d), "_diff"), levels=[0.0],
              colors="k", linewidths=3.4)
    a.plot(s.delta, s.U, ".", ms=3.4, color="0.25", zorder=5)
    Us = np.sort(s.U.unique())
    rd = [peak(g.delta.values, g.dtot_N.values)
          for g in (s[s.U == U].sort_values("delta") for U in Us)]
    a.plot(rd, Us, "o", ms=9.5, mfc="none", mec="w", mew=3.4, zorder=6)
    a.plot(rd, Us, "o", ms=9.5, mfc="none", mec="k", mew=1.8, zorder=7)
    a.set_title(rf"$L={L}$", fontsize=25, pad=12)
    a.set_xlabel(r"anisotropy  $\delta$", labelpad=8)
    a.set_xticks([0.2, 0.3, 0.4, 0.5, 0.6])
    a.tick_params(direction="out", top=False, right=False, length=6, width=1.4)

ax[0].set_ylabel(r"$U/t$", labelpad=8)
H = [plt.Line2D([], [], color=CD, ls="-",  lw=2.2),
     plt.Line2D([], [], color=CX, ls="--", lw=2.2),
     plt.Line2D([], [], color="k", ls="-", lw=3.4),
     plt.Line2D([], [], color="k", ls="none", marker="o", ms=9.5, mfc="none", mew=1.8)]
fig.legend(H, [r"$d_{x^2-y^2}$", r"$d_{xy}$", r"$d_{xy}=d_{x^2-y^2}$",
               r"ridge of $\Delta_{\rm tot}$"], loc="lower center", ncol=4,
           frameon=False, bbox_to_anchor=(0.5, -0.21), handlelength=2.4,
           columnspacing=2.6)
cb = fig.colorbar(im, ax=ax.tolist(), pad=0.014, fraction=0.060/len(LS))
cb.set_label(r"$\Delta_{\rm tot}/N$", fontsize=24, labelpad=12)
cb.ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%.3f"))
cb.locator = matplotlib.ticker.MaxNLocator(6); cb.update_ticks()
cb.ax.tick_params(labelsize=17)

tag = "".join(str(L) for L in LS)
out = os.path.join(HERE, "..", "manuscript", "figures", f"fig_dchannels_eqtime_L{tag}")
fig.savefig(out + ".pdf"); fig.savefig(out + ".png")
print(f"wrote {out}.png / .pdf   panels: L = {LS}")
