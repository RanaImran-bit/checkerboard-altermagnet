"""The two d channels together: which is attractive, which repulsive, and where they cross.

Row 1  both channels against delta, one panel per L. The zero line separates
       attractive (above) from repulsive (below), and the crossing is visible.
Row 2  the difference d_xy - d_x2-y2 over the (delta, U) plane, with the zero
       contour drawn heavy. That contour IS the crossover line.

Diverging data needs a diverging map with a neutral midpoint, not a rainbow:
red where d_xy wins, blue where d_x2-y2 wins, white where they are equal.
"""
import os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from gridinterp import reggrid

plt.rcParams.update({"font.family":"serif","font.serif":["DejaVu Serif"],
    "mathtext.fontset":"dejavuserif","font.size":9,"axes.labelsize":11,
    "legend.fontsize":8.5,"xtick.labelsize":9,"ytick.labelsize":9,
    "axes.linewidth":0.9,"lines.linewidth":1.5,"lines.markersize":4.5,
    "figure.dpi":300,"savefig.bbox":"tight","savefig.pad_inches":0.03})

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
d = pd.read_csv(os.path.join(D, "fortran_L8_L10_L12_dedup.csv"))
d = d[np.isclose(d.n, 1.0) & (d.U >= 2)]
LS = sorted(d.L.unique())
CD, CX = "#0072B2", "#D55E00"          # d_x2-y2 blue, d_xy vermillion

fig, ax = plt.subplots(2, len(LS), figsize=(4.3*len(LS), 7.2),
                       gridspec_kw={"hspace":0.34, "wspace":0.22})

# ---- row 1: both channels vs delta
Us = sorted(d.U.unique())
sh = np.linspace(0.35, 1.0, len(Us))
for j, L in enumerate(LS):
    a = ax[0, j]; s = d[d.L == L]
    for k, U in enumerate(Us):
        r = s[s.U == U].sort_values("delta")
        a.plot(r.delta, r.d,   "-o", color=CD, alpha=sh[k], mfc="white", mew=1.0)
        a.plot(r.delta, r.dxy, "-s", color=CX, alpha=sh[k], mfc="white", mew=1.0)
    a.axhline(0, color="0.35", lw=1.1)
    a.text(0.115, 0.012, "attractive", fontsize=7.5, color="0.35")
    a.text(0.115, -0.030, "repulsive", fontsize=7.5, color="0.35")
    a.set_title(rf"$L={L}$", fontsize=12)
    a.set_xlabel(r"anisotropy  $\delta$")
    if j == 0:
        a.set_ylabel(r"pairing vertex")
        a.plot([], [], "-o", color=CD, mfc="white", label=r"$d_{x^2-y^2}$")
        a.plot([], [], "-s", color=CX, mfc="white", label=r"$d_{xy}$")
        a.legend(frameon=False, loc="upper left")
        a.text(0.02, 0.62, "lighter = smaller $U$", transform=a.transAxes,
               fontsize=7.5, color="0.45")
    a.tick_params(direction="in", top=True, right=True)
    a.set_ylim(-0.06, 0.48)

# ---- row 2: the difference, with the crossover contour
gx, gy = np.meshgrid(np.linspace(d.delta.min(), d.delta.max(), 300),
                     np.linspace(d.U.min(), d.U.max(), 300))
vmax = np.abs(d.dxy - d.d).max()
lev = np.linspace(-vmax, vmax, 41)
for j, L in enumerate(LS):
    a = ax[1, j]; s = d[d.L == L]
    Z = reggrid(s.assign(_diff=s.dxy - s.d), "delta", "U", "_diff", gx, gy)
    im = a.contourf(gx, gy, Z, levels=lev, cmap="RdBu_r", extend="both")
    z0 = a.contour(gx, gy, Z, levels=[0.0], colors="k", linewidths=2.0)
    a.clabel(z0, inline=True, fontsize=8, fmt=lambda v: "equal")
    a.set_xlabel(r"anisotropy  $\delta$")
    if j == 0: a.set_ylabel(r"$U/t$")
    a.tick_params(direction="in", top=True, right=True, colors="0.25")
    a.set_title(rf"$L={L}$", fontsize=12)
cb = fig.colorbar(im, ax=ax[1, :].tolist(), pad=0.015, fraction=0.020)
cb.set_label(r"$d_{xy} - d_{x^2-y^2}$", fontsize=11)
cb.ax.tick_params(labelsize=8.5)
fig.text(0.5, 0.475, r"red: $d_{xy}$ dominates      blue: $d_{x^2-y^2}$ dominates"
         "      heavy line: the two are equal", ha="center", fontsize=9, color="0.3")

out = os.path.join(D, "..", "manuscript", "figures", "fig_d_channels_together")
fig.savefig(out + ".pdf"); fig.savefig(out + ".png")
print(f"wrote {out}.png / .pdf\n")
print(f"{'L':>3}{'U':>6}  crossover delta (linear interp between grid points)")
for L in LS:
    s = d[d.L == L]
    for U in Us:
        r = s[s.U == U].sort_values("delta")
        diff = (r.dxy - r.d).to_numpy(); dl = r.delta.to_numpy()
        k = np.where(np.diff(np.sign(diff)))[0]
        xc = (dl[k[0]] - diff[k[0]]*(dl[k[0]+1]-dl[k[0]])/(diff[k[0]+1]-diff[k[0]])
              ) if len(k) else np.nan
        print(f"{L:3d}{U:6.1f}   {xc:.3f}" if np.isfinite(xc) else f"{L:3d}{U:6.1f}   none")
