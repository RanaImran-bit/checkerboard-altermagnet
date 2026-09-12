"""Delta_tot over the complete L=12 half-filling grid, as a smooth colormap.

7 U x 7 delta = 49 cells, every one measured, so shading='gouraud' interpolates
between real data points rather than across gaps. This is why the map is only
honest now: on the partial grid it would have filled empty cells with invented
colour.

Same data in both panels, axes swapped, so either control parameter can be read
along the horizontal:

  (a) delta on x -- follow a row to see what changes as delta changes
  (b) U     on x -- follow a row to see what changes as U changes

Delta_tot = sum_k |n_up(k) - n_dn(k)| / L^2, PRL Eq. (2). Absolute value inside
the sum, so it is finite despite the net magnetisation vanishing.

U = 0 is DROPPED from the map. It is a control that is identically zero at every
delta, and the U grid is unevenly spaced (0, 2, 3, 3.5, 4, 4.5, 5), so the gap
between 0 and 2 is four times the spacing used higher up. Gouraud smooths across
it, which put a wide interpolated gradient with no measured points inside over
roughly the bottom 40% of the axis. Keeping U >= 2 leaves only the densely and
evenly sampled region, where the shading is interpolating between neighbours
that are actually close together.
"""
import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt

mpl.rcParams.update({
    "font.family": "serif", "font.serif": ["DejaVu Serif"],
    "mathtext.fontset": "dejavuserif", "axes.linewidth": 1.1,
    "xtick.direction": "out", "ytick.direction": "out",
    "axes.labelsize": 16, "xtick.labelsize": 13, "ytick.labelsize": 13,
    "figure.dpi": 130})

CMAP = "inferno"          # one-line switch: "jet" matches the older Delta_tot figures

D = "/Users/liujiaxin/Desktop/checkerboard-altermagnet"
d = pd.read_csv(f"{D}/data/fortran_L12_all49.csv")
d = d[np.isclose(d.n, 1.0) & (d.U >= 2)]      # see docstring: U = 0 dropped
P = d.pivot_table(index="U", columns="delta", values="dtot_N")
US, DS, Z = P.index.values, P.columns.values, P.values
assert not np.isnan(Z).any(), "grid has a hole, gouraud would invent colour there"
vmax = Z.max()

fig, ax = plt.subplots(1, 2, figsize=(13.6, 5.0))

im0 = ax[0].pcolormesh(DS, US, Z, cmap=CMAP, shading="gouraud",
                       vmin=0, vmax=vmax)
ax[0].set_xlabel(r"anisotropy  $\delta$"); ax[0].set_ylabel(r"$U/t$")
ax[0].set_title("(a)", loc="left", fontsize=15)
cb0 = fig.colorbar(im0, ax=ax[0], pad=0.02)
cb0.set_label(r"$\Delta_{\mathrm{tot}}$", fontsize=15)

im1 = ax[1].pcolormesh(US, DS, Z.T, cmap=CMAP, shading="gouraud",
                       vmin=0, vmax=vmax)
ax[1].set_xlabel(r"$U/t$"); ax[1].set_ylabel(r"anisotropy  $\delta$")
ax[1].set_title("(b)", loc="left", fontsize=15)
cb1 = fig.colorbar(im1, ax=ax[1], pad=0.02)
cb1.set_label(r"$\Delta_{\mathrm{tot}}$", fontsize=15)

fig.tight_layout()
fig.savefig(f"{D}/figs/fig_dtot_49_map.png", dpi=600, bbox_inches="tight",
            facecolor="white")

print(f"grid {Z.shape[0]} U x {Z.shape[1]} delta, no missing cells")
print(f"Delta_tot range 0 to {vmax:.4f}")
print("\nridge -- delta of maximum Delta_tot at each U:")
for i, U in enumerate(US):
    j = int(np.argmax(Z[i]))
    print(f"   U={U:g}: delta={DS[j]:.1f}  ({Z[i, j]:.4f})")
