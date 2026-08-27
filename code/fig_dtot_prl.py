"""Delta_tot with the PRL normalisation, and the Delta n(k) difference maps.

PRL Eq. (2) and the reference plotting code define

    Delta_tot = sum_k |n_up(k) - n_dn(k)| / L^2

i.e. the momentum sum divided by the number of sites. The absolute value is
INSIDE the sum, so Delta_tot is nonzero even though the net magnetisation
sum_k [n_up(k) - n_dn(k)] vanishes: compensated but spin-split.

Our earlier version left out the 1/L^2, which put the numbers an order of two of
magnitude above the PRL's scale. With it, L=14 half filling runs 0 to 0.078
against their 0 to 0.50 at L=16.

Panel (a) is the direct analogue of PRL Fig. 1(a). The white region is real:
the (U, delta) grid holds 19 of 35 cells, only delta=0.2 has a full U column,
and delta=0 and 0.4 have one point each, so the convex hull is a diamond and
griddata returns NaN outside it. Open circles mark the measured points.
"""
import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt
from scipy.interpolate import griddata
from gridinterp import reggrid

mpl.rcParams.update({
    "font.family": "serif", "font.serif": ["DejaVu Serif"],
    "mathtext.fontset": "dejavuserif", "axes.linewidth": 1.2,
    "xtick.direction": "out", "ytick.direction": "out",
    "xtick.top": False, "ytick.right": False,
    "axes.labelsize": 17, "xtick.labelsize": 13, "ytick.labelsize": 13,
    "figure.dpi": 130})

D = "/Users/liujiaxin/Desktop/checkerboard-altermagnet/data"
d = pd.read_csv(f"{D}/fortran_dnk.csv")
d["dtot_N"] = d.delta_tot / d.L ** 2          # <-- the PRL 1/L^2
z = np.load(f"{D}/fortran_dnk_maps.npz", allow_pickle=True)
maps, meta = z["maps"], z["meta"]             # meta: L, n, U, delta

# ---------------- (1) Delta_tot heatmaps ----------------
fig, ax = plt.subplots(1, 2, figsize=(14.4, 5.6))

h = d[(d.L == 14) & np.isclose(d.n, 1.0)]
gx, gy = np.meshgrid(np.linspace(0, 0.4, 300), np.linspace(0, 5, 300))
gi = reggrid(h, 'delta', 'U', 'dtot_N', gx, gy)
im = ax[0].contourf(gx, gy, gi, levels=100, cmap="jet", extend="both")
cb = fig.colorbar(im, ax=ax[0], pad=0.02)
cb.set_label(r"$\Delta_{tot}$", fontsize=17, rotation=0, labelpad=-38, y=1.10)
cb.ax.tick_params(labelsize=11)
ax[0].plot(h.delta, h.U, "o", ms=5, mfc="none", mec="k", mew=1.0)
ax[0].set_xlabel(r"$t_A$" "  " r"(our $\delta$)"); ax[0].set_ylabel(r"$U/t$")
ax[0].set_title(r"(a) $L=14$, half filling", fontsize=17)

g10 = d[(d.L == 10) & (d.U == 4.0)]
gx2, gy2 = np.meshgrid(np.linspace(0.5, 0.98, 300), np.linspace(0, 0.4, 300))
gi2 = reggrid(g10, 'n', 'delta', 'dtot_N', gx2, gy2)
im2 = ax[1].contourf(gx2, gy2, gi2, levels=100, cmap="jet", extend="both")
cb2 = fig.colorbar(im2, ax=ax[1], pad=0.02)
cb2.set_label(r"$\Delta_{tot}$", fontsize=17, rotation=0, labelpad=-38, y=1.10)
cb2.ax.tick_params(labelsize=11)
ax[1].plot(g10.n, g10.delta, "o", ms=5, mfc="none", mec="k", mew=1.0)
ax[1].set_xlabel(r"filling  $n$"); ax[1].set_ylabel(r"anisotropy  $\delta$")
ax[1].set_title(r"(b) $L=10$, $U=4$", fontsize=17)
fig.tight_layout()
fig.savefig(f"{D}/../figs/fig_dtot_prl.png", dpi=600, bbox_inches="tight",
            facecolor="white")

# ---------------- (2) the difference over the Brillouin zone ----------------
want = [(14, 1.0, 4.0, 0.0), (14, 1.0, 4.0, 0.2), (14, 1.0, 4.0, 0.4)]
sel = []
for L_, n_, U_, dl_ in want:
    hit = np.where((meta[:, 0] == L_) & np.isclose(meta[:, 1], n_) &
                   np.isclose(meta[:, 2], U_) & np.isclose(meta[:, 3], dl_))[0]
    if len(hit): sel.append((hit[0], dl_))

fig, ax = plt.subplots(1, len(sel), figsize=(5.0 * len(sel), 4.8))
# one shared symmetric scale, so the delta = 0 panel reads as genuinely blank
# rather than being stretched to fill its own range
vmax = max(np.abs(maps[i][:, 2]).max() for i, _ in sel)
for j, (i, dl_) in enumerate(sel):
    a = ax[j] if len(sel) > 1 else ax
    kx, ky, dn = maps[i][:, 0], maps[i][:, 1], maps[i][:, 2]
    # tolerance matters: the file writes pi as 3.141592741 (float32), larger
    # than numpy's float64 pi, so a bare > np.pi sends the +pi points to -pi
    kxf = np.where(kx > np.pi + 1e-6, kx - 2 * np.pi, kx)
    kyf = np.where(ky > np.pi + 1e-6, ky - 2 * np.pi, ky)
    gx, gy = np.meshgrid(np.linspace(-np.pi, np.pi, 260),
                         np.linspace(-np.pi, np.pi, 260))
    gz = griddata((kxf, kyf), dn, (gx, gy), method="linear")
    im = a.contourf(gx, gy, gz, levels=80, cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    a.set_title(rf"$\delta = {dl_:g}$", fontsize=17)
    a.set_xlabel(r"$k_x$"); a.set_aspect("equal")
    a.set_xticks([-np.pi, 0, np.pi]); a.set_xticklabels([r"$-\pi$", "0", r"$\pi$"])
    a.set_yticks([-np.pi, 0, np.pi])
    a.set_yticklabels([r"$-\pi$", "0", r"$\pi$"] if j == 0 else [])
    if j == 0: a.set_ylabel(r"$k_y$")
cb = fig.colorbar(im, ax=ax, pad=0.02, fraction=0.03)
cb.set_label(r"$\Delta n(\mathbf{k}) = n_\uparrow(\mathbf{k}) - n_\downarrow(\mathbf{k})$",
             fontsize=15)
fig.suptitle(r"$L=14$, half filling, $U=4$", fontsize=17, y=1.02)
fig.savefig(f"{D}/../figs/fig_dnk_prl.png", dpi=600, bbox_inches="tight",
            facecolor="white")

print("Delta_tot = sum_k |dn(k)| / L^2   (PRL definition)\n")
print(h.pivot_table(index="U", columns="delta", values="dtot_N").round(4).to_string())
print(f"\nrange {h.dtot_N.min():.4f} to {h.dtot_N.max():.4f}"
      f"   (PRL Fig 1a: 0 to 0.50 at L=16)")
