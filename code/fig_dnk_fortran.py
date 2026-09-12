"""Delta_tot and the Delta n(k) maps, from the Fortran CPQMC output.

Delta_tot = sum_k |n_up(k) - n_dn(k)|, the altermagnetic order parameter of the
PRL. The absolute value is inside the sum, so it is nonzero even though the net
magnetisation vanishes: compensated but spin-split.

The Fortran runs use a SYMMETRY-BROKEN trial (wfup.txt and wfdn.txt are separate
determinants), so the spin degeneracy is lifted before the QMC starts and the
difference is physical. Our spin-restricted Python runs build both trials from
the same matrix, which is why n_up(k) = n_dn(k) identically there.

Two limits confirm the emergent-altermagnetism picture, and both come out right:
    U = 0, delta = 0.2  ->  Delta_tot = 0.0000   no interaction, no splitting
    U = 4, delta = 0    ->  Delta_tot = 0.1552   no anisotropy, no splitting
The splitting needs BOTH, exactly as the mean-field schematic requires both a
staggered moment M and a finite delta.
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
z = np.load(f"{D}/fortran_dnk_maps.npz", allow_pickle=True)
maps, meta = z["maps"], z["meta"]        # meta cols: L, n, U, delta

# ---------------- FIGURE 1: Delta_tot heatmaps ----------------
fig, ax = plt.subplots(1, 2, figsize=(14.4, 5.6))

# (a) (U, delta) at L = 14, half filling -- the PRL Fig. 1(a) grid
h = d[(d.L == 14) & np.isclose(d.n, 1.0)]
_gx, _gy = np.meshgrid(np.linspace(0, 0.4, 300), np.linspace(0, 5, 300))
gi = reggrid(h, 'delta', 'U', 'delta_tot', _gx, _gy)
im = ax[0].contourf(np.linspace(0, 0.4, 300), np.linspace(0, 5, 300), gi,
                    levels=100, cmap="jet", extend="both")
cb = fig.colorbar(im, ax=ax[0], pad=0.02)
cb.set_label(r"$\Delta_{tot}$", fontsize=16); cb.ax.tick_params(labelsize=11)
ax[0].plot(h.delta, h.U, "o", ms=5, mfc="none", mec="k", mew=1.0)
ax[0].set_xlabel(r"anisotropy  $\delta$"); ax[0].set_ylabel(r"$U/t$")
ax[0].set_title(r"(a) $L=14$, half filling", fontsize=17)

# (b) (n, delta) at L = 10, U = 4
g10 = d[(d.L == 10) & (d.U == 4.0)]
_gx2, _gy2 = np.meshgrid(np.linspace(0.5, 0.98, 300), np.linspace(0, 0.4, 300))
gi2 = reggrid(g10, 'n', 'delta', 'delta_tot', _gx2, _gy2)
im2 = ax[1].contourf(np.linspace(0.5, 0.98, 300), np.linspace(0, 0.4, 300), gi2,
                     levels=100, cmap="jet", extend="both")
cb2 = fig.colorbar(im2, ax=ax[1], pad=0.02)
cb2.set_label(r"$\Delta_{tot}$", fontsize=16); cb2.ax.tick_params(labelsize=11)
ax[1].plot(g10.n, g10.delta, "o", ms=5, mfc="none", mec="k", mew=1.0)
ax[1].set_xlabel(r"filling  $n$"); ax[1].set_ylabel(r"anisotropy  $\delta$")
ax[1].set_title(r"(b) $L=10$, $U=4$", fontsize=17)
fig.tight_layout()
fig.savefig(f"{D}/../figs/fig_dtot_heatmap.png", dpi=600, bbox_inches="tight",
            facecolor="white")

# ---------------- FIGURE 2: Delta n(k) over the Brillouin zone ----------------
# Pick L=14 half filling at increasing delta, matching the schematic's panels.
want = [(14, 1.0, 4.0, 0.0), (14, 1.0, 4.0, 0.2), (14, 1.0, 4.0, 0.4)]
sel = []
for L_, n_, U_, dl_ in want:
    hit = np.where((meta[:, 0] == L_) & np.isclose(meta[:, 1], n_) &
                   np.isclose(meta[:, 2], U_) & np.isclose(meta[:, 3], dl_))[0]
    if len(hit): sel.append((hit[0], dl_))

fig, ax = plt.subplots(1, len(sel), figsize=(5.0 * len(sel), 4.8))
vmax = max(np.abs(maps[i][:, 2]).max() for i, _ in sel)
for j, (i, dl_) in enumerate(sel):
    a = ax[j] if len(sel) > 1 else ax
    kx, ky, dn = maps[i][:, 0], maps[i][:, 1], maps[i][:, 2]
    # fold to [-pi, pi] so the zone is centred on Gamma
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
cb.set_label(r"$\Delta n(\mathbf{k}) = n_\uparrow - n_\downarrow$", fontsize=15)
fig.suptitle(r"$L=14$, half filling, $U=4$", fontsize=17, y=1.02)
fig.savefig(f"{D}/../figs/fig_dnk_maps.png", dpi=600, bbox_inches="tight",
            facecolor="white")

print("Delta_tot, L=14 half filling (rows U, cols delta):")
print(h.pivot_table(index="U", columns="delta", values="delta_tot").round(3).to_string())
print("\ncorrelation of dn(k) with sin(kx)sin(ky):")
print(h.pivot_table(index="U", columns="delta", values="corr_dxy").round(3).to_string())
