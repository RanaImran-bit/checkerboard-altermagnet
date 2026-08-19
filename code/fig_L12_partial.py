"""L=12 Fortran grid, half filling -- what has finished so far.

28 of 49 (U, delta) cells are complete. A run counts as complete only when
nohup.out ends with the timing line AND both n_up.dat and n_dn.dat exist;
partial runs write valid-looking files with unconverged statistics, so they are
excluded rather than averaged in.

Two things this set adds over the L=14 grid, which stopped at delta = 0.4:

  1. delta reaches 0.6, where the diagonal hopping t+ = t' + delta = 0.9 is
     comparable to t. This is the regime where the dxy channel was expected to
     take over, and it does.

  2. Delta_tot has a sharp onset in U that MOVES UP with delta. Panel (c) shows
     why that onset must be read with care: it tracks the staggered moment of
     the UHF trial, not something the QMC found on its own. Below the trial's
     own Stoner threshold the trial is very nearly spin-degenerate, and the
     constrained path cannot manufacture a splitting the trial does not have.
     Above it, the QMC value is NOT proportional to the trial moment -- along
     delta = 0.1 the trial moment rises 0.15 -> 0.68 while Delta_tot FALLS
     0.061 -> 0.028 -- so the trend above onset is QMC, the onset itself is not.

Delta_tot = sum_k |n_up(k) - n_dn(k)| / L^2   (PRL Eq. 2, absolute value inside
the sum, so it is finite even though the net magnetisation vanishes).
Pairing values are the vertex contributions at k = (0,0), full - L*Unpair, the
Fortran convention.
"""
import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt

mpl.rcParams.update({
    "font.family": "serif", "font.serif": ["DejaVu Serif"],
    "mathtext.fontset": "dejavuserif", "axes.linewidth": 1.2,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True, "legend.frameon": False,
    "axes.labelsize": 15, "xtick.labelsize": 12, "ytick.labelsize": 12,
    "figure.dpi": 130})

D = "/Users/liujiaxin/Desktop/checkerboard-altermagnet"
d = pd.read_csv(f"{D}/data/fortran_L12_partial.csv")
d = d[np.isclose(d.n, 1.0)].sort_values(["U", "delta"])
US, DS = sorted(d.U.unique()), sorted(d.delta.unique())
CU = plt.cm.viridis(np.linspace(0.05, 0.88, len(US)))
CD = plt.cm.plasma(np.linspace(0.05, 0.82, len(DS)))
MK = ["o", "s", "^", "D", "v", "P", "X"]

# ============================ FIGURE 1: Delta_tot ============================
fig, ax = plt.subplots(1, 3, figsize=(16.2, 4.9))

for i, U in enumerate(US):                                   # (a) against delta
    s = d[d.U == U]
    ax[0].plot(s.delta, s.dtot_N, marker=MK[i % len(MK)], ms=6.5, lw=1.6,
               color=CU[i], label=f"$U={U:g}$")
ax[0].set_xlabel(r"anisotropy  $\delta$")
ax[0].set_ylabel(r"$\Delta_{\mathrm{tot}}$")
ax[0].legend(fontsize=10, ncol=2, loc="best")
ax[0].set_title("(a)", loc="left", fontsize=15)

for i, dl in enumerate(DS):                                      # (b) against U
    s = d[d.delta == dl]
    if len(s) < 2: continue
    ax[1].plot(s.U, s.dtot_N, marker=MK[i % len(MK)], ms=6.5, lw=1.6,
               color=CD[i], label=rf"$\delta={dl:g}$")
ax[1].set_xlabel(r"$U/t$"); ax[1].set_ylabel(r"$\Delta_{\mathrm{tot}}$")
ax[1].legend(fontsize=10, ncol=2, loc="best")
ax[1].set_title("(b)", loc="left", fontsize=15)

# (c) the trial-dependence diagnostic. Colour by delta so the onset structure
# is visible: the small-moment cluster near the origin is entirely delta >= 0.4.
sc = ax[2].scatter(d.m_trial, d.dtot_N, c=d.delta, cmap="plasma", s=52,
                   edgecolor="k", linewidth=0.6, zorder=3)
cb = fig.colorbar(sc, ax=ax[2], pad=0.02); cb.set_label(r"$\delta$", fontsize=14)
ax[2].axvspan(0, 0.12, color="0.85", zorder=0)
ax[2].text(0.055, 0.9 * d.dtot_N.max(), "trial below\nits own\nStoner point",
           ha="center", va="top", fontsize=9.5, color="0.35")
ax[2].set_xlabel(r"UHF trial staggered moment  $m_{\mathrm{trial}}$")
ax[2].set_ylabel(r"$\Delta_{\mathrm{tot}}$  (QMC)")
ax[2].set_title("(c)", loc="left", fontsize=15)

fig.tight_layout()
fig.savefig(f"{D}/figs/fig_L12_dtot_partial.png", dpi=600, bbox_inches="tight",
            facecolor="white")

# ======================= FIGURE 2: pairing, dxy vs dx2-y2 ====================
# Two colours only, as agreed: blue for dx2-y2, red for dxy. Extended s is drawn
# in grey for scale -- it is the largest channel throughout and would otherwise
# make the d-channel comparison look larger than it is.
fig2, ax2 = plt.subplots(1, 2, figsize=(11.6, 4.9))

for i, U in enumerate([2.0, 3.0, 3.5, 4.0]):                # (a) channels vs delta
    s = d[(d.U == U) & (d.U > 0)].sort_values("delta")
    if len(s) < 2: continue
    al = 0.35 + 0.65 * i / 3
    ax2[0].plot(s.delta, s.d,   "-o", ms=5.5, lw=1.6, color="tab:blue",
                alpha=al, label=rf"$d_{{x^2-y^2}}$, $U={U:g}$")
    ax2[0].plot(s.delta, s.dxy, "-s", ms=5.5, lw=1.6, color="tab:red",
                alpha=al, label=rf"$d_{{xy}}$, $U={U:g}$")
ax2[0].axhline(0, color="k", lw=0.8, ls=":")
ax2[0].set_xlabel(r"anisotropy  $\delta$")
ax2[0].set_ylabel("pairing vertex at $\\mathbf{k}=0$")
ax2[0].legend(fontsize=8.5, ncol=2, loc="best")
ax2[0].set_title("(a)", loc="left", fontsize=15)

# (b) the crossover itself: which d channel wins, cell by cell
w = d[d.U > 0].copy()
w["diff"] = w.dxy - w.d
v = np.abs(w["diff"]).max()
sc2 = ax2[1].scatter(w.delta, w.U, c=w["diff"], cmap="coolwarm",
                     vmin=-v, vmax=v, s=280, marker="s",
                     edgecolor="k", linewidth=0.7)
cb2 = fig2.colorbar(sc2, ax=ax2[1], pad=0.02)
cb2.set_label(r"$d_{xy} - d_{x^2-y^2}$", fontsize=13)
ax2[1].set_xlabel(r"anisotropy  $\delta$"); ax2[1].set_ylabel(r"$U/t$")
ax2[1].set_xlim(0.03, 0.67); ax2[1].set_ylim(1.5, 5.5)
ax2[1].set_title("(b)  red = $d_{xy}$ wins", loc="left", fontsize=15)

fig2.tight_layout()
fig2.savefig(f"{D}/figs/fig_L12_pairing_partial.png", dpi=600,
             bbox_inches="tight", facecolor="white")

# ------------------------------- what it says -------------------------------
print(f"complete cells: {len(d)} of 49\n")
print("Delta_tot (rows U, cols delta):")
print(d.pivot_table(index="U", columns="delta", values="dtot_N").round(4).to_string())
print("\nd_xy - d_x2y2  (positive = dxy wins):")
print(d[d.U > 0].assign(dif=lambda f: f.dxy - f.d)
       .pivot_table(index="U", columns="delta", values="dif").round(3).to_string())
