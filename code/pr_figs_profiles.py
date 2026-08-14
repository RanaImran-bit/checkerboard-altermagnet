"""Rebuilt Fig 1 and Fig 3.

Fig 1 was on a log axis out to R = 8.5, where the vertex is pure noise and the
panel became a tangle of crossing lines. The signal lives at R <= 4.5, and the
SIGN is the whole point, so it is now linear in v(R) over that range. R = 0 is
left out because it is 20x larger than everything else and already carries its
own panel in Fig 2.

Fig 3 dropped on-site s: it is large and negative, it is not a competitor for
long-range pairing, and at U = 2 it stretched the axis so the three real channels
were squashed into a band a few pixels tall.
"""
import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt
exec(open("figs.py").read().split("# ============================ FIGURE 1")[0])

CH3 = [c for c in CH if c[0] != "son"]
DL = [0.0, 0.2, 0.4, 0.7]
RMAX = 4.55

# ------------------------- FIGURE 1 rebuilt -------------------------------
fig, ax = plt.subplots(1, 4, figsize=(15.5, 4.2), sharey=True)
for j, dl in enumerate(DL):
    a = ax[j]; t = prof(4.0, dl)
    a.axhline(0, color="k", lw=1.0, ls="-", zorder=1)
    a.axvspan(2.0, RMAX, color="#f2f2f2", zorder=0)
    for c, lab, col, mk in CH3:
        s = t[(t.channel == c) & (t.R >= 1) & (t.R <= RMAX)].sort_values("R")
        a.errorbar(s.R, s["mean"] * 1e3, yerr=s["sem"] * 1e3, color=col, marker=mk,
                   ms=6.5, lw=1.7, capsize=2.5, label=lab, zorder=3)
    a.set_xlabel(r"separation  $R$")
    a.set_title(rf"$\delta = {dl:g}$", fontsize=13.5)
    a.set_xlim(0.7, RMAX); a.set_ylim(-1.6, 5.7)
    a.set_xticks([1, np.sqrt(2), 2, 3, 4])
    a.set_xticklabels(["1", r"$\sqrt{2}$", "2", "3", "4"])
for a, txt in zip(ax, ["bond and diagonal\nboth modest",
                       r"$d_{xy}$ diagonal turns on",
                       r"$d_{xy}$ dominates",
                       r"$d_{xy}$ weakening"]):
    a.text(0.97, 0.96, txt, transform=a.transAxes, ha="right", va="top",
           fontsize=10.5, color="0.25", style="italic")
ax[0].set_ylabel(r"per-pair vertex  $v(R)\times 10^{3}$")
ax[0].legend(fontsize=11, loc="center right")
ax[3].text(3.2, 4.6, "shaded = long range", fontsize=9.5, color="0.45")
fig.suptitle(r"Real-space pairing vertex at $U=4$, half filling, $L=12$"
             "\n" r"$R=1$ is the nearest-neighbour bond, $R=\sqrt{2}$ the plaquette "
             r"diagonal that $\delta$ acts on", fontsize=13.5, y=1.08)
fig.tight_layout(); fig.savefig("fig_pr_profiles.png", dpi=300, bbox_inches="tight",
                                facecolor="white")

# ------------------------- FIGURE 3 rebuilt -------------------------------
DLS = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
US = [2.0, 4.0, 6.0, 8.0]
fig, ax = plt.subplots(1, 4, figsize=(15.5, 4.0), sharey=True)
for j, U in enumerate(US):
    a = ax[j]; t = {dl: local_long(U, dl) for dl in DLS}
    a.axhline(0, color="k", lw=1.0, ls="--")
    for c, lab, col, mk in CH3:
        y = [t[dl][c][2] * 1e3 for dl in DLS]
        e = [t[dl][c][3] * 1e3 for dl in DLS]
        a.errorbar(DLS, y, yerr=e, color=col, marker=mk, ms=6.5, lw=1.7,
                   capsize=3, label=lab)
    a.set_xlabel(r"anisotropy  $\delta$")
    a.set_title(rf"$U = {U:g}$", fontsize=13.5)
ax[0].set_ylabel(r"long-range vertex  $\langle v\rangle_{R>2}\times 10^{3}$")
ax[0].set_ylim(-0.045, 0.185)
ax[0].legend(fontsize=11, loc="upper right")
fig.suptitle(r"Long-range pairing vertex against anisotropy, half filling, $L=12$"
             "\n" r"on-site $s$ omitted: strongly repulsive and not a competitor here",
             fontsize=13.5, y=1.07)
fig.tight_layout(); fig.savefig("fig_pr_longrange_U.png", dpi=300, bbox_inches="tight",
                                facecolor="white")
print("rebuilt")
