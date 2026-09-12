"""Doped fillings for P(R), and the delta-polarisation plot at every filling."""
import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt
exec(open("figs.py").read().split("# ============================ FIGURE 1")[0])

CH3 = [c for c in CH if c[0] != "son"]
DLS = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
NUPS = [36, 40, 48, 56, 64, 72]
NLAB = {k: 2 * k / 144 for k in NUPS}

# ---------------- FIGURE 4: long-range vertex vs delta, per filling -------
fig, ax = plt.subplots(2, 3, figsize=(14.2, 7.6), sharex=True)
for j, nu in enumerate(NUPS):
    a = ax.flat[j]
    t = {dl: local_long(4.0, dl, nu) for dl in DLS}
    a.axhline(0, color="k", lw=1.0, ls="--")
    for c, lab, col, mk in CH3:
        y = [t[dl][c][2] * 1e3 for dl in DLS]
        e = [t[dl][c][3] * 1e3 for dl in DLS]
        a.errorbar(DLS, y, yerr=e, color=col, marker=mk, ms=6, lw=1.6,
                   capsize=3, label=lab)
    a.set_title(rf"$n = {NLAB[nu]:.3f}$", fontsize=13)
    if j >= 3: a.set_xlabel(r"anisotropy  $\delta$")
    if j % 3 == 0: a.set_ylabel(r"$\langle v\rangle_{R>2}\times 10^{3}$")
ax.flat[0].legend(fontsize=10.5, loc="best")
fig.suptitle(r"Long-range pairing vertex against anisotropy at every filling "
             r"($U=4$, $L=12$)", fontsize=14, y=1.0)
fig.tight_layout(); fig.savefig("fig_pr_fillings.png", dpi=300, bbox_inches="tight",
                                facecolor="white")

# ---------------- FIGURE 5: local and long range vs filling at fixed delta -
fig, ax = plt.subplots(1, 2, figsize=(11.6, 4.3))
for k, (idx, ttl) in enumerate([(0, r"local  $v(R=0)$"),
                                (2, r"long range  $\langle v\rangle_{R>2}$")]):
    a = ax[k]; a.axhline(0, color="k", lw=1.0, ls="--")
    for c, lab, col, mk in CH3:
        y, e = [], []
        for nu in NUPS:
            t = local_long(4.0, 0.4, nu)[c]
            y.append(t[idx] * 1e3); e.append(t[idx + 1] * 1e3)
        a.errorbar([NLAB[k2] for k2 in NUPS], y, yerr=e, color=col, marker=mk,
                   ms=6.5, lw=1.7, capsize=3, label=lab)
    a.set_xlabel(r"filling  $n$"); a.set_title(ttl, fontsize=13)
    a.set_ylabel(r"$v\times 10^{3}$")
ax[0].legend(fontsize=10.5, loc="best")
fig.suptitle(r"Where in filling each channel pairs, at $\delta=0.4$, $U=4$",
             fontsize=14, y=1.02)
fig.tight_layout(); fig.savefig("fig_pr_vs_filling.png", dpi=300, bbox_inches="tight",
                                facecolor="white")

print("delta=0.4 U=4  long-range vertex x1000 vs filling")
for nu in NUPS:
    t = local_long(4.0, 0.4, nu)
    print(f"  n={NLAB[nu]:.3f}: sext {t['sext'][2]*1e3:7.3f}  d {t['d'][2]*1e3:7.3f}"
          f"  dxy {t['dxy'][2]*1e3:7.3f}   winner "
          f"{max(['sext','d','dxy'], key=lambda c: t[c][2])}")
print()
print("long-range winner map (rows n, cols delta), U=4")
for nu in NUPS:
    w = [max(["sext", "d", "dxy"], key=lambda c: local_long(4.0, dl, nu)[c][2])
         for dl in DLS]
    print(f"  n={NLAB[nu]:.3f}: " + " ".join(f"{x:>4}" for x in w))
