#!/usr/bin/env python3
"""Does ISOTROPIC t' enhance d-wave? Compare the spin-INDEPENDENT (isotropic, conventional)
NNN t' = tp against the spin-DEPENDENT (anisotropic, altermagnet) t' = t1, for the d-wave
pair susceptibility chi_d^maxk on the SAME 4x4 U=4 doped lattice. Two methods: T=0 CPQMC and
finite-T DQMC (beta=2). Conclusion: only the ANISOTROPIC t' gives a strong, monotonic
enhancement (it carries a d-wave form factor); the isotropic t' gives at most a weak,
non-monotonic bump (it couples to d-wave only indirectly, via the band/Fermi surface)."""
import os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

root = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", ".."))
D = os.path.join(root, "results/dqmc_scan")
def dq(name): 
    a = np.array([l.split(",") for l in open(os.path.join(D, name)) if l[0].isdigit()], float)
    return a[:, 0], a[:, 3]                                  # knob, sc_maxk
def cp(name, col):
    a = np.array([l.split() for l in open(os.path.join(D, name)) if l[0].isdigit()], float)
    return a[:, 0], a[:, col]

fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
# CPQMC T=0
k, t1 = cp("cpqmc_t1_4x4_n0875.txt", 2); _, tp = cp("cpqmc_tp_4x4_n0875.txt", 2)
ax[0].plot(k, t1, "C2-o", lw=2.2, ms=8, label=r"anisotropic $t'$ ($t_1$, spin-dep.)")
ax[0].plot(k, tp, "C1--s", lw=2.2, ms=8, mfc="none", label=r"isotropic $t'$ ($t_p$, spin-indep.)")
ax[0].set_title(r"T=0 CPQMC")
# DQMC beta=2
k1, x1 = dq("dqmc_t1_4x4_b2.csv"); kp, xp = dq("dqmc_tp_4x4_b2.csv")
ax[1].plot(k1, x1, "C2-o", lw=2.2, ms=8, label=r"anisotropic $t'$ ($t_1$)")
ax[1].plot(kp, xp, "C1--s", lw=2.2, ms=8, mfc="none", label=r"isotropic $t'$ ($t_p$)")
ax[1].set_title(r"finite-$T$ DQMC ($\beta{=}2$)")
for a in ax:
    a.set_xlabel(r"$t'$ strength"); a.set_ylabel(r"$d$-wave susceptibility $\chi_d^{\max k}$")
    a.grid(alpha=.3)
ax[0].legend(frameon=False, loc="center right"); ax[1].legend(frameon=False, loc="lower left")
ax[0].annotate("anisotropic:\nstrong ENHANCE", (0.20, 9.5), fontsize=9, color="C2")
ax[0].annotate("isotropic: weak\nbump, turns over", (0.03, 7.4), fontsize=9, color="C1")
fig.suptitle(r"Only ANISOTROPIC $t'$ enhances $d$-wave; isotropic $t'$ does not "
             r"(4$\times$4, $U{=}4$, doped)", fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.95])
out = os.path.join(root, "docs", "iso_vs_aniso_tprime.png")
fig.savefig(out, dpi=130); print("wrote", out)
