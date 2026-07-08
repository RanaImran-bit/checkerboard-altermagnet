#!/usr/bin/env python3
"""Definitive three-method cross-validation capstone, vs anisotropic NNN t' on the SAME
4x4 U=4 doped lattice (n~0.88), all <sign>=1. THREE methods:
  - T=0 CPQMC (constrained path),
  - finite-T DQMC (beta=2, exact, sign-problem-free here),
  - finite-T CONSTRAINED-PATH AFQMC (FT-CPMC, beta=2, sign-controlled).
LEFT: tau-integrated d-wave SUSCEPTIBILITY rises with t'.  RIGHT: equal-time d-wave
CORRELATION falls with t'.  The two finite-T curves (DQMC exact, FT-CPMC constrained)
overlap across the whole scan -- the new sign-controlled method reproduces the exact one --
and the T=0 CPQMC shows the same (steeper) trends.  This is the susceptibility-vs-correlation
distinction that resolves the manuscript discrepancy, now confirmed by three methods."""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

root = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", ".."))
D = os.path.join(root, "results/dqmc_scan")


def arr(name, sep, cols):
    rows = [l.replace(",", " ").split() for l in open(os.path.join(D, name)) if l[0].isdigit()]
    a = np.array(rows, float)
    return {k: a[:, c] for k, c in cols.items()}


# correlation (equal-time, q=0) and susceptibility (tau-integrated, peak-q) for each method
ft_c = arr("ftcpmc_t1_4x4_b2.txt", " ", {"t": 0, "cd": 3})                 # S_d
ft_x = arr("ftcpmc_susc_t1_4x4_b2.txt", " ", {"t": 0, "xd": 3})            # chi_d
dq = arr("dqmc_t1_4x4_b2.csv", ",", {"t": 0, "xd": 3, "cd": 6})           # sc_maxk, cr_maxk
cp = arr("cpqmc_t1_4x4_n0875.txt", " ", {"t": 0, "xd": 2, "cd": 5})       # suscF_maxk, corrF_maxk

fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.6))
# --- SUSCEPTIBILITY (normalized) ---
ax[0].plot(cp["t"], cp["xd"] / cp["xd"][0], "C0-o", lw=2, ms=8, label="T=0 CPQMC")
ax[0].plot(dq["t"], dq["xd"] / dq["xd"][0], "C3--s", lw=2, ms=9, mfc="none", label=r"DQMC $\beta{=}2$ (exact)")
ax[0].plot(ft_x["t"], ft_x["xd"] / ft_x["xd"][0], "C2:^", lw=2, ms=8, label=r"FT-CPMC $\beta{=}2$ (constr.)")
ax[0].axhline(1, color="0.7", lw=.8, ls=":")
ax[0].set_title(r"$d$-wave SUSCEPTIBILITY $\chi_d(t')$  $\Rightarrow$ UP")
ax[0].set_ylabel(r"normalized $\chi_d$")
# --- CORRELATION (normalized) ---
ax[1].plot(cp["t"], cp["cd"] / cp["cd"][0], "C0-o", lw=2, ms=8, label="T=0 CPQMC")
ax[1].plot(dq["t"], dq["cd"] / dq["cd"][0], "C3--s", lw=2, ms=9, mfc="none", label=r"DQMC $\beta{=}2$ (exact)")
ax[1].plot(ft_c["t"], ft_c["cd"] / ft_c["cd"][0], "C2:^", lw=2, ms=8, label=r"FT-CPMC $\beta{=}2$ (constr.)")
ax[1].axhline(1, color="0.7", lw=.8, ls=":")
ax[1].set_title(r"equal-time CORRELATION $C_d(t')$  $\Rightarrow$ DOWN")
ax[1].set_ylabel(r"normalized $C_d$")
ax[0].annotate("DQMC & FT-CPMC\noverlap", (0.13, 1.02), fontsize=8.5, color="0.3")
for a in ax:
    a.set_xlabel(r"$t_1$  (anisotropic NNN $t'$)"); a.grid(alpha=.3); a.legend(frameon=False, fontsize=9)
fig.suptitle(r"Three-method cross-validation vs $t'$: SUSCEPTIBILITY up, CORRELATION down "
             r"(4$\times$4, $U{=}4$, $n{\approx}0.88$, all $\langle$sign$\rangle{=}1$)", fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.95])
out = os.path.join(root, "docs", "xval_three_method_t1.png")
fig.savefig(out, dpi=130); print("wrote", out)
