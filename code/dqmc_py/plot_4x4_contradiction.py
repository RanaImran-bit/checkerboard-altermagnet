#!/usr/bin/env python3
"""Resolve the apparent 4x4 'susc UP / corr DOWN' contradiction by separating the bubble.
DQMC 4x4 U=4 n~0.88 beta=2 (sign=1), q=0 d-wave. LEFT (FULL = vertex+bubble): suscF rises
but corrF falls -- the apparent contradiction. RIGHT (VERTEX = true pairing): suscV AND corrV
both FALL -> consistent. The full susceptibility rises only because its non-interacting
BUBBLE rises with t1 (band reshaping), not pairing."""
import os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
root = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", ".."))
a = np.array([l.split(",") for l in open(os.path.join(root,"results/dqmc_scan","dqmc_vertex_4x4_b2.csv")) if l[0].isdigit()], float)
# cols: t1, dens, sign, suscF_maxk, suscV_maxk, corrF_maxk, corrV_maxk, suscF_k0, suscV_k0, corrF_k0, corrV_k0
t = a[:,0]; sF=a[:,7]; sV=a[:,8]; cF=a[:,9]; cV=a[:,10]; bub = sF - sV
fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.6))
ax[0].plot(t, sF/sF[0], "C0-o", lw=2.2, ms=8, label=r"susceptibility $\chi_d$ (full)")
ax[0].plot(t, cF/cF[0], "C3-s", lw=2.2, ms=8, label=r"correlation $C_d$ (full, $\tau{=}0$)")
ax[0].axhline(1, color="0.7", lw=.8, ls=":")
ax[0].set_title(r"FULL (= pairing + bubble): the apparent contradiction")
ax[0].annotate("susc UP", (0.28, 1.045), color="C0", fontsize=10)
ax[0].annotate("corr DOWN", (0.25, 0.945), color="C3", fontsize=10)
ax[1].plot(t, sV/sV[0], "C0-o", lw=2.2, ms=8, label=r"vertex susceptibility $\chi_d^{V}$")
ax[1].plot(t, cV/cV[0], "C3-s", lw=2.2, ms=8, label=r"vertex correlation $C_d^{V}$ ($\tau{=}0$)")
ax[1].plot(t, bub/bub[0], "C2--^", lw=1.6, ms=7, mfc="none", label=r"bubble (= full $-$ vertex)")
ax[1].axhline(1, color="0.7", lw=.8, ls=":")
ax[1].set_title(r"VERTEX (true pairing): both DOWN $\Rightarrow$ consistent")
ax[1].annotate("bubble UP\n(drives full susc)", (0.21, 1.04), color="C2", fontsize=9)
ax[1].annotate("pairing DOWN", (0.05, 0.80), color="C0", fontsize=10)
for x in ax:
    x.set_xlabel(r"$t_1$ (anisotropic NNN $t'$)"); x.set_ylabel("normalized to $t_1{=}0$")
    x.grid(alpha=.3); x.legend(frameon=False, fontsize=9)
fig.suptitle(r"4$\times$4 ($U{=}4$, $n{\approx}0.88$, $\beta{=}2$): the susc-up/corr-down 'contradiction' is the BUBBLE", fontsize=11)
fig.tight_layout(rect=[0,0,1,0.95])
out = os.path.join(root,"docs","dqmc_4x4_contradiction.png"); fig.savefig(out, dpi=130); print("wrote", out)
