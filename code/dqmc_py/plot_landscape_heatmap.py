#!/usr/bin/env python3
"""(t', tam) parameter-landscape heatmaps of the finite-T DQMC d-wave pairing on the 4x4
U=4 doped lattice (n~0.88, beta=2, all <sign>=1). LEFT: tau-integrated SUSCEPTIBILITY
chi_d^maxk -- bright at high t', low tam. RIGHT: equal-time CORRELATION C_d^maxk -- bright
at high tam, low t'. The two are anti-correlated: t' enhances the susceptibility, tam
enhances the correlation. One picture of the whole cross-validation conclusion."""
import os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

root = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", ".."))
rows = [l.split() for l in open(os.path.join(root, "results/dqmc_scan",
        "dqmc_tam_t1_landscape_4x4_b2.txt")) if l[0].isdigit()]
a = np.array(rows, float)
tam = np.unique(a[:, 0]); t1 = np.unique(a[:, 1])
susc = a[:, 4].reshape(len(tam), len(t1))      # rows=tam, cols=t1
corr = a[:, 5].reshape(len(tam), len(t1))
ext = [t1.min() - 0.05, t1.max() + 0.05, tam.min() - 0.05, tam.max() + 0.05]

fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.8))
for a_, M, ttl, cm in ((ax[0], susc, r"$d$-wave SUSCEPTIBILITY  $\chi_d^{\max k}$", "viridis"),
                       (ax[1], corr, r"equal-time CORRELATION  $C_d^{\max k}$", "magma")):
    im = a_.imshow(M, origin="lower", extent=ext, cmap=cm, aspect="auto", interpolation="bilinear")
    for i, tv in enumerate(tam):
        for j, t1v in enumerate(t1):
            a_.text(t1v, tv, f"{M[i, j]:.1f}", ha="center", va="center", fontsize=7.5,
                    color="w" if M[i, j] < M.mean() else "k")
    a_.set_xlabel(r"$t_1$  (anisotropic NNN $t'$)"); a_.set_ylabel(r"$t_{am}$  (NN altermagnet)")
    a_.set_title(ttl); fig.colorbar(im, ax=a_, fraction=0.046, pad=0.04)
ax[0].annotate("", xy=(0.4, 0.0), xytext=(0.0, 0.0), arrowprops=dict(arrowstyle="->", color="w", lw=1.5))
ax[1].annotate("", xy=(0.0, 0.4), xytext=(0.0, 0.0), arrowprops=dict(arrowstyle="->", color="w", lw=1.5))
fig.suptitle(r"Finite-$T$ DQMC $d$-wave pairing landscape vs ($t'$, $t_{am}$) "
             r"(4$\times$4, $U{=}4$, $\beta{=}2$, $n{\approx}0.88$, $\langle$sign$\rangle{=}1$)", fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.95])
out = os.path.join(root, "docs", "dqmc_landscape_heatmap.png")
fig.savefig(out, dpi=130); print("wrote", out)
