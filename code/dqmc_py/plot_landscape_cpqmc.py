#!/usr/bin/env python3
"""T=0 CPQMC (t', tam) d-wave pairing landscape (4x4 U=4, nup=ndn=7, n=0.875), the
ground-state companion to docs/dqmc_landscape_heatmap.png. LEFT: SUSCEPTIBILITY suscF^maxk
-- enhanced by t', suppressed by tam, with a sharp t'-vs-tam COMPETITION boundary (need
enough t' to overcome tam). RIGHT: equal-time CORRELATION corrF^maxk -- the mirror. Same
trends as the finite-T DQMC landscape but with a sharper T=0 structure: temperature-robust."""
import os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

root = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", ".."))
rows = [l.split() for l in open(os.path.join(root, "results/dqmc_scan",
        "cpqmc_tam_t1_landscape_4x4.txt")) if l[0].isdigit()]
a = np.array(rows, float)
tam = np.unique(a[:, 0]); t1 = np.unique(a[:, 1])
susc = a[:, 3].reshape(len(tam), len(t1)); corr = a[:, 4].reshape(len(tam), len(t1))
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
fig.suptitle(r"T=0 CPQMC $d$-wave pairing landscape vs ($t'$, $t_{am}$) "
             r"(4$\times$4, $U{=}4$, $n{=}0.875$) -- $t'$ vs $t_{am}$ competition", fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.95])
out = os.path.join(root, "docs", "cpqmc_landscape_heatmap.png")
fig.savefig(out, dpi=130); print("wrote", out)
