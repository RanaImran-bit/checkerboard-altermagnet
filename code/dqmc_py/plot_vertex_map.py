#!/usr/bin/env python3
"""Appeal heatmaps: d-wave VERTEX (connected, true-pairing) susceptibility and equal-time
correlation over the (t', tam) plane, all 4 reductions (maxk, k0, r0, r>2 avg), from the
validated finite-T DQMC at the appeal filling (6x6, U=4, beta=2, n~0.97, <sign>=1).
Top row: vertex SUSCEPTIBILITY suscV. Bottom row: vertex CORRELATION corrV (equal-time).
Result: tam ENHANCES the equal-time pairing vertex (corrV up); both knobs suppress the
dynamic suscV. A second figure shows the FULL (bubble-contaminated) quantities."""
import os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

root = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", ".."))
csv = sys.argv[1] if len(sys.argv) > 1 else "dqmc_vertex2d_6x6_b2.csv"
tag = csv.replace("dqmc_vertex2d_", "").replace(".csv", "")
a = np.array([l.split(",") for l in open(os.path.join(root, "results/dqmc_scan", csv)) if l[0].isdigit()], float)
tam = np.unique(a[:, 0]); t1 = np.unique(a[:, 1]); nt, n1 = len(tam), len(t1)
ext = [t1.min()-0.05, t1.max()+0.05, tam.min()-0.05, tam.max()+0.05]
reds = ["maxk", "k0", "r0", "rgt"]; rlab = ["peak-$q$", "$q{=}0$", "$R{=}0$", "$|R|{>}2$ avg"]
# column offsets: suscV 4, suscF 8, corrV 12, corrF 16
off = {"suscV": 4, "suscF": 8, "corrV": 12, "corrF": 16}

def make(fields, fname, suptitle):
    fig, ax = plt.subplots(len(fields), 4, figsize=(15, 3.6*len(fields)))
    if len(fields) == 1: ax = ax[None, :]
    for ri, (fld, flab, cm) in enumerate(fields):
        for ci, (rd, rl) in enumerate(zip(reds, rlab)):
            M = a[:, off[fld]+ci].reshape(nt, n1)
            im = ax[ri, ci].imshow(M, origin="lower", extent=ext, cmap=cm, aspect="auto",
                                   interpolation="bilinear")
            for i, tv in enumerate(tam):
                for j, t1v in enumerate(t1):
                    ax[ri, ci].text(t1v, tv, f"{M[i,j]:.2f}", ha="center", va="center",
                                    fontsize=6.5, color="w" if M[i,j] < M.mean() else "k")
            ax[ri, ci].set_title(f"{flab}: {rl}", fontsize=9.5)
            ax[ri, ci].set_xlabel(r"$t_1$ ($t'$)"); ax[ri, ci].set_ylabel(r"$t_{am}$")
            fig.colorbar(im, ax=ax[ri, ci], fraction=0.046, pad=0.04)
    fig.suptitle(suptitle, fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    out = os.path.join(root, "docs", fname); fig.savefig(out, dpi=120); print("wrote", out)

make([("suscV", r"vertex $\chi_d$ (dynamic)", "viridis"),
      ("corrV", r"vertex $C_d$ (equal-time)", "magma")],
     f"dqmc_vertex_map_{tag}.png",
     r"$d$-wave VERTEX (true pairing) vs ($t'$, $t_{am}$) -- DQMC 6$\times$6 $n{\approx}0.97$, $\beta{=}2$, $\langle$sign$\rangle{=}1$")
make([("suscF", r"FULL $\chi_d$", "viridis"), ("corrF", r"FULL $C_d$", "magma")],
     f"dqmc_full_map_{tag}.png",
     r"$d$-wave FULL (bubble-contaminated) vs ($t'$, $t_{am}$) -- same runs")
