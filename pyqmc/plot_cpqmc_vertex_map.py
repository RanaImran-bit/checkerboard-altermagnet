#!/usr/bin/env python3
"""CPQMC (T=0) counterpart of code/dqmc_py/plot_vertex_map.py: d-wave VERTEX (connected,
true-pairing) and FULL susceptibility + equal-time correlation over the (t', tam) plane,
all 4 reductions (maxk, k0, R=0, |R|>2 avg), from the T=0 constrained-path CPQMC at 6x6,
U=4, n=0.944 (nup=ndn=17, the nearest closed canonical filling to the appeal n~0.97 --
the bare appeal filling is open-shell for the free-particle trial). Directly comparable to
the finite-T DQMC maps (dqmc_vertex_map_6x6_b2.png / dqmc_full_map_6x6_b2.png)."""
import os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

root = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", ""))
csv = sys.argv[1] if len(sys.argv) > 1 else "cpqmc_vertex2d_6x6.csv"
tag = csv.replace("cpqmc_vertex2d_", "").replace(".csv", "")
a = np.array([l.split(",") for l in open(os.path.join(root, "results/dqmc_scan", csv)) if l[0].isdigit()], float)
dens = a[0, 2]; latt = tag.split("_")[0]
tam = np.unique(a[:, 0]); t1 = np.unique(a[:, 1]); nt, n1 = len(tam), len(t1)
ext = [t1.min()-0.05, t1.max()+0.05, tam.min()-0.05, tam.max()+0.05]
reds = ["maxk", "k0", "r0", "rgt"]; rlab = ["peak-$q$", "$q{=}0$", "$R{=}0$", "$|R|{>}2$ avg"]
off = {"suscV": 4, "suscF": 8, "corrV": 12, "corrF": 16}


def make(fields, fname, suptitle):
    plt.rcParams.update({"font.size": 13})
    fig, ax = plt.subplots(len(fields), 4, figsize=(22, 5.6*len(fields)))
    if len(fields) == 1:
        ax = ax[None, :]
    for ri, (fld, flab, cm) in enumerate(fields):
        for ci, (rd, rl) in enumerate(zip(reds, rlab)):
            M = a[:, off[fld]+ci].reshape(nt, n1)
            im = ax[ri, ci].imshow(M, origin="lower", extent=ext, cmap=cm, aspect="auto",
                                   interpolation="bilinear")
            ax[ri, ci].contour(t1, tam, M, 5, colors="w", linewidths=0.5, alpha=0.4)
            fmt = "%.0f" if np.nanmax(np.abs(M)) >= 5 else "%.2f"
            for i, tv in enumerate(tam):
                for j, t1v in enumerate(t1):
                    ax[ri, ci].text(t1v, tv, fmt % M[i, j], ha="center", va="center",
                                    fontsize=13, fontweight="bold",
                                    color="w" if M[i, j] < 0.55*np.nanmax(M) else "k")
            ax[ri, ci].set_title(f"{flab}: {rl}", fontsize=15)
            ax[ri, ci].set_xlabel(r"$t_1$ ($t'$)", fontsize=13); ax[ri, ci].set_ylabel(r"$t_{am}$", fontsize=13)
            fig.colorbar(im, ax=ax[ri, ci], fraction=0.046, pad=0.04)
    fig.suptitle(suptitle, fontsize=17)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out = os.path.join(root, "docs", fname); fig.savefig(out, dpi=140); print("wrote", out)


make([("suscV", r"vertex $\chi_d$ (dynamic)", "viridis"),
      ("corrV", r"vertex $C_d$ (equal-time)", "magma")],
     f"cpqmc_vertex_map_{tag}.png",
     rf"$d$-wave VERTEX (true pairing) vs ($t'$, $t_{{am}}$) -- T=0 CPQMC {latt} $n{{=}}{dens:.3f}$")
make([("suscF", r"FULL $\chi_d$", "viridis"), ("corrF", r"FULL $C_d$", "magma")],
     f"cpqmc_full_map_{tag}.png",
     rf"$d$-wave FULL (bubble-contaminated) vs ($t'$, $t_{{am}}$) -- T=0 CPQMC {latt}, same runs")
