#!/usr/bin/env python3
"""(n, tam) phase-diagram heatmap of the local d-wave vertex N_pp(R=0) = C_d(0),
the pyqmc analog of the manuscript's central figure. Reads the pair_rspace --csv dump:
lx ly nup ndn U tam t1 seed Npp_R0 e Npp_Rgt e qsum e."""
import sys, argparse, collections
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.interpolate import griddata


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("infile")
    ap.add_argument("-o", "--out", default="docs/pd_heatmap_Cd0.png"); a = ap.parse_args()
    pts = {}  # (n,tam) -> Npp_R0
    for l in open(a.infile):
        p = l.split()
        if len(p) >= 14 and p[0].isdigit() and p[1].isdigit():
            lx, ly = int(p[0]), int(p[1]); N = int(p[2]); tam = round(float(p[5]), 3)
            n = 2 * N / (lx * ly); pts[(round(n, 4), tam)] = float(p[8])
    ns = sorted({k[0] for k in pts}); tams = sorted({k[1] for k in pts})
    Z = np.full((len(tams), len(ns)), np.nan)
    for i, t in enumerate(tams):
        for j, n in enumerate(ns):
            if (n, t) in pts: Z[i, j] = pts[(n, t)]
    fig, ax = plt.subplots(figsize=(7.6, 5.2))
    # smooth interpolated background + raw cell annotations
    NN, TT = np.meshgrid(ns, tams)
    mask = ~np.isnan(Z)
    gx, gy = np.meshgrid(np.linspace(min(ns), max(ns), 200), np.linspace(min(tams), max(tams), 200))
    Zi = griddata((NN[mask], TT[mask]), Z[mask], (gx, gy), method="cubic")
    im = ax.imshow(Zi, origin="lower", aspect="auto", cmap="inferno",
                   extent=[min(ns), max(ns), min(tams), max(tams)])
    ax.scatter(NN[mask], TT[mask], c=Z[mask], cmap="inferno", edgecolors="white",
               s=60, linewidths=0.6, vmin=np.nanmin(Zi), vmax=np.nanmax(Zi))
    for i, t in enumerate(tams):
        for j, n in enumerate(ns):
            if mask[i, j]:
                ax.text(n, t, f"{Z[i, j]:.2f}", ha="center", va="center", fontsize=6.5,
                        color="white")
    ax.set_xlabel(r"density $n$", fontsize=13)
    ax.set_ylabel(r"altermagnet anisotropy $t_A$", fontsize=13)
    ax.set_title(r"pyqmc d-wave vertex $N_{pp}(R{=}0)\equiv C_d(0)$   (8×8, $U=4$)", fontsize=12)
    fig.colorbar(im, ax=ax, label=r"$C_d(0)$")
    fig.tight_layout(); fig.savefig(a.out, dpi=150); print("wrote", a.out)


if __name__ == "__main__":
    main()
