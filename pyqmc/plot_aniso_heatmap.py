#!/usr/bin/env python3
"""Plot the (tam, t1) d-wave pairing-susceptibility connected-vertex heatmap from the
wide cluster scan CSVs. Diverging colormap centred at 0; the zero contour is the
d-wave pairing boundary. Reads the collector dump (===CSV6x6=== / ===CSV8x8=== blocks)
or the raw chid_tam_scan CSV lines, aggregates chi_d_vtx over seeds.

    python pyqmc/plot_aniso_heatmap.py <collector_output.txt> -o docs/anisotropy_heatmap.png
"""
import sys, argparse, collections
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load(path, COL=10):
    agg = collections.defaultdict(list)   # (lx,tam,t1) -> [chi_d_vtx,...]
    for ln in open(path):
        p = ln.split()
        if len(p) >= 18 and p[0].isdigit() and p[1].isdigit():
            lx = int(p[0]); tam = round(float(p[5]), 3); t1 = round(float(p[6]), 3)
            agg[(lx, tam, t1)].append(float(p[COL]))
    return agg


def grid(agg, lx):
    tams = sorted({k[1] for k in agg if k[0] == lx})
    t1s = sorted({k[2] for k in agg if k[0] == lx})
    Z = np.full((len(tams), len(t1s)), np.nan)
    for i, tam in enumerate(tams):
        for j, t1 in enumerate(t1s):
            v = agg.get((lx, tam, t1))
            if v:
                Z[i, j] = np.mean(v)
    return tams, t1s, Z


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("infile")
    ap.add_argument("-o", "--out", default="docs/anisotropy_heatmap.png")
    ap.add_argument("--col", type=int, default=10, help="CSV col: 10=chi_d vertex (susc), 18=equal-time vertex")
    ap.add_argument("--quantity", default=r"$\chi_d$ connected vertex")
    a = ap.parse_args()
    agg = load(a.infile, a.col)
    lxs = sorted({k[0] for k in agg})
    fig, axes = plt.subplots(1, len(lxs), figsize=(5.4 * len(lxs), 4.8), squeeze=False)
    vmax = max(abs(np.nanmin([np.nanmin(grid(agg, lx)[2]) for lx in lxs])),
              abs(np.nanmax([np.nanmax(grid(agg, lx)[2]) for lx in lxs])))
    for ax, lx in zip(axes[0], lxs):
        tams, t1s, Z = grid(agg, lx)
        im = ax.imshow(Z, origin="lower", aspect="auto", cmap="RdBu_r",
                       vmin=-vmax, vmax=vmax,
                       extent=[min(t1s) - 0.05, max(t1s) + 0.05,
                               min(tams) - 0.05, max(tams) + 0.05])
        # zero contour = pairing boundary
        T1, TAM = np.meshgrid(t1s, tams)
        try:
            cs = ax.contour(T1, TAM, Z, levels=[0.0], colors="k", linewidths=2.0, linestyles="--")
            ax.clabel(cs, fmt={0.0: "vertex = 0"}, fontsize=8)
        except Exception:
            pass
        for i, tam in enumerate(tams):
            for j, t1 in enumerate(t1s):
                if not np.isnan(Z[i, j]):
                    ax.text(t1, tam, f"{Z[i, j]:.2f}", ha="center", va="center",
                            fontsize=7.5, color="k")
        ax.set_xlabel(r"$t_1$  (NNN nesting-breaker)")
        ax.set_ylabel(r"$t_{am}$  (altermagnet NN anisotropy)")
        ax.set_title(f"{lx}×{lx}  half-filling, $U=6$")
        ax.set_xticks(t1s); ax.set_yticks(tams)
        fig.colorbar(im, ax=ax, label=a.quantity)
    fig.suptitle(a.quantity + r"  (red = pairing-favorable, blue = suppressed)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(a.out, dpi=150)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
