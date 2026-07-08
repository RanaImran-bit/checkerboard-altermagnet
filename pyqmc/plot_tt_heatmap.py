#!/usr/bin/env python3
"""Two-panel (t1, tam) heatmap of the local d-wave vertex C_d(0)=N_pp(R=0) at two
fillings. Reads collector dump with ===CSV_N31=== / ===CSV_N26=== separators (or just
all pair_rspace --csv lines, grouped by N)."""
import sys, argparse, collections
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("infile")
    ap.add_argument("-o", "--out", default="docs/tt_heatmap_Cd0.png"); a = ap.parse_args()
    byN = collections.defaultdict(dict)   # N -> {(tam,t1): Cd0}
    for l in open(a.infile):
        p = l.split()
        if len(p) >= 14 and p[0].isdigit() and p[1].isdigit():
            lx = int(p[0]); N = int(p[2]); tam = round(float(p[5]), 3); t1 = round(float(p[6]), 3)
            byN[N][(tam, t1)] = float(p[8])
    Ns = sorted(byN)
    vmax = max(v for d in byN.values() for v in d.values())
    fig, axes = plt.subplots(1, len(Ns), figsize=(6.0 * len(Ns), 5.0), squeeze=False)
    for ax, N in zip(axes[0], Ns):
        d = byN[N]; tams = sorted({k[0] for k in d}); t1s = sorted({k[1] for k in d})
        Z = np.full((len(tams), len(t1s)), np.nan)
        for i, tam in enumerate(tams):
            for j, t1 in enumerate(t1s):
                if (tam, t1) in d: Z[i, j] = d[(tam, t1)]
        n = 2 * N / 64
        im = ax.imshow(Z, origin="lower", aspect="auto", cmap="inferno", vmin=0, vmax=vmax,
                       extent=[min(t1s) - 0.05, max(t1s) + 0.05, min(tams) - 0.05, max(tams) + 0.05])
        for i, tam in enumerate(tams):
            for j, t1 in enumerate(t1s):
                if not np.isnan(Z[i, j]):
                    ax.text(t1, tam, f"{Z[i, j]:.2f}", ha="center", va="center",
                            fontsize=7, color="white")
        ax.set_xlabel(r"$t_1$  (NNN nesting-breaker)", fontsize=12)
        ax.set_ylabel(r"$t_A$  (altermagnet NN)", fontsize=12)
        ax.set_title(rf"$n={n:.3f}$ (N={N})", fontsize=12)
        ax.set_xticks(t1s); ax.set_yticks(tams)
        fig.colorbar(im, ax=ax, label=r"$C_d(0)=N_{pp}(R{=}0)$")
    fig.suptitle(r"pyqmc local d-wave vertex $C_d(0)$ over $(t_1, t_A)$  (8×8, $U=4$)", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.96]); fig.savefig(a.out, dpi=150); print("wrote", a.out)


if __name__ == "__main__":
    main()
