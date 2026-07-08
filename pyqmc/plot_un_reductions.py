#!/usr/bin/env python3
"""Four figures from the k-resolved unified_scan CSV, each a 2x2 (rows: correlation /
susceptibility connected VERTEX; cols: the two dopings) of one reduction over (tam,t1):
  (a) maxk = max peak in k-space   (b) k0 = value at q=0
  (c) r0  = value at R=0           (d) rgt = average over |R|>2
CSV col map (0-idx): 5 tam, 6 t1, 14 corrV_maxk 15 k0 16 r0 17 rgt ; 22 suscV_maxk 23 k0 24 r0 25 rgt.

    python pyqmc/plot_un_reductions.py un.csv --outdir docs --tag t1 --xaxis t1
"""
import argparse, collections, os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

RED = {"maxk": ("max peak in k", 14, 22), "k0": ("value at q=0", 15, 23),
       "r0": ("value at R=0", 16, 24), "rgt": (r"avg over $|R|>2$", 17, 25)}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("infile")
    ap.add_argument("--outdir", default="docs"); ap.add_argument("--tag", default="t1")
    ap.add_argument("--xaxis", default="t1", help="x-axis knob label (t1 or tp)")
    a = ap.parse_args()
    # parse: key (N, tam, knob) -> [corrV..., suscV...]
    rows = collections.defaultdict(dict)
    for l in open(a.infile):
        p = l.split()
        if len(p) < 26 or not (p[0].isdigit() and p[1].isdigit()):
            continue
        N = int(p[2]); tam = round(float(p[5]), 3); knob = round(float(p[6]), 3)
        rows[N][(tam, knob)] = [float(x) for x in p]
    Ns = sorted(rows, reverse=True)
    dens = {N: 2 * N / 64 for N in Ns}
    for rk, (label, ci_c, ci_s) in RED.items():
        tams = sorted({k[0] for k in rows[Ns[0]]}); knobs = sorted({k[1] for k in rows[Ns[0]]})
        fig, axes = plt.subplots(2, len(Ns), figsize=(5.6 * len(Ns), 9.0), squeeze=False)
        panels = [("correlation vertex", ci_c), ("susceptibility vertex", ci_s)]
        for ri, (rowlab, ci) in enumerate(panels):
            Zs = {}
            for N in Ns:
                Z = np.full((len(tams), len(knobs)), np.nan)
                for i, tam in enumerate(tams):
                    for j, kn in enumerate(knobs):
                        if (tam, kn) in rows[N]: Z[i, j] = rows[N][(tam, kn)][ci]
                Zs[N] = Z
            vmax = max(np.nanmax(np.abs(Z)) for Z in Zs.values()) or 1.0
            for ci2, N in enumerate(Ns):
                ax = axes[ri][ci2]; Z = Zs[N]
                im = ax.imshow(Z, origin="lower", aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                               extent=[min(knobs) - .05, max(knobs) + .05, min(tams) - .05, max(tams) + .05])
                for i, tam in enumerate(tams):
                    for j, kn in enumerate(knobs):
                        if not np.isnan(Z[i, j]):
                            ax.text(kn, tam, f"{Z[i, j]:.2f}", ha="center", va="center", fontsize=6.5, color="k")
                ax.set_xlabel(rf"${a.xaxis}$"); ax.set_ylabel(r"$t_A$")
                ax.set_title(rf"{rowlab}  |  $n={dens[N]:.3f}$", fontsize=10)
                ax.set_xticks(knobs); ax.set_yticks(tams); fig.colorbar(im, ax=ax)
        fig.suptitle(rf"d-wave connected vertex — {label}  over $(t_A,\,{a.xaxis})$  (8×8, $U=4$)", fontsize=13)
        fig.tight_layout(rect=[0, 0, 1, 0.97])
        out = os.path.join(a.outdir, f"un_{a.tag}_{rk}.png")
        fig.savefig(out, dpi=150); plt.close(fig); print("wrote", out)


if __name__ == "__main__":
    main()
