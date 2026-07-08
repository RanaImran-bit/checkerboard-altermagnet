#!/usr/bin/env python3
"""Plot the matched R-resolved d-wave vertex vs tam (pyqmc, U=4, closed shell):
(a) local N_pp(R=0) [= manuscript C_d(0)] and (b) long-range <N_pp>_{|R|>2}.
Reads the collector CSV dump (lx ly nup ndn U tam t1 seed Npp_R0 e Npp_Rgt e qsum e)."""
import sys, argparse, collections, statistics as st, math
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("infile")
    ap.add_argument("-o", "--out", default="docs/rspace_vertex_vs_tam.png"); a = ap.parse_args()
    rows = [l.split() for l in open(a.infile) if l[:1].isdigit() and len(l.split()) >= 14]
    agg = collections.defaultdict(lambda: collections.defaultdict(list))
    for r in rows:
        lx = int(r[0]); tam = round(float(r[5]), 3)
        agg[lx][tam].append((float(r[8]), float(r[10])))   # R0, Rgt
    def ms(v): return st.mean(v), (st.pstdev(v) / math.sqrt(len(v)) if len(v) > 1 else 0.0)
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(10.5, 4.4))
    col = {6: "#1f77b4", 8: "#d62728"}
    for lx in sorted(agg):
        nn = {6: 17, 8: 31}[lx]; n = 2 * nn / (lx * lx)
        tams = sorted(agg[lx])
        r0 = [ms([x[0] for x in agg[lx][t]]) for t in tams]
        rg = [ms([x[1] for x in agg[lx][t]]) for t in tams]
        lab = rf"{lx}×{lx}, $n$={n:.3f}"
        ax0.errorbar(tams, [m for m, e in r0], yerr=[e for m, e in r0], fmt="o-",
                     c=col[lx], capsize=3, ms=5, lw=1.6, label=lab)
        ax1.errorbar(tams, [m for m, e in rg], yerr=[e for m, e in rg], fmt="s-",
                     c=col[lx], capsize=3, ms=5, lw=1.6, label=lab)
    ax0.set_xlabel(r"$t_A$", fontsize=12); ax0.set_ylabel(r"$N_{pp}(R{=}0)\;\equiv\;C_d(0)$", fontsize=12)
    ax0.set_title("(a) Local d-wave vertex", fontsize=12)
    ax1.set_xlabel(r"$t_A$", fontsize=12); ax1.set_ylabel(r"$\langle N_{pp}(R)\rangle_{|R|>2}$", fontsize=12)
    ax1.set_title("(b) Long-range d-wave vertex", fontsize=12)
    ax1.axhline(0, color="gray", lw=0.8, ls=":")
    for ax in (ax0, ax1): ax.legend(fontsize=9); ax.grid(alpha=0.3)
    fig.suptitle(r"pyqmc d-wave pairing vertex vs altermagnet anisotropy $t_A$  "
                 r"($U=4$, closed shell) — reproduces manuscript $C_d(0)$ increase", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.95]); fig.savefig(a.out, dpi=150)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
