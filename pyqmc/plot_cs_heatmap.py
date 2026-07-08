#!/usr/bin/env python3
"""2x2 (tam, t1) heatmaps from the run DB (campaign with chid data at two fillings):
rows = {equal-time connected vertex Cd0_vtx, susceptibility chi_d_vtx}, cols = fillings.
Both are the q=0-sum connected d-wave vertex (consistent projection).

    python pyqmc/plot_cs_heatmap.py --campaign 7 -o docs/cs_heatmaps.png
"""
import argparse, sqlite3, os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

DB = os.path.join(os.path.dirname(__file__), "..", "results", "runs.db")


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--campaign", type=int, required=True)
    ap.add_argument("-o", "--out", default="docs/cs_heatmaps.png"); a = ap.parse_args()
    c = sqlite3.connect(DB)
    rows = c.execute("SELECT nup, density, tam, t1, AVG(chi_d_vtx), AVG(Cd0_vtx) "
                     "FROM points WHERE campaign_id=? GROUP BY nup,tam,t1", (a.campaign,)).fetchall()
    Ns = sorted({r[0] for r in rows}, reverse=True)   # higher N (near half) first
    d = {(N, tam, t1): (chi, cd0) for N, n, tam, t1, chi, cd0 in rows}
    dens = {r[0]: r[1] for r in rows}
    tams = sorted({r[2] for r in rows}); t1s = sorted({r[3] for r in rows})
    quant = [("Cd0_vtx", 1, "equal-time vertex correlation  $C_d(0)^{\\rm vtx}$"),
             ("chi_d_vtx", 0, "vertex susceptibility  $\\chi_d^{\\rm vtx}$")]
    fig, axes = plt.subplots(2, len(Ns), figsize=(5.6 * len(Ns), 9.0), squeeze=False)
    for ri, (key, idx, label) in enumerate(quant):
        Zs = {}
        for N in Ns:
            Z = np.full((len(tams), len(t1s)), np.nan)
            for i, tam in enumerate(tams):
                for j, t1 in enumerate(t1s):
                    if (N, tam, t1) in d: Z[i, j] = d[(N, tam, t1)][idx]
            Zs[N] = Z
        vmax = max(np.nanmax(np.abs(Z)) for Z in Zs.values())
        for ci, N in enumerate(Ns):
            ax = axes[ri][ci]; Z = Zs[N]
            im = ax.imshow(Z, origin="lower", aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                           extent=[min(t1s) - .05, max(t1s) + .05, min(tams) - .05, max(tams) + .05])
            for i, tam in enumerate(tams):
                for j, t1 in enumerate(t1s):
                    if not np.isnan(Z[i, j]):
                        ax.text(t1, tam, f"{Z[i, j]:.2f}", ha="center", va="center", fontsize=6.5, color="k")
            ax.set_xlabel(r"$t_1$"); ax.set_ylabel(r"$t_A$")
            ax.set_title(rf"{label}  |  $n={dens[N]:.3f}$ (N={N})", fontsize=10)
            ax.set_xticks(t1s); ax.set_yticks(tams)
            fig.colorbar(im, ax=ax)
    fig.suptitle(r"q=0 connected d-wave vertex over $(t_1, t_A)$  (8×8, $U=4$)  "
                 r"— red=positive, blue=negative", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.97]); fig.savefig(a.out, dpi=150); print("wrote", a.out)


if __name__ == "__main__":
    main()
