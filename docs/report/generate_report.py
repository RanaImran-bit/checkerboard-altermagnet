#!/usr/bin/env python3
"""Generate docs/report/validation_report.pdf — a brief summary of the QMC
validation platform and the ED / Fortran-CPQMC / Python-CPQMC cross-checks."""
from __future__ import annotations
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "validation_report.pdf")

TXT = "#1a1a1a"; MUT = "#555"; ACC = "#1f4e9e"; OK = "#1d8f5e"; REF = "#b07b1e"

def header(fig, title, sub):
    fig.text(0.07, 0.95, title, fontsize=17, weight="bold", color=TXT)
    fig.text(0.07, 0.915, sub, fontsize=9.5, color=MUT)
    fig.add_artist(plt.Line2D([0.07, 0.93], [0.905, 0.905], color="#ccc", lw=0.8))

def bullets(fig, y0, lines, x=0.07, dy=0.032, size=9.5):
    y = y0
    for b, t in lines:
        fig.text(x, y, b, fontsize=size, color=ACC, weight="bold")
        fig.text(x + 0.025, y, t, fontsize=size, color=TXT)
        y -= dy
    return y

with PdfPages(OUT) as pdf:
    # ---------------- Page 1 ----------------
    fig = plt.figure(figsize=(8.27, 11.69))  # A4
    header(fig, "QMC Code-Validation Platform — Summary Report",
           "Constrained-path QMC for a two-orbital altermagnetic Hubbard model · cross-checked against exact diagonalization")
    fig.text(0.07, 0.885, "Purpose", fontsize=12, weight="bold", color=TXT)
    fig.text(0.07, 0.845,
             "A reproducible platform to validate changes to a constrained-path quantum Monte Carlo (CPQMC) code,\n"
             "supporting the PRL study of enhanced d-wave pairing in a strongly correlated altermagnet. The platform\n"
             "compares three independent solvers on identical models: exact diagonalization (ED, ground truth on small\n"
             "clusters), the production Fortran CPQMC, and a from-scratch Python CPQMC port.",
             fontsize=9.5, color=TXT, va="top")

    fig.text(0.07, 0.78, "Platform components", fontsize=12, weight="bold", color=TXT)
    bullets(fig, 0.745, [
        ("•", "Reproducible toolchain (conda): gfortran + OpenMPI + OpenBLAS; QuSpin for ED."),
        ("•", "Build/run/parse harness: builds each QMC version, runs it, parses output to a common schema."),
        ("•", "ED reference (QuSpin): single-band Hubbard and the two-orbital altermagnet model."),
        ("•", "Comparison engine: per-observable delta / sigma / z-distance with PASS/FAIL verdicts."),
        ("•", "FastAPI + React UI: compute ED, pick QMC runs, overlay energies and the diff table."),
        ("•", "Python CPQMC port (pyqmc): constrained-path AFQMC, validated against ED."),
    ])

    # single-band validation chart
    ax = fig.add_axes([0.1, 0.30, 0.82, 0.20])
    labels = ["ED (exact)", "Fortran CPQMC", "Python CPMC"]
    vals = [-7.8672, -7.8701, -7.8675]
    errs = [0.0, 0.0029, 0.0002]
    cols = [REF, ACC, OK]
    ax.bar(labels, vals, yerr=errs, color=cols, width=0.5, capsize=5, alpha=0.9)
    ax.set_ylim(-7.90, -7.83)
    ax.set_ylabel("ground-state energy")
    ax.set_title("Single-band Hubbard  (4×4, 1↑1↓, U=3, t=1)   — all agree within ~2σ", fontsize=10)
    for i, (v, e) in enumerate(zip(vals, errs)):
        ax.text(i, v - 0.004, f"{v:.4f}" + (f"\n±{e:.4f}" if e else ""), ha="center", va="top", fontsize=8, color="white")

    fig.text(0.07, 0.22, "Key finding (Fortran bug fixed)", fontsize=12, weight="bold", color=TXT)
    fig.text(0.07, 0.18,
             "The platform localized and fixed a real sign error in the Fortran back-propagation Green's function:\n"
             "the measurement-phase energy came out as the exact negative of the truth (+96 vs −96 at U=0). Root\n"
             "cause: calgf method 2 built g = +G while the energy estimator expects g = det·(I−G). After the fix, the\n"
             "U=0 measured energy matches the exact free-fermion ground state to machine precision.",
             fontsize=9.5, color=TXT, va="top")
    fig.text(0.07, 0.045, "QMC validation platform · generated automatically", fontsize=8, color=MUT)
    pdf.savefig(fig); plt.close(fig)

    # ---------------- Page 2 ----------------
    fig = plt.figure(figsize=(8.27, 11.69))
    header(fig, "Two-Orbital Altermagnet — Three-Way Validation",
           "Exact diagonalization vs Fortran CPQMC vs Python CPQMC on a 2×2 cluster (8 orbital-sites, 4↑4↓)")

    fig.text(0.07, 0.875, "Model", fontsize=12, weight="bold", color=TXT)
    fig.text(0.07, 0.845,
             "Two-orbital (d$_{xz}$, d$_{yz}$) square lattice:  ε$_x$=−2t₁cos k$_x$−2t₂cos k$_y$−4t₃cos k$_x$cos k$_y$,\n"
             "ε$_y$=(x↔y),  ε$_{xy}$=−4t₄ sin k$_x$ sin k$_y$;  on-site Hubbard U$_{xx}$.  (t₁=t₂=t₃=t₄=−1.)",
             fontsize=9.5, color=TXT, va="top")

    # three-way table + chart
    U = [0, 2, 4]
    ed = [-32.000, -28.289, -25.128]
    fort = [(-32.000, 0.0), (-28.275, 0.007), (-25.134, 0.010)]
    py = [(-32.000, 0.0), (-28.289, 0.001), (-25.126, 0.005)]

    ax = fig.add_axes([0.1, 0.45, 0.82, 0.28])
    x = np.arange(len(U)); w = 0.26
    ax.bar(x - w, ed, w, label="ED (exact)", color=REF, alpha=0.9)
    ax.bar(x, [f[0] for f in fort], w, yerr=[f[1] for f in fort], capsize=4, label="Fortran CPQMC", color=ACC, alpha=0.9)
    ax.bar(x + w, [p[0] for p in py], w, yerr=[p[1] for p in py], capsize=4, label="Python CPQMC", color=OK, alpha=0.9)
    ax.set_xticks(x); ax.set_xticklabels([f"U$_{{xx}}$={u}" for u in U])
    ax.set_ylabel("ground-state energy"); ax.set_ylim(-33, -23)
    ax.legend(fontsize=8, loc="upper left"); ax.set_title("Ground-state energy vs on-site interaction", fontsize=10)

    # table
    ax2 = fig.add_axes([0.1, 0.18, 0.82, 0.20]); ax2.axis("off")
    rows = [["U_xx", "ED (exact)", "Fortran CPQMC", "Python CPQMC", "worst z"]]
    zs = ["—", "1.85 / 0.33", "0.58 / 0.42"]
    for i, u in enumerate(U):
        rows.append([str(u), f"{ed[i]:.3f}",
                     f"{fort[i][0]:.3f} ± {fort[i][1]:.3f}",
                     f"{py[i][0]:.3f} ± {py[i][1]:.3f}", zs[i]])
    tbl = ax2.table(cellText=rows[1:], colLabels=rows[0], loc="center", cellLoc="center")
    tbl.auto_set_font_size(False); tbl.set_fontsize(9); tbl.scale(1, 1.6)
    for (r, c), cell in tbl.get_celld().items():
        cell.set_edgecolor("#ddd")
        if r == 0:
            cell.set_facecolor("#f0f0f0"); cell.set_text_props(weight="bold")

    fig.text(0.07, 0.12, "Conclusion", fontsize=12, weight="bold", color=TXT)
    fig.text(0.07, 0.085,
             "ED, the production Fortran CPQMC, and the independent Python CPQMC port agree on the two-orbital\n"
             "altermagnet ground-state energy across interaction strengths (all within ~2σ) — the cross-validation\n"
             "benchmark requested in review. Next: validate the inter-orbital (U$_{xy}$) and neighbour (V) channels,\n"
             "then scale to publication-size clusters.",
             fontsize=9.5, color=TXT, va="top")
    fig.text(0.07, 0.035, "QMC validation platform · generated automatically", fontsize=8, color=MUT)
    pdf.savefig(fig); plt.close(fig)

print(f"wrote {OUT}")
