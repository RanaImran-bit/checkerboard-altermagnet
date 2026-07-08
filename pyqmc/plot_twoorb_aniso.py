#!/usr/bin/env python3
"""Figure: two-orbital altermagnet d-wave (and ext-s) pairing VERTEX susceptibility vs
altermagnetic anisotropy alpha = 1 - t2/t1, across methods (CPQMC T=0 full model;
DQMC / CP-DQMC finite-T reduced U-only model). Reads the campaign CSVs in
results/dqmc_scan/twoorb_*.csv and overlays the curves.

  python pyqmc/plot_twoorb_aniso.py -o results/dqmc_scan/twoorb_aniso.png
"""
from __future__ import annotations
import os, sys, csv, argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load(path):
    if not os.path.exists(path):
        return None
    with open(path) as f:
        rows = list(csv.DictReader(f))
    return rows


def col(rows, key):
    return np.array([float(r[key]) for r in rows if r.get(key) not in (None, "", "nan")])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="results/dqmc_scan")
    ap.add_argument("-o", "--out", default="results/dqmc_scan/twoorb_aniso.png")
    a = ap.parse_args()
    D = a.dir
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))

    # --- panel 1: d-wave vertex susceptibility vs alpha, all methods ---
    ax = axes[0]
    series = [
        ("twoorb_cpqmc_L4_half.csv", "chi_d_vertex", "chi_d_vertex_err", "CPQMC T=0  L4 full model", "o-", "C0"),
        ("twoorb_cpqmc_L3_half_full.csv", "chi_d_vertex", "chi_d_vertex_err", "CPQMC T=0  L3 full model", "s--", "C1"),
    ]
    for fn, yk, ek, lab, fmt, c in series:
        rows = load(os.path.join(D, fn))
        if not rows:
            continue
        al = col(rows, "alpha"); y = col(rows, yk)
        try:
            e = col(rows, ek)
        except Exception:
            e = None
        ax.errorbar(al, y, yerr=e if e is not None and len(e) == len(y) else None,
                    fmt=fmt, color=c, label=lab, capsize=3)
    # finite-T DQMC / CP-DQMC reduced
    for fn, lab, fmt, c in (("twoorb_ft_dqmc_L4_half.csv", "DQMC  L4 U-only (sign=1)", "^-", "C2"),
                            ("twoorb_ft_cpdqmc_L4_half.csv", "CP-DQMC  L4 U-only", "v-", "C3")):
        rows = load(os.path.join(D, fn))
        if not rows:
            continue
        rows = [r for r in rows if r.get("method", "").startswith(fn.split("_")[2][:2]) or True]
        al = col(rows, "alpha"); y = col(rows, "chi_d_vertex")
        ax.plot(al, y, fmt, color=c, label=lab)
    ax.axhline(0, color="k", lw=0.6, ls=":")
    ax.set_xlabel(r"altermagnetic anisotropy  $\alpha = 1 - t_2/t_1$")
    ax.set_ylabel(r"$d_{x^2-y^2}$ pairing vertex susceptibility  $\chi_d^{\rm vertex}$")
    ax.set_title("d-wave pairing vertex vs anisotropy (half-filling)")
    ax.legend(fontsize=8)

    # --- panel 2: d vs ext-s (DQMC exact, half-filling) ---
    ax = axes[1]
    rows = load(os.path.join(D, "twoorb_ft_dqmc_L4_half.csv"))
    if rows:
        al = col(rows, "alpha")
        ax.plot(al, col(rows, "chi_d_vertex"), "o-", color="C0", label=r"$\chi_d^{\rm vertex}$ (d-wave)")
        ax.plot(al, col(rows, "chi_s_vertex"), "s-", color="C3", label=r"$\chi_s^{\rm vertex}$ (ext-s)")
        ax.set_xlabel(r"anisotropy  $\alpha$")
        ax.set_ylabel(r"vertex susceptibility")
        ax.set_title("DQMC L4 U-only (sign=1, exact):  d $\\gg$ ext-s, d peaks at weak $\\alpha$")
        ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(a.out, dpi=140)
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
