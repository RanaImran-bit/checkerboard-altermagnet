#!/usr/bin/env python3
"""Benchmark figure for the momentum-resolved magnetic susceptibility chi_zz(q)
(docs/PLAN_chi_spin_AHE.md W5): reads results/chi_spin_bench/*.json and renders

  (a) finite-T 2x2: chi_zz(q) per q -- ED vs DQMC vs CP-DQMC(free) vs CP-DQMC(CP)
  (b) finite-T 2x2 tam=0 vs tam=0.3: the tam effect on the chi_zz(q) grid (ED + DQMC)
  (c) T=0 CPQMC vs ED: windowed chi_zz(q) per q (closed-shell gates)
  (d) T=0 2x2 AM point: C_q(tau) decay curves, CPMC (points+err) vs ED Lehmann (lines)

    python pyqmc/plot_chi_spin_bench.py            # -> docs/chi_spin_bench.png
"""
from __future__ import annotations
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.join(os.path.dirname(__file__), "..")
BEN = os.path.join(ROOT, "results", "chi_spin_bench")
OUT = os.path.join(ROOT, "docs", "chi_spin_bench.png")


def load(name):
    p = os.path.join(BEN, name)
    if not os.path.exists(p):
        return None
    with open(p) as f:
        return json.load(f)


def qlabels(lx, ly):
    return [f"({kx},{ky})" for kx in range(lx) for ky in range(ly)]


def main():
    fig, axes = plt.subplots(2, 2, figsize=(11, 8.5))
    (ax_a, ax_b), (ax_c, ax_d) = axes

    # ---- (a) finite-T three-way, 2x2 tam=0.3 ----------------------------------
    ft = load("ft_2x2_tam03.json") or load("ft_2x2_tam0.json")
    if ft:
        lx, ly = ft["params"]["lx"], ft["params"]["ly"]
        labs = qlabels(lx, ly); x = np.arange(len(labs))
        w = 0.2
        for k, (tag, off) in enumerate((("ED", -1.5), ("DQMC", -0.5),
                                        ("CPfree", 0.5), ("CPcons", 1.5))):
            if tag not in ft["chi"]:
                continue
            v = np.array(ft["chi"][tag]).reshape(-1)
            ax_a.bar(x + off * w, v, w, label=tag)
        ax_a.set_xticks(x); ax_a.set_xticklabels(labs)
        ax_a.set_ylabel(r"$\chi_{zz}(q)$ per site")
        p = ft["params"]
        ax_a.set_title(f"(a) finite-T three-way  {lx}x{ly}  U={p['U']} "
                       f"$\\beta$={p['beta']} tam={p['tam']}  "
                       f"[CP $\\langle s\\rangle$={ft['signs'].get('CPcons', 1):.3f}]")
        ax_a.legend(fontsize=8)

    # ---- (b) tam effect on the ED/DQMC grids ----------------------------------
    f0 = load("ft_2x2_tam0.json"); f3 = load("ft_2x2_tam03.json")
    if f0 and f3:
        labs = qlabels(2, 2); x = np.arange(len(labs))
        for ftag, r, mk in (("tam=0", f0, "o"), ("tam=0.3", f3, "s")):
            ax_b.plot(x, np.array(r["chi"]["ED"]).reshape(-1), mk + "-",
                      label=f"ED {ftag}")
            ax_b.plot(x, np.array(r["chi"]["DQMC"]).reshape(-1), mk + "--", alpha=0.6,
                      label=f"DQMC {ftag}")
        ax_b.set_xticks(x); ax_b.set_xticklabels(labs)
        ax_b.set_ylabel(r"$\chi_{zz}(q)$ per site")
        ax_b.set_title("(b) altermagnet knob: (pi,pi) response vs tam (2x2, $\\beta$=2)")
        ax_b.legend(fontsize=8)

    # ---- (c) T=0 CPQMC vs ED (closed-shell gates) ------------------------------
    shown = False
    for name, lab in (("t0v2_4x2_U0_cs.json", "4x2 U=0 (exact gate)"),
                      ("t0v2_4x2_U4_cs.json", "4x2 U=4 closed shell"),
                      ("t0v2_4x4_U3.json", "4x4 U=3 (validate_chi pt)")):
        r = load(name)
        if not r:
            continue
        p = r["params"]; lx, ly = p["lx"], p["ly"]
        ed = np.array(r["ed_chi_win"]).reshape(-1)
        qc = np.array(r["qmc"]["chi_q"]).reshape(-1)
        qe = np.array(r["qmc"]["chi_q_err"]).reshape(-1)
        ln, = ax_c.plot(range(len(ed)), ed, "-", lw=1.2, label=f"ED {lab}")
        ax_c.errorbar(range(len(qc)), qc, yerr=qe, fmt="o", ms=4,
                      color=ln.get_color(), label=f"CPMC {lab}")
        shown = True
    if shown:
        ax_c.set_xlabel("q index (kx*ly+ky)")
        ax_c.set_ylabel(r"windowed $\chi_{zz}(q)$")
        ax_c.set_title("(c) T=0 CPQMC vs Lehmann ED (windowed Kubo)")
        ax_c.legend(fontsize=7)

    # ---- (d) C_q(tau) curves, T=0 2x2 AM point ---------------------------------
    r = load("t0v2_2x2_U4_am.json") or load("t0_2x2_U4_am.json")
    if r:
        p = r["params"]; lx, ly = p["lx"], p["ly"]
        taus = np.array(r["qmc"]["taus"])
        Cq = np.array(r["qmc"]["Ctau_q"]); Ce = np.array(r["qmc"]["Cerr_q"])
        Ced = np.array(r["ed_Ctau_q"])
        for (kx, ky) in ((1, 1), (0, 1), (1, 0)):
            ln, = ax_d.plot(taus, Ced[:, kx, ky], "-", lw=1.2,
                            label=f"ED q=({kx},{ky})")
            ax_d.errorbar(taus[::2], Cq[::2, kx, ky], yerr=Ce[::2, kx, ky],
                          fmt="o", ms=3, color=ln.get_color())
        ax_d.set_xlabel(r"$\tau$"); ax_d.set_ylabel(r"$C_q(\tau)$ per site")
        ax_d.set_title(f"(d) T=0 CPMC $C_q(\\tau)$ vs ED  "
                       f"({lx}x{ly} U={p['U']} tam={p['tam']} t1={p['t1']})")
        ax_d.legend(fontsize=8)

    fig.suptitle("chi_zz(q) benchmark: ED vs DQMC vs CP-DQMC vs CPQMC "
                 "(results/chi_spin_bench)", y=0.995)
    fig.tight_layout()
    fig.savefig(OUT, dpi=160)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
