#!/usr/bin/env python3
"""CP-DQMC τ-integrated pairing eigenvalue at SC-dome doping (n≈0.87).

Sign-controlled (constraint → sign=1) check: does d-wave lead in the
τ-integrated pairing eigenvalue at the doping where SC actually lives?
This is the cleanest positive-or-negative SC statement — DQMC can't reach
this regime cleanly due to the sign problem.

Step 1 (auto): binary-search μ to hit target filling for each (β, tam).
Step 2: measure τ-integrated pair matrix → leading eigenvalue + d/s overlap.

Demo (L=4, ~2 min):
    python code/ftcpqmc_py/cpdqmc_dome_scan.py --lx 4 --ly 4
Cluster (L=8):
    python code/ftcpqmc_py/cpdqmc_dome_scan.py --lx 8 --ly 8 --betas 4,5,6,8 --nmeas 300 --nw 20
Output CSV:
    python code/ftcpqmc_py/cpdqmc_dome_scan.py --lx 6 --ly 6 --csv > results/cpdqmc_dome_L6.csv
"""
from __future__ import annotations
import os, sys, argparse
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "dqmc_py"))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "pyqmc"))
from ftcpmc import FTCPMC


def _density(lx, ly, U, mu, beta, tam, seed=42):
    q = FTCPMC(lx, ly, U, mu, beta, dt=0.125, tam=tam, t1=0.0, seed=seed,
               nw=6, constrained=True, stab=True)
    return q.run_fb_stab(30, kres=False, chi=False)["dens"]


def find_mu(lx, ly, U, beta, tam, target_n, tol=0.01):
    """Binary search μ ∈ [0.5, U/2] for the target filling."""
    lo, hi = 0.5, U / 2
    mu = 0.5 * (lo + hi)
    n = target_n
    for _ in range(16):
        mu = 0.5 * (lo + hi)
        n = _density(lx, ly, U, mu, beta, tam)
        if abs(n - target_n) < tol:
            return mu, n
        if n < target_n:
            lo = mu
        else:
            hi = mu
    return mu, n


def run_point(lx, ly, U, mu, beta, tam, nmeas, nw, seed):
    q = FTCPMC(lx, ly, U, mu, beta, dt=0.125, tam=tam, t1=0.0, seed=seed,
               nw=nw, constrained=True, stab=True)
    q.paireig_tau = True
    r = q.run_fb_stab(nmeas, kres=False, chi=True)
    return dict(dens=r["dens"], sign=r["sign"],
                pe_lam=r.get("pe_lam", float("nan")),
                pe_dov=r.get("pe_dov", float("nan")),
                pe_sov=r.get("pe_sov", float("nan")))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lx", type=int, default=4)
    ap.add_argument("--ly", type=int, default=4)
    ap.add_argument("--U", type=float, default=4.0)
    ap.add_argument("--target-n", type=float, default=0.875, dest="target_n")
    ap.add_argument("--tams", type=str, default="0.0,0.1,0.2,0.3")
    ap.add_argument("--betas", type=str, default="3,4,5")
    ap.add_argument("--nmeas", type=int, default=120)
    ap.add_argument("--nw", type=int, default=12)
    ap.add_argument("--nseed", type=int, default=2)
    ap.add_argument("--csv", action="store_true")
    args = ap.parse_args()

    betas = [float(b) for b in args.betas.split(",")]
    tams = [float(t) for t in args.tams.split(",")]

    if args.csv:
        print("beta,tam,mu,dens,sign,pe_lam,pe_dov,pe_sov,verdict")
    else:
        print(f"CP-DQMC dome pairing eigenvalue  L={args.lx}x{args.ly}  "
              f"target_n={args.target_n}  U={args.U}")
        print(f"{'beta':>5} {'tam':>5} {'mu':>6} {'dens':>6} {'sign':>6} "
              f"{'λ_τ':>9} {'d_ov':>7} {'s_ov':>7} verdict")
        print("-" * 72)

    for beta in betas:
        for tam in tams:
            mu, _ = find_mu(args.lx, args.ly, args.U, beta, tam, args.target_n)
            rows = [run_point(args.lx, args.ly, args.U, mu, beta, tam,
                              args.nmeas, args.nw, s)
                    for s in range(1, args.nseed + 1)]
            avg = {k: np.mean([r[k] for r in rows]) for k in rows[0]}
            verdict = "d-WAVE" if avg["pe_dov"] > avg["pe_sov"] else "ext-s"
            if args.csv:
                print(f"{beta:.1f},{tam:.2f},{mu:.4f},{avg['dens']:.4f},"
                      f"{avg['sign']:.3f},{avg['pe_lam']:.5f},"
                      f"{avg['pe_dov']:.3f},{avg['pe_sov']:.3f},{verdict}")
            else:
                print(f"{beta:>5.1f} {tam:>5.2f} {mu:>6.3f} {avg['dens']:>6.4f} "
                      f"{avg['sign']:>6.3f} {avg['pe_lam']:>9.5f} "
                      f"{avg['pe_dov']:>7.3f} {avg['pe_sov']:>7.3f} {verdict}")

    if not args.csv:
        print()
        print("Verdict: d-WAVE means d-wave eigenvector leads (d_ov > s_ov).")
        print("         λ_τ growing with β at fixed tam → coherent SC pairing.")


if __name__ == "__main__":
    main()
