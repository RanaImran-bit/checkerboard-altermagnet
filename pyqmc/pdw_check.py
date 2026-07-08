#!/usr/bin/env python3
"""PDW check: does the finite-q peak of χ_d(q) GROW with β at fixed tam > 0?

A genuine PDW requires χ_d(q*) at finite q* to diverge as β→∞ (or at least
grow faster than the q=0 channel). So far at L≤8 this does not happen —
this script makes the check reproducible and explicit.

Demo (L=4, ~1 min):
    python pyqmc/pdw_check.py --lx 4 --ly 4 --tam 0.2 --betas 2,3,4
Cluster (L=8):
    python pyqmc/pdw_check.py --lx 8 --ly 8 --tam 0.2 --betas 3,4,5,6,8 --nmeas 300 --nw 20
Output CSV:
    python pyqmc/pdw_check.py --lx 6 --ly 6 --tam 0.2 --betas 3,4,5,6 --csv > results/pdw_L6_tam02.csv
"""
from __future__ import annotations
import os, sys, argparse
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "code", "ftcpqmc_py"))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "code", "dqmc_py"))
from ftcpmc import FTCPMC


def run_one(lx, ly, U, mu, beta, tam, nmeas, nw, seed):
    q = FTCPMC(lx, ly, U, mu, beta, dt=0.125, tam=tam, t1=0.0, seed=seed,
               nw=nw, constrained=True, stab=True)
    r = q.run_fb_stab(nmeas, kres=True, chi=True)
    sv = r["suscV_d"]
    pq = sv["Pq"]                           # (lx, ly) q-resolved vertex susceptibility
    maxk = float(pq.max())
    k0 = float(sv["k0"])                    # q=0 uniform channel
    ix = np.unravel_index(pq.argmax(), pq.shape)
    qx = 2 * np.pi * ix[0] / lx
    qy = 2 * np.pi * ix[1] / ly
    return dict(maxk=maxk, k0=k0, qx=qx, qy=qy,
                sign=r["sign"], dens=r["dens"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lx", type=int, default=4)
    ap.add_argument("--ly", type=int, default=4)
    ap.add_argument("--tam", type=float, default=0.2)
    ap.add_argument("--mu", type=float, default=2.0)    # half-filling: mu=U/2=2
    ap.add_argument("--U", type=float, default=4.0)
    ap.add_argument("--betas", type=str, default="2,3,4")
    ap.add_argument("--nmeas", type=int, default=120)
    ap.add_argument("--nw", type=int, default=12)
    ap.add_argument("--nseed", type=int, default=3)
    ap.add_argument("--csv", action="store_true")
    args = ap.parse_args()

    betas = [float(b) for b in args.betas.split(",")]

    if args.csv:
        print("beta,sign,dens,maxk,k0,maxk_over_k0,qx,qy")
    else:
        print(f"PDW check  L={args.lx}x{args.ly}  tam={args.tam}  U={args.U}  mu={args.mu}")
        print(f"{'beta':>5} {'sign':>6} {'dens':>6} {'chi(q*)':>10} {'chi(0)':>10} "
              f"{'ratio':>8} {'q*':>16}")
        print("-" * 70)

    for beta in betas:
        rows = [run_one(args.lx, args.ly, args.U, args.mu, beta, args.tam,
                        args.nmeas, args.nw, s) for s in range(1, args.nseed + 1)]
        avg = {k: np.mean([r[k] for r in rows]) for k in rows[0]}
        ratio = avg["maxk"] / avg["k0"] if avg["k0"] != 0 else float("nan")
        if args.csv:
            print(f"{beta:.1f},{avg['sign']:.3f},{avg['dens']:.4f},"
                  f"{avg['maxk']:.5f},{avg['k0']:.5f},{ratio:.4f},"
                  f"{avg['qx']:.4f},{avg['qy']:.4f}")
        else:
            q_str = f"({avg['qx']:.2f},{avg['qy']:.2f})"
            print(f"{beta:>5.1f} {avg['sign']:>6.3f} {avg['dens']:>6.4f} "
                  f"{avg['maxk']:>10.5f} {avg['k0']:>10.5f} {ratio:>8.4f} {q_str:>16}")

    if not args.csv:
        print()
        print("PDW verdict:")
        print("  chi(q*) GROWING with beta AND ratio > 1 → PDW candidate.")
        print("  chi(q*) flat/falling or ratio ≈ 1       → no PDW at this L/T.")


if __name__ == "__main__":
    main()
