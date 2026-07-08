#!/usr/bin/env python3
"""Finite-T two-orbital (d_xz,d_yz) altermagnet d-wave / ext-s pairing SUSCEPTIBILITY
vs anisotropy alpha = 1 - t2/t1, via DQMC (free, sign-tracked) AND CP-DQMC (constrained,
sign-controlled). REDUCED interaction (on-site U only -- DQMC/CP-DQMC do not implement
the inter-orbital uxy / neighbour v of the full manuscript model; that is CPQMC's job).
With t3=t4=0 the U-only two-orbital model is two decoupled anisotropic Hubbard bands and
is sign-problem-free at half-filling -- a clean finite-T cross-check of the trend.

  chi_a = int_0^beta <Delta_a(tau) Delta_a^dag(0)> dtau,  vertex = chi_a - bubble (a=s,d)

  python code/dqmc_py/two_orb_ft_scan.py --lx 4 --ly 4 --U 4 --beta 5 \
      --method both --t2list -1.0 -0.6 -0.2 -o results/dqmc_scan/twoorb_ft_L4.csv
"""
from __future__ import annotations
import os, sys, csv, argparse
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "ftcpqmc_py"))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "ed"))
from dqmc import DQMC
from altermagnet_ed import build_hopping


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lx", type=int, default=4); ap.add_argument("--ly", type=int, default=4)
    ap.add_argument("--U", type=float, default=4.0); ap.add_argument("--beta", type=float, default=5.0)
    ap.add_argument("--dt", type=float, default=0.125)
    ap.add_argument("--t1", type=float, default=-1.0); ap.add_argument("--t3", type=float, default=0.0)
    ap.add_argument("--t4", type=float, default=0.0)
    ap.add_argument("--uxy", type=float, default=0.0); ap.add_argument("--v", type=float, default=0.0)
    ap.add_argument("--method", choices=["dqmc", "cpdqmc", "both"], default="both")
    ap.add_argument("--nwarm", type=int, default=120); ap.add_argument("--nmeas", type=int, default=400)
    ap.add_argument("--nw", type=int, default=48)           # CP-DQMC walkers
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--t2list", type=float, nargs="+", default=[-1.0, -0.6, -0.2])
    ap.add_argument("-o", "--out", default="")
    a = ap.parse_args()
    nsites = 2 * a.lx * a.ly
    mu = 0.5 * a.U                                          # half-filling (PH symmetric point)
    rows = []
    print(f"# two-orbital finite-T {a.lx}x{a.ly} nsites={nsites} U={a.U} uxy={a.uxy} v={a.v} "
          f"beta={a.beta} t1={a.t1} t3={a.t3} t4={a.t4} mu=U/2")
    print(f"# {'method':>7} {'alpha':>6} {'dens':>6} {'sign':>6} {'chi_d':>9} {'chi_dV':>9} {'chi_s':>9} {'chi_sV':>9}")
    for t2 in a.t2list:
        alpha = 1 - t2 / a.t1
        K = build_hopping(a.lx, a.ly, a.t1, t2, a.t3, a.t4)
        res = {}
        if a.method in ("dqmc", "both"):
            q = DQMC(a.lx, a.ly, a.U, mu, a.beta, a.dt, Kmat=K, seed=a.seed, uxy=a.uxy, v=a.v)
            r = q.run(nwarm=a.nwarm, nmeas=a.nmeas, chi=True)
            res["dqmc"] = r
            print(f"  {'dqmc':>7} {alpha:>6.2f} {r['dens']:>6.3f} {r['sign']:>6.3f} "
                  f"{r['chid']:>9.4f} {r['chidV']:>9.4f} {r['chis']:>9.4f} {r['chisV']:>9.4f}")
            sys.stdout.flush()
        if a.method in ("cpdqmc", "both"):
            from ftcpmc import FTCPMC
            q = FTCPMC(a.lx, a.ly, a.U, mu, a.beta, a.dt, Kmat=K, seed=a.seed, nw=a.nw, stab=True,
                       uxy=a.uxy, v=a.v)
            r = q.run_fb_stab(nmeas=a.nmeas, chi=True)
            res["cpdqmc"] = r
            print(f"  {'cpdqmc':>7} {alpha:>6.2f} {r['dens']:>6.3f} {r['sign']:>6.3f} "
                  f"{r['chid']:>9.4f} {r.get('chidV', float('nan')):>9.4f} {r['chis']:>9.4f} {r.get('chisV', float('nan')):>9.4f}")
            sys.stdout.flush()
        for meth, r in res.items():
            rows.append(dict(method=meth, t2=t2, alpha=alpha, nsites=nsites, U=a.U, uxy=a.uxy, v=a.v,
                             beta=a.beta, t3=a.t3, t4=a.t4, dens=r["dens"], sign=r["sign"],
                             chi_d=r["chid"], chi_d_vertex=r.get("chidV", float("nan")),
                             chi_s=r["chis"], chi_s_vertex=r.get("chisV", float("nan")),
                             Sd=r["Sd"], SdV=r.get("SdV", float("nan")),
                             Ss=r["Ss"], SsV=r.get("SsV", float("nan")),
                             S_AM=r.get("S_AM", float("nan")), energy=r["energy"]))
    if a.out:
        os.makedirs(os.path.dirname(a.out), exist_ok=True)
        fields = list(rows[0].keys())
        with open(a.out, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields); w.writeheader()
            for row in rows:
                w.writerow(row)
        print(f"# wrote {a.out}")


if __name__ == "__main__":
    main()
