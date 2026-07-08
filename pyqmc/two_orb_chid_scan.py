#!/usr/bin/env python3
"""Two-orbital (d_xz,d_yz) altermagnet: CPQMC (T=0 constrained-path) d-wave / ext-s
PAIRING SUSCEPTIBILITY vs the altermagnetic anisotropy alpha = 1 - t2/t1.

This is the LARGER-LATTICE companion to the exact ED scan (ed_aniso_scan.py): same
two-orbital model (build_hopping t1,t2,t3,t4) with the FULL interaction set
(on-site uxx, inter-orbital uxy, neighbour v), measured with the back-propagated
unequal-time pairing susceptibility run_bp_chid (full + connected VERTEX, s & d).

  chi_a = int_0^inf <Delta_a(tau) Delta_a^dag(0)> dtau    (a = s, d_{x2-y2})
  C_a(tau0) = equal-time pair structure factor (consistency vs ED S_a)

Output: one CSV row per anisotropy point with chi_d/chi_d_vertex/chi_s/... + errors.

  python pyqmc/two_orb_chid_scan.py --lx 4 --ly 4 --nup 8 --ndn 8 \
      --U 1 --uxy 1 --v 0.1 --t3 -1 --t4 -1 --bp 24 --nw 480 \
      --t2list -1.0 -0.8 -0.6 -0.4 -0.2 -o results/dqmc_scan/twoorb_cpqmc_L4.csv
"""
from __future__ import annotations
import os, sys, csv, argparse
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "ed"))
from cpqmc import CPMC
from altermagnet_ed import build_hopping


def scan_point(lx, ly, nup, ndn, t1, t2, t3, t4, U, uxy, v, dt, nw, seed,
               nequil, nblocks, bp):
    K = build_hopping(lx, ly, t1, t2, t3, t4)            # two-orbital, n=2*lx*ly
    qmc = CPMC(lx, ly, nup, ndn, U=U, dt=dt, nwalkers=nw, seed=seed,
               K=K, uxy=uxy, v=v)
    r = qmc.run_bp_chid(nequil=nequil, nblocks=nblocks, bp=bp)
    e, eerr = qmc.run_bp(nequil=40, nblocks=max(nblocks // 2, 8), bp=bp)
    r["E"] = e; r["Eerr"] = eerr
    return r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lx", type=int, default=4); ap.add_argument("--ly", type=int, default=4)
    ap.add_argument("--nup", type=int, default=8); ap.add_argument("--ndn", type=int, default=8)
    ap.add_argument("--t1", type=float, default=-1.0); ap.add_argument("--t3", type=float, default=-1.0)
    ap.add_argument("--t4", type=float, default=-1.0)
    ap.add_argument("--U", type=float, default=1.0); ap.add_argument("--uxy", type=float, default=1.0)
    ap.add_argument("--v", type=float, default=0.1)
    ap.add_argument("--dt", type=float, default=0.02); ap.add_argument("--nw", type=int, default=480)
    ap.add_argument("--nequil", type=int, default=180); ap.add_argument("--nblocks", type=int, default=32)
    ap.add_argument("--bp", type=int, default=24); ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--t2list", type=float, nargs="+",
                    default=[-1.0, -0.8, -0.6, -0.4, -0.2])
    ap.add_argument("-o", "--out", default="")
    a = ap.parse_args()
    nsites = 2 * a.lx * a.ly
    fields = ["t2", "alpha", "nsites", "nup", "ndn", "U", "uxy", "v", "t3", "t4",
              "E", "Eerr", "chi_d", "chi_d_err", "chi_d_vertex", "chi_d_vertex_err",
              "chi_s", "chi_s_err", "chi_s_vertex", "chi_s_vertex_err",
              "Cd_tau0", "Cs_tau0"]
    rows = []
    print(f"# two-orbital CPQMC {a.lx}x{a.ly} nsites={nsites} nup={a.nup} ndn={a.ndn} "
          f"U={a.U} uxy={a.uxy} v={a.v} t1={a.t1} t3={a.t3} t4={a.t4} bp={a.bp} nw={a.nw}")
    print(f"# {'alpha':>6} {'E':>9} {'chi_d':>9} {'chi_dV':>9} {'chi_s':>9} {'chi_sV':>9} {'Cd(0)':>9} {'Cs(0)':>9}")
    for t2 in a.t2list:
        r = scan_point(a.lx, a.ly, a.nup, a.ndn, a.t1, t2, a.t3, a.t4, a.U, a.uxy, a.v,
                       a.dt, a.nw, a.seed, a.nequil, a.nblocks, a.bp)
        alpha = 1 - t2 / a.t1
        row = dict(t2=t2, alpha=alpha, nsites=nsites, nup=a.nup, ndn=a.ndn, U=a.U,
                   uxy=a.uxy, v=a.v, t3=a.t3, t4=a.t4,
                   E=r["E"], Eerr=r["Eerr"],
                   chi_d=r["chi_d"], chi_d_err=r["chi_d_err"],
                   chi_d_vertex=r["chi_d_vertex"], chi_d_vertex_err=r["chi_d_vertex_err"],
                   chi_s=r["chi_s"], chi_s_err=r["chi_s_err"],
                   chi_s_vertex=r["chi_s_vertex"], chi_s_vertex_err=r["chi_s_vertex_err"],
                   Cd_tau0=r["Cdtau"][0], Cs_tau0=r["Cstau"][0])
        rows.append(row)
        print(f"  {alpha:>6.2f} {r['E']:>9.3f} {r['chi_d']:>9.4f} {r['chi_d_vertex']:>9.4f} "
              f"{r['chi_s']:>9.4f} {r['chi_s_vertex']:>9.4f} {r['Cdtau'][0]:>9.4f} {r['Cstau'][0]:>9.4f}")
        sys.stdout.flush()
    if a.out:
        os.makedirs(os.path.dirname(a.out), exist_ok=True)
        with open(a.out, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields); w.writeheader()
            for row in rows:
                w.writerow(row)
        print(f"# wrote {a.out}")


if __name__ == "__main__":
    main()
