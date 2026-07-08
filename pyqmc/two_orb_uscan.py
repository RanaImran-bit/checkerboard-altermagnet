#!/usr/bin/env python3
"""Two-orbital altermagnet CPQMC: scan on-site U at FIXED anisotropic bands
(leeb2024spontaneous / Fig.S9: t1,t2,t3,t4 fixed). Tracks whether the altermagnetic
order and d-wave pairing turn on together as the interaction grows.

Back-propagated equal-time observables (run_bp_pairmag):
  pair_dwave / pair_swave           : FULL pair structure factors S_d, S_s
  pair_dwave_vertex / _swave_vertex : CONNECTED (vertex) pair structure factors
  S_AM                              : altermagnetic order SF (orbital-staggered Sz)
  Sq_pipi, moment2                  : Neel SF and local-moment reference

  python pyqmc/two_orb_uscan.py --lx 4 --ly 4 --nup 16 --ndn 16 \
      --t1 1 --t2 1.75 --t3 0.85 --t4 0.65 --Ulist 0 2 4 8 -o out.csv
"""
from __future__ import annotations
import os, sys, csv, argparse
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")
import numpy as np
sys.path.insert(0, os.path.dirname(__file__)); sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "ed"))
from cpqmc import CPMC
from altermagnet_ed import build_hopping


def scan_point(lx, ly, nup, ndn, t1, t2, t3, t4, U, dt, nw, seed, nequil, nblocks, bp, uxy=0.0, v=0.0):
    K = build_hopping(lx, ly, t1, t2, t3, t4)
    qmc = CPMC(lx, ly, nup, ndn, U=U, dt=dt, nwalkers=nw, seed=seed, K=K, uxy=uxy, v=v)
    return qmc.run_bp_pairmag(nequil=nequil, nblocks=nblocks, bp=bp)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lx", type=int, default=4); ap.add_argument("--ly", type=int, default=4)
    ap.add_argument("--nup", type=int, default=16); ap.add_argument("--ndn", type=int, default=16)
    ap.add_argument("--t1", type=float, default=1.0); ap.add_argument("--t2", type=float, default=1.75)
    ap.add_argument("--t3", type=float, default=0.85); ap.add_argument("--t4", type=float, default=0.65)
    ap.add_argument("--uxy", type=float, default=0.0); ap.add_argument("--v", type=float, default=0.0)
    ap.add_argument("--dt", type=float, default=0.04); ap.add_argument("--nw", type=int, default=160)
    ap.add_argument("--nequil", type=int, default=120); ap.add_argument("--nblocks", type=int, default=24)
    ap.add_argument("--bp", type=int, default=16); ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--Ulist", type=float, nargs="+", default=[0,2,4,6,8])
    ap.add_argument("-o", "--out", default="")
    a = ap.parse_args()
    nsites = 2 * a.lx * a.ly
    fields = ["U","nsites","nup","ndn","t1","t2","t3","t4",
              "pair_dwave","pair_dwave_err","pair_dwave_vertex","pair_dwave_vertex_err",
              "pair_swave","pair_swave_vertex","S_AM","S_AM_err","Sq_pipi","moment2"]
    rows = []
    print(f"# CPQMC two-orbital U-scan {a.lx}x{a.ly} nsites={nsites} nup={a.nup} "
          f"t1={a.t1} t2={a.t2} t3={a.t3} t4={a.t4} (fixed bands; scanning U)")
    print(f"# {'U':>5} {'pairD_V':>9} {'pairS_V':>9} {'S_AM':>9} {'Sq_pipi':>9} {'mom2':>8}")
    for U in a.Ulist:
        r = scan_point(a.lx, a.ly, a.nup, a.ndn, a.t1, a.t2, a.t3, a.t4, U,
                       a.dt, a.nw, a.seed, a.nequil, a.nblocks, a.bp, uxy=a.uxy, v=a.v)
        g = lambda k: r[k]["value"]; ge = lambda k: r[k]["error"]
        rows.append(dict(U=U, nsites=nsites, nup=a.nup, ndn=a.ndn, t1=a.t1, t2=a.t2, t3=a.t3, t4=a.t4,
                         pair_dwave=g("pair_dwave"), pair_dwave_err=ge("pair_dwave"),
                         pair_dwave_vertex=g("pair_dwave_vertex"), pair_dwave_vertex_err=ge("pair_dwave_vertex"),
                         pair_swave=g("pair_swave"), pair_swave_vertex=g("pair_swave_vertex"),
                         S_AM=g("S_AM"), S_AM_err=ge("S_AM"), Sq_pipi=g("Sq_pipi"), moment2=g("moment2")))
        print(f"  {U:>5.1f} {g('pair_dwave_vertex'):>9.4f} {g('pair_swave_vertex'):>9.4f} "
              f"{g('S_AM'):>9.4f} {g('Sq_pipi'):>9.4f} {g('moment2'):>8.4f}")
        sys.stdout.flush()
    if a.out:
        os.makedirs(os.path.dirname(a.out), exist_ok=True)
        with open(a.out, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields); w.writeheader()
            for row in rows: w.writerow(row)
        print(f"# wrote {a.out}")


if __name__ == "__main__":
    main()
