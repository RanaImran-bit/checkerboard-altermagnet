#!/usr/bin/env python3
"""CPMC (constrained-path, the manuscript's method) tam-scan of the single-band
spin-dependent altermagnet on a SQUARE cluster at (near) half filling: equal-time
FULL + CONNECTED/VERTEX d-wave & s-wave pairing and the AFM S(pi,pi), vs anisotropy
tam. One (tam, seed) per invocation for trivial parallelism (GNU parallel).

    python pyqmc/cpmc_tam_scan.py --lx 4 --ly 4 --nup 8 --ndn 8 --U 6 --tam 0.3 --seed 1
"""
from __future__ import annotations
import argparse, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from cpqmc import CPMC, am_hopping


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lx", type=int, default=4); ap.add_argument("--ly", type=int, default=4)
    ap.add_argument("--nup", type=int, default=8); ap.add_argument("--ndn", type=int, default=8)
    ap.add_argument("--t0", type=float, default=1.0); ap.add_argument("--U", type=float, default=6.0)
    ap.add_argument("--tam", type=float, default=0.0); ap.add_argument("--dt", type=float, default=0.01)
    ap.add_argument("--nw", type=int, default=240); ap.add_argument("--nequil", type=int, default=250)
    ap.add_argument("--nblocks", type=int, default=40); ap.add_argument("--bp", type=int, default=12)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    Ku, Kd = am_hopping(a.lx, a.ly, a.t0, a.tam, a.t1)
    q = CPMC(a.lx, a.ly, a.nup, a.ndn, U=a.U, dt=a.dt, nwalkers=a.nw, seed=a.seed, K=Ku, K_dn=Kd)
    r = q.run_bp_pairmag(nequil=a.nequil, nblocks=a.nblocks, bp=a.bp)
    # one CSV line: tam seed d_full d_full_err d_vertex d_vertex_err Sq Sq_err s_vertex
    def ve(k): return r[k]["value"], r[k]["error"]
    df, dfe = ve("pair_dwave"); dv, dve = ve("pair_dwave_vertex")
    sq, sqe = ve("Sq_pipi"); sv, sve = ve("pair_swave_vertex")
    print(f"{a.tam:.3f},{a.seed},{df:.6f},{dfe:.6f},{dv:.6f},{dve:.6f},"
          f"{sq:.6f},{sqe:.6f},{sv:.6f},{sve:.6f}", flush=True)


if __name__ == "__main__":
    main()
