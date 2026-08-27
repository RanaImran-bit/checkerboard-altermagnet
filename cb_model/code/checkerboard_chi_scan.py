#!/usr/bin/env python3
"""One (checkerboard delta, seed) point of the DYNAMIC singlet pairing susceptibility
chi_a = integral C_a(tau) dtau (full + connected/VERTEX) for a in {s(ext-s), d(dx2-y2),
dxy}, via the ED-validated CP-AFQMC estimator. Embarrassingly parallel (one process per
core: xargs -P). Spin-independent checkerboard => AM is emergent.

Model matches mc2duph.f90/GetK: t0=-1 (NN), t1=-t'(=+0.3), t2=-delta.

Prints one CSV line (header with --header):
  lx ly nup ndn U delta seed  chi_s chi_s_err chi_s_vtx chi_s_vtx_err
  chi_d chi_d_err chi_d_vtx chi_d_vtx_err  chi_dxy chi_dxy_err chi_dxy_vtx chi_dxy_vtx_err

  python pyqmc/checkerboard_chi_scan.py --lx 6 --ly 6 --nup 11 --ndn 11 --U 4 --delta 0.2 --seed 1
"""
from __future__ import annotations
import os
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("MKL_NUM_THREADS", "1")
import argparse, sys
sys.path.insert(0, os.path.dirname(__file__))
from cpqmc import CPMC
import checkerboard as cb

HEADER = ("lx ly nup ndn U delta seed chi_s chi_s_err chi_s_vtx chi_s_vtx_err "
          "chi_d chi_d_err chi_d_vtx chi_d_vtx_err chi_dxy chi_dxy_err chi_dxy_vtx chi_dxy_vtx_err")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lx", type=int, default=6); ap.add_argument("--ly", type=int, default=6)
    ap.add_argument("--nup", type=int, default=11); ap.add_argument("--ndn", type=int, default=11)
    ap.add_argument("--U", type=float, default=4.0)
    ap.add_argument("--t0", type=float, default=-1.0)   # production NN (physical +1)
    ap.add_argument("--t1", type=float, default=0.3)    # = -t'
    ap.add_argument("--delta", type=float, default=0.2) # t2 = -delta
    ap.add_argument("--dt", type=float, default=0.05); ap.add_argument("--nw", type=int, default=160)
    ap.add_argument("--nequil", type=int, default=60); ap.add_argument("--nblocks", type=int, default=40)
    ap.add_argument("--bp", type=int, default=16); ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--header", action="store_true", help="print the CSV header and exit")
    a = ap.parse_args()
    if a.header:
        print(HEADER); return

    K = cb.checkerboard_hopping(a.lx, a.ly, a.t0, a.t1, -a.delta)   # spin-independent
    Fs, Fd = cb.nn_bond_factors(a.lx, a.ly)
    Fdxy = cb.diag_bond_factors(a.lx, a.ly)
    Ffac = {"s": Fs, "d": Fd, "dxy": Fdxy}

    q = CPMC(a.lx, a.ly, a.nup, a.ndn, U=a.U, dt=a.dt, nwalkers=a.nw, seed=a.seed, K=K, K_dn=None)
    r = cb.run_bp_chid_cb(q, Ffac, nequil=a.nequil, nblocks=a.nblocks, bp=a.bp)

    def g(tag, kind=""):   # kind = "" (full) or "_vertex"
        return r[f"chi_{tag}{kind}"], r[f"chi_{tag}{kind}_err"]
    vals = []
    for tag in ("s", "d", "dxy"):
        vals += list(g(tag)) + list(g(tag, "_vertex"))
    print("%d %d %d %d %g %g %d " % (a.lx, a.ly, a.nup, a.ndn, a.U, a.delta, a.seed)
          + " ".join(f"{v:.6f}" for v in vals))


if __name__ == "__main__":
    main()
