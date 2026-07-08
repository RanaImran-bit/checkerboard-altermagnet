#!/usr/bin/env python3
"""One (lattice, tam, seed) point of the SINGLE-BAND altermagnet dynamic pairing
susceptibility, for embarrassingly-parallel cluster sweeps (xargs -P / one process
per core). Computes the unequal-time singlet d-wave and s-wave pairing
susceptibilities chi_a = integral C_a(tau) dtau (run_bp_chid) for the spin-
dependent anisotropic-hopping model (am_hopping; tam = altermagnet anisotropy).

chi_d is the d-wave PAIRING SUSCEPTIBILITY -- the tau-resolved/integrated
thermodynamic diagnostic the referee asked for (vs the equal-time vertex the
manuscript reported). C_a(tau=0) is the equal-time pair structure factor.

Prints one CSV line:  lx ly nup ndn U tam seed chi_d chi_d_err chi_s chi_s_err Cd0 Cs0

    python pyqmc/chid_tam_scan.py --lx 6 --ly 6 --nup 18 --ndn 18 --U 6 --tam 0.3 --seed 1
"""
from __future__ import annotations
import os
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("MKL_NUM_THREADS", "1")
import argparse, sys
sys.path.insert(0, os.path.dirname(__file__))
from cpqmc import CPMC, am_hopping


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lx", type=int, default=6); ap.add_argument("--ly", type=int, default=6)
    ap.add_argument("--nup", type=int, default=18); ap.add_argument("--ndn", type=int, default=18)
    ap.add_argument("--U", type=float, default=6.0); ap.add_argument("--t0", type=float, default=1.0)
    ap.add_argument("--tam", type=float, default=0.0)
    ap.add_argument("--t1", type=float, default=0.0,
                    help="NNN spin-dependent (d_xy) anisotropy — the nesting-BREAKING knob "
                         "that reproduces the manuscript trend (AFM down, d-wave up)")
    ap.add_argument("--dt", type=float, default=0.05); ap.add_argument("--nw", type=int, default=160)
    ap.add_argument("--nequil", type=int, default=60); ap.add_argument("--nblocks", type=int, default=40)
    ap.add_argument("--bp", type=int, default=16); ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    Ku, Kd = am_hopping(a.lx, a.ly, a.t0, a.tam, a.t1)
    q = CPMC(a.lx, a.ly, a.nup, a.ndn, U=a.U, dt=a.dt, nwalkers=a.nw, seed=a.seed, K=Ku, K_dn=Kd)
    r = q.run_bp_chid(nequil=a.nequil, nblocks=a.nblocks, bp=a.bp)
    # CSV: lx ly nup ndn U tam t1 seed  chi_d chi_d_err  chi_d_vtx chi_d_vtx_err
    #      chi_s chi_s_err  chi_s_vtx chi_s_vtx_err  Cd0 Cs0
    # ...Cd0 Cs0 = equal-time FULL pair structure factor; Cd0_vtx Cs0_vtx = equal-time
    # CONNECTED vertex (tau=0 slice of the connected correlator)
    print("%d %d %d %d %g %g %g %d %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f" % (
        a.lx, a.ly, a.nup, a.ndn, a.U, a.tam, a.t1, a.seed,
        r["chi_d"], r["chi_d_err"], r["chi_d_vertex"], r["chi_d_vertex_err"],
        r["chi_s"], r["chi_s_err"], r["chi_s_vertex"], r["chi_s_vertex_err"],
        r["Cdtau"][0], r["Cstau"][0], r["Cdtau_vertex"][0], r["Cstau_vertex"][0]), flush=True)


if __name__ == "__main__":
    main()
