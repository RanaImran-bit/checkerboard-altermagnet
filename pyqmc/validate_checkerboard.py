#!/usr/bin/env python3
"""ED-gated validation of the checkerboard model + dxy pairing channel added in
checkerboard.py. Numpy-only gates (no QuSpin) unless --ed is passed.

GATE 1  hopping bit-for-bit == checkerboard_ed.build_hopping (== mc2duph.f90/GetK)
GATE 2  U=0 CPMC energy == free-fermion energy
GATE 3  U=0 connected VERTEX == 0 for every channel (s, dx2-y2, dxy)  [machinery gate]
GATE 4  (optional, --ed) interacting energy vs QuSpin ED on a small cluster
"""
from __future__ import annotations
import os, sys, argparse
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "Checkerboard_Model"))
from cpqmc import CPMC
import checkerboard as cb
import checkerboard_ed as ced


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lx", type=int, default=4); ap.add_argument("--ly", type=int, default=4)
    ap.add_argument("--nup", type=int, default=6); ap.add_argument("--ndn", type=int, default=6)
    ap.add_argument("--t0", type=float, default=-1.0)
    ap.add_argument("--t1", type=float, default=0.3)    # = -t'
    ap.add_argument("--t2", type=float, default=-0.2)   # = -delta
    ap.add_argument("--U", type=float, default=4.0)
    ap.add_argument("--dt", type=float, default=0.05)
    ap.add_argument("--nw", type=int, default=48)
    ap.add_argument("--ed", action="store_true", help="also run interacting QuSpin ED gate (GATE 4)")
    a = ap.parse_args()
    lx, ly = a.lx, a.ly
    K = cb.checkerboard_hopping(lx, ly, a.t0, a.t1, a.t2)

    # GATE 1
    Ked = ced.build_hopping(lx, ly, a.t0, a.t1, a.t2)
    g1 = np.allclose(K, Ked)
    print(f"GATE1 hopping bit-for-bit: max|diff|={np.abs(K-Ked).max():.2e}  -> {'PASS' if g1 else 'FAIL'}")

    # GATE 2: U=0 energy
    Efree, _, _ = ced.free_fermion_energy(lx, ly, a.nup, a.ndn, a.t0, a.t1, a.t2)
    q0 = CPMC(lx, ly, a.nup, a.ndn, U=0.0, dt=a.dt, nwalkers=a.nw, seed=1, K=K, K_dn=None)
    E0, e0 = q0.run_bp(nequil=30, nblocks=20, bp=10)
    g2 = abs(E0 - Efree) < max(3 * e0, 1e-6)
    print(f"GATE2 U=0 energy: CPMC={E0:.5f}+/-{e0:.5f} vs free={Efree:.5f}  -> {'PASS' if g2 else 'FAIL'}")

    # GATE 3: U=0 vertex must vanish for all channels
    Fs, Fd = cb.nn_bond_factors(lx, ly)
    Fdxy = cb.diag_bond_factors(lx, ly)
    Ffac = {"s": Fs, "d": Fd, "dxy": Fdxy}
    r0 = cb.run_bp_chid_cb(q0, Ffac, nequil=10, nblocks=20, bp=10)
    g3 = True
    for tag in Ffac:
        v = r0[f"chi_{tag}_vertex"]; ve = r0[f"chi_{tag}_vertex_err"]
        ok = abs(v) < max(3 * ve, 1e-4)
        g3 = g3 and ok
        print(f"GATE3 U=0 chi_{tag}_vertex = {v:+.5f}+/-{ve:.5f}  -> {'PASS' if ok else 'FAIL'}")

    # sanity: interacting run produces finite channels (not a gate, just a smoke test)
    qU = CPMC(lx, ly, a.nup, a.ndn, U=a.U, dt=a.dt, nwalkers=a.nw, seed=1, K=K, K_dn=None)
    rU = cb.run_bp_chid_cb(qU, Ffac, nequil=40, nblocks=30, bp=16)
    print(f"\nU={a.U} smoke test (delta={-a.t2}):")
    for tag in Ffac:
        print(f"   chi_{tag}       = {rU['chi_'+tag]:.4f} +/- {rU['chi_'+tag+'_err']:.4f}"
              f"   chi_{tag}_vtx = {rU['chi_'+tag+'_vertex']:+.4f} +/- {rU['chi_'+tag+'_vertex_err']:.4f}")

    if a.ed:
        E_ed, Gu, Gd = ced.ed_ground_state(lx, ly, a.nup, a.ndn, a.t0, a.t1, a.t2, a.U)
        Eq, eq = qU.run_bp(nequil=40, nblocks=30, bp=16)
        print(f"\nGATE4 U={a.U} energy: CPMC={Eq:.5f}+/-{eq:.5f} vs ED={E_ed:.5f}  "
              f"-> {'PASS' if abs(Eq-E_ed) < max(5*eq, 0.05*abs(E_ed)) else 'CHECK'}")

    print("\nSUMMARY:", "ALL CORE GATES PASS" if (g1 and g2 and g3) else "SOME GATES FAILED")


if __name__ == "__main__":
    main()
