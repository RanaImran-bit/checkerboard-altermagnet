#!/usr/bin/env python3
"""Demonstrate the occupation-driven active-space selector + ED-free accuracy
estimate as correlation strength / doping varies. For each (U, filling) it runs a
short STABLE fixed-trial CPMC, forms the natural orbitals, and reports:
  entropy   : occupation entanglement entropy (multireference-ness)
  n_active  : auto-selected active size (grows where the state is multireference)
  dE_trunc  : |E_cas(n_active+1)-E_cas(n_active)| -- ED-free accuracy/ceiling proxy

This isolates the SELECTION + ACCURACY feature from the multi-determinant
propagation (which is separately ill-conditioned at very large k).

    python pyqmc/scan_active.py --lx 4 --ly 4 --nup 2 --ndn 2 --Us 1 2 4 8
"""
from __future__ import annotations
import argparse, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(__file__)); sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "ed"))
from cpqmc import CPMC, square_hopping
from casci import (ensemble_rdm, natural_orbitals, occupation_entropy,
                   auto_active_space, cas_accuracy_estimate)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lx", type=int, default=4); ap.add_argument("--ly", type=int, default=4)
    ap.add_argument("--nup", type=int, default=2); ap.add_argument("--ndn", type=int, default=2)
    ap.add_argument("--Us", type=float, nargs="+", default=[1, 2, 4, 8])
    ap.add_argument("--occthr", type=float, default=0.05); ap.add_argument("--maxact", type=int, default=8)
    ap.add_argument("--nw", type=int, default=200); ap.add_argument("--nequil", type=int, default=150)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    K = square_hopping(a.lx, a.ly, 1.0)
    print(f"# {a.lx}x{a.ly} nup={a.nup} ndn={a.ndn}  occthr={a.occthr} maxact={a.maxact}")
    print(f"{'U':>5} {'entropy':>8} {'n_active':>9} {'k(dets)':>8} {'dE_trunc':>9} {'occ (top)':>26}")
    for U in a.Us:
        q = CPMC(a.lx, a.ly, a.nup, a.ndn, t=1.0, U=U, dt=0.01, nwalkers=a.nw,
                 seed=a.seed, trial="fixed")
        q.run(nequil=a.nequil, nmeas=0)
        ru, rd = ensemble_rdm(q)
        occ, W = natural_orbitals(ru + rd)
        S = occupation_entropy(occ)
        nc, na = auto_active_space(occ, a.nup, a.ndn, thr=a.occthr, max_active=a.maxact)
        from math import comb
        nau, nad = a.nup - nc, a.ndn - nc
        k = comb(na, max(nau, 0)) * comb(na, max(nad, 0))
        _, dE, _ = cas_accuracy_estimate(K, U, W, a.nup, a.ndn, nc, na)
        print(f"{U:>5.1f} {S:>8.3f} {na:>9} {k:>8} {dE:>9.4f}   {np.round(occ[:6], 2)}")


if __name__ == "__main__":
    main()
