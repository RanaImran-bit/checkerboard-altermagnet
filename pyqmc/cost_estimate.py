#!/usr/bin/env python3
"""Per-step CPMC cost: fixed vs adaptive vs multidet(k) trials. The multidet
(CASSCF) propagation evaluates the full k-determinant overlap for every candidate
HS field, so its per-step cost is ~k x the single-determinant (fixed/adaptive)
cost. This measures the constant and confirms the ratio is set by k (the number
of CAS determinants), ~independent of lattice size -> the basis for the L=20x20
time estimate.

    python pyqmc/cost_estimate.py
"""
from __future__ import annotations
import os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(__file__)); sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "ed"))
from cpqmc import CPMC, square_hopping
from casci import transform_hubbard, casci_solve, build_multidet_trial


def time_steps(q, nstep=10):
    for _ in range(3):
        q.step()
    t0 = time.time()
    for _ in range(nstep):
        q.step()
    return (time.time() - t0) / nstep


def main():
    print(f"{'config':>16} {'trial':>10} {'k':>5} {'ms/step':>9} {'x fixed':>8}")
    for (lx, ly, nup, ndn) in [(4, 4, 2, 2), (4, 4, 4, 4), (6, 6, 2, 2)]:
        n = lx * ly
        tag = f"{lx}x{ly} {nup}+{ndn}"
        qf = CPMC(lx, ly, nup, ndn, t=1.0, U=4, dt=0.01, nwalkers=60, seed=1, trial="fixed")
        tf = time_steps(qf)
        print(f"{tag:>16} {'fixed':>10} {1:>5} {tf*1e3:>9.2f} {1.0:>8.1f}")
        qa = CPMC(lx, ly, nup, ndn, t=1.0, U=4, dt=0.01, nwalkers=60, seed=1, trial="adaptive", trial_every=10)
        ta = time_steps(qa)
        print(f"{'':>16} {'adaptive':>10} {1:>5} {ta*1e3:>9.2f} {ta/tf:>8.1f}")
        # k determinants for the TIMING only (quality irrelevant): random
        # orthonormal Slater matrices -- avoids the exponential full-CI build.
        rng = np.random.default_rng(0)
        def rand_det(ne):
            Q, _ = np.linalg.qr(rng.standard_normal((n, ne))); return Q
        for k in [4, 16, 36]:
            ups = [rand_det(nup) for _ in range(k)]
            dns = [rand_det(ndn) for _ in range(k)]
            coef = np.ones(k) / np.sqrt(k)
            qm = CPMC(lx, ly, nup, ndn, t=1.0, U=4, dt=0.01, nwalkers=60, seed=1, trial="multidet")
            qm.trial.set_multidet(ups, dns, coef); qm.trial._refresh_olp(qm.walkers)
            tm = time_steps(qm)
            print(f"{'':>16} {'multidet':>10} {k:>5} {tm*1e3:>9.2f} {tm/tf:>8.1f}")


if __name__ == "__main__":
    main()
