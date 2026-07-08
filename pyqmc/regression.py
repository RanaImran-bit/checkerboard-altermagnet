#!/usr/bin/env python3
"""Regression gate for the pyqmc refactor: runs a few deterministic (fixed-seed)
CPMC configurations and prints energies / a correlation checksum. Because the
seed is fixed and the algorithm unchanged, the numbers must be reproduced
bit-for-bit after the OOP refactor.

    python pyqmc/regression.py            # prints the canonical numbers
"""
from __future__ import annotations
import os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "ed"))
from cpqmc import CPMC
from altermagnet_ed import build_hopping


def main():
    out = []

    # 1) single-band Hubbard, mixed estimator
    q = CPMC(2, 2, 2, 2, t=1.0, U=4.0, dt=0.01, nwalkers=60, seed=1)
    e, err = q.run(nequil=40, nmeas=60, ortho=10, pc=10, meas_every=5)
    out.append(("hubbard_mixed", e, err))

    # 2) altermagnet uxx=2, mixed estimator
    K = build_hopping(2, 2, -1, -1, -1, -1)
    q = CPMC(2, 2, 4, 4, U=2.0, dt=0.01, nwalkers=60, seed=1, K=K)
    e, err = q.run(nequil=40, nmeas=60, ortho=10, pc=10, meas_every=5)
    out.append(("alt_uxx_mixed", e, err))

    # 3) altermagnet uxx=2, back-propagated energy
    q = CPMC(2, 2, 4, 4, U=2.0, dt=0.01, nwalkers=60, seed=1, K=K)
    e, err = q.run_bp(nequil=40, nblocks=12, bp=10, ortho=10, pc=10)
    out.append(("alt_uxx_bp", e, err))

    # 4) altermagnet uxy=2, back-propagated energy
    q = CPMC(2, 2, 4, 4, U=0.0, dt=0.01, nwalkers=60, seed=1, K=K, uxy=2.0)
    e, err = q.run_bp(nequil=40, nblocks=12, bp=10, ortho=10, pc=10)
    out.append(("alt_uxy_bp", e, err))

    # 5) back-propagated correlations checksum (sum of |green_up| + |szsz|)
    q = CPMC(2, 2, 4, 4, U=2.0, dt=0.01, nwalkers=60, seed=1, K=K)
    rec = q.run_bp_obs(nequil=40, nblocks=12, bp=10, ortho=10, pc=10)
    chk = float(np.abs(np.array(rec["green_up"])).sum()
                + np.abs(np.array(rec["szsz"])).sum())
    out.append(("corr_checksum", chk, 0.0))

    for name, v, e in out:
        print(f"{name:18s} {v:+.10f}  {e:+.10f}")


if __name__ == "__main__":
    main()
