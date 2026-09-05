#!/usr/bin/env python3
"""One pairing-susceptibility point per invocation, on one core.

The stock checkerboard_chi_grid.py fans a whole grid out over a multiprocessing
Pool inside a single job. For this campaign every parameter point gets its own
core through Slurm instead, so the work is one point per process and the array
scheduler does the fan-out.

It also drops the stock L=6: that size is barred for this project after it gave a
5.3 sigma wrong answer on the d_x2-y2 trend.

Writes one CSV row with EVERY quantity run_bp_chid_cb returns, not just the two
vertex channels the grid driver kept, so the full and connected pieces are both
on disk.

    chi_one.py <L> <U> <delta> <nup> <seed> <outcsv>

Sampling is set by env: NW, NEQ, NBLK, BP, DT.
"""
import os, sys
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
sys.path.insert(0, os.path.expanduser("~/qmc-platform-master/pyqmc"))
import numpy as np, pandas as pd
from cpqmc import CPMC
import checkerboard as cb

L     = int(sys.argv[1]);   U    = float(sys.argv[2])
delta = float(sys.argv[3]); nup  = int(sys.argv[4])
seed  = int(sys.argv[5]);   out  = sys.argv[6]

T0, T1 = -1.0, 0.3
NW   = int(os.environ.get("NW",   320))
NEQ  = int(os.environ.get("NEQ",  100))
NBLK = int(os.environ.get("NBLK",  60))
BP   = int(os.environ.get("BP",    20))
DT   = float(os.environ.get("DT", 0.05))

K    = cb.checkerboard_hopping(L, L, T0, T1, -delta)
Fd   = cb.nn_bond_factors(L, L)[1]
Fdxy = cb.diag_bond_factors(L, L)

q = CPMC(L, L, nup, nup, U=U, dt=DT, nwalkers=NW, seed=seed, K=K, K_dn=None)
r = cb.run_bp_chid_cb(q, {"d": Fd, "dxy": Fdxy}, nequil=NEQ, nblocks=NBLK, bp=BP)

row = {"L": L, "U": U, "delta": delta, "nup": nup, "ndn": nup,
       "n": 2 * nup / (L * L), "seed": seed,
       "NW": NW, "NEQ": NEQ, "NBLK": NBLK, "BP": BP, "DT": DT}
for k, v in r.items():
    row[k] = v if np.isscalar(v) else np.asarray(v).tolist()
pd.DataFrame([row]).to_csv(out, index=False)
print("wrote", out, "keys:", sorted(r.keys()), flush=True)
