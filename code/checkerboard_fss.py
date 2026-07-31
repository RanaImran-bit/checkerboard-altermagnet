#!/usr/bin/env python3
"""Finite-size scan (meeting point #12): does the channel swap survive larger L?
At a fixed doped filling (n ~ 0.78) and U, sweep delta and measure the connected-vertex
pairing susceptibilities chi_d, chi_dxy. Run once per lattice size L and compare.

SELF-CONTAINED: needs only cpqmc.py + checkerboard.py in the SAME directory (pairing only,
no run_bp_chi_spin / unified_scan). Import path is this script's own folder, so it runs on
any worker node once those 3 files are copied there.

    NPROC=32 python checkerboard_fss.py 8      # L=8
    NPROC=32 python checkerboard_fss.py 10     # L=10
    U=4 NTARGET=0.78 NSEED=3 python checkerboard_fss.py 12
Writes fss_L{L}.csv (cols: L,nup,delta,seed,chi_d,chi_dxy,n).
"""
import os, sys
os.environ["OMP_NUM_THREADS"] = "1"; os.environ["MKL_NUM_THREADS"] = "1"
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # local cpqmc + checkerboard
import numpy as np, pandas as pd
from multiprocessing import Pool
from cpqmc import CPMC
import checkerboard as cb

L = int(sys.argv[1]) if len(sys.argv) > 1 else 8
T0, T1 = -1.0, 0.3
NW, NEQ, NBLK, BP, DT = 160, 60, 40, 16, 0.05
NPROC = int(os.environ.get("NPROC", min(os.cpu_count(), 32)))
U = float(os.environ.get("U", 4.0))
N_TARGET = float(os.environ.get("NTARGET", 0.78))
DELTAS = [0.0, 0.1, 0.2, 0.3, 0.4]
NSEED = int(os.environ.get("NSEED", 3)); SEEDS = list(range(1, NSEED + 1))
NUP = round(N_TARGET * L * L / 2)                # doped filling nearest n_target


def run(args):
    delta, seed = args
    K = cb.checkerboard_hopping(L, L, T0, T1, -delta)
    Fd = cb.nn_bond_factors(L, L)[1]; Fdxy = cb.diag_bond_factors(L, L)
    q = CPMC(L, L, NUP, NUP, U=U, dt=DT, nwalkers=NW, seed=seed, K=K, K_dn=None)
    r = cb.run_bp_chid_cb(q, {"d": Fd, "dxy": Fdxy}, nequil=NEQ, nblocks=NBLK, bp=BP)
    return (L, NUP, delta, seed, r["chi_d_vertex"], r["chi_dxy_vertex"])


if __name__ == "__main__":
    jobs = [(d, s) for d in DELTAS for s in SEEDS]
    print(f"FSS L={L} nup={NUP} n={2*NUP/(L*L):.3f} U={U}: {len(jobs)} points on {NPROC} cores",
          flush=True)
    with Pool(NPROC) as pool:
        rows = pool.map(run, jobs)
    df = pd.DataFrame(rows, columns=["L", "nup", "delta", "seed", "chi_d", "chi_dxy"])
    df["n"] = 2 * df["nup"] / (L * L)
    df.to_csv(f"fss_L{L}.csv", index=False)
    print(f"saved fss_L{L}.csv  rows: {len(df)}", flush=True)
