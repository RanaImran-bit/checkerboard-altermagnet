#!/usr/bin/env python3
"""FOCUSED finite-size scan (#12) -- publication quality, low noise.
Spends the whole compute budget on the ONE filling the finite-size figure uses (n~0.78, doped)
at HIGH seed count, instead of a broad grid at low seeds. delta-sweep, U=4, periodic BC, per L.

Cheap AND clean: 1 filling x 5 delta x 12 seeds = 60 points/L (2 waves on 32 cores).
  L=6 ~12m, L=8 ~30m, L=10 ~60m, L=12 ~110m  -> ~3.5 h for all four sizes.

  NPROC=32 NSEED=12 nohup python checkerboard_fss_focused.py > fss_focused.log 2>&1 &
Env: SIZES="6,8,10,12", NTARGET=0.78, NSEED, U. Writes fss_L{L}.csv (+ fss_all.csv).
Self-contained: needs cpqmc.py + checkerboard.py in the same directory.
"""
import os, sys, time
os.environ["OMP_NUM_THREADS"] = "1"; os.environ["MKL_NUM_THREADS"] = "1"
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
from multiprocessing import Pool
from cpqmc import CPMC
import checkerboard as cb

T0, T1 = -1.0, 0.3
NW, NEQ, NBLK, BP, DT = 160, 60, 40, 16, 0.05
NPROC = int(os.environ.get("NPROC", min(os.cpu_count(), 32)))
U = float(os.environ.get("U", 4.0))
NSEED = int(os.environ.get("NSEED", 12)); SEEDS = list(range(1, NSEED + 1))
SIZES = [int(x) for x in os.environ.get("SIZES", "6,8,10,12").split(",")]
NTARGET = float(os.environ.get("NTARGET", 0.78))
DELTAS = [0.0, 0.1, 0.2, 0.3, 0.4]


def run(args):
    L, nup, delta, seed = args
    K = cb.checkerboard_hopping(L, L, T0, T1, -delta)
    Fd = cb.nn_bond_factors(L, L)[1]; Fdxy = cb.diag_bond_factors(L, L)
    q = CPMC(L, L, nup, nup, U=U, dt=DT, nwalkers=NW, seed=seed, K=K, K_dn=None)
    r = cb.run_bp_chid_cb(q, {"d": Fd, "dxy": Fdxy}, nequil=NEQ, nblocks=NBLK, bp=BP)
    return (L, nup, delta, seed, r["chi_d_vertex"], r["chi_dxy_vertex"])


if __name__ == "__main__":
    print(f"focused FSS: sizes={SIZES} n~{NTARGET} U={U} NSEED={NSEED} on {NPROC} cores", flush=True)
    frames = []
    for L in SIZES:
        nup = round(NTARGET * L * L / 2)
        jobs = [(L, nup, d, s) for d in DELTAS for s in SEEDS]
        print(f"[L={L}] nup={nup} n={2*nup/(L*L):.3f}  {len(jobs)} points ...", flush=True)
        t0 = time.time()
        with Pool(NPROC) as pool:
            rows = pool.map(run, jobs)
        df = pd.DataFrame(rows, columns=["L", "nup", "delta", "seed", "chi_d", "chi_dxy"])
        df["n"] = 2 * df["nup"] / (L * L)
        df.to_csv(f"fss_L{L}.csv", index=False); frames.append(df)
        print(f"[L={L}] saved fss_L{L}.csv  rows: {len(df)}  ({(time.time()-t0)/60:.1f} min)", flush=True)
    pd.concat(frames, ignore_index=True).to_csv("fss_all.csv", index=False)
    print("saved fss_all.csv", flush=True)
