#!/usr/bin/env python3
"""Combined finite-size scan for L = 8, 10, 12 (meeting point #12) -- OVERNIGHT run.
Full (n, delta) grid at each L (same target fillings as the L=6 grid), fixed U, periodic BC,
connected-vertex pairing susceptibilities chi_d, chi_dxy. Tests whether the channel swap
(dx2-y2 attractive in the doped window, dxy suppressed) survives to larger lattices.

Robust for unattended running:
  - sizes ordered SMALLEST-FIRST (L=8 -> 10 -> 12), so the cheap ones finish first,
  - writes fss_L{L}.csv the moment EACH size completes (partial progress is never lost),
  - writes fss_all.csv at the very end.

Self-contained: needs only cpqmc.py + checkerboard.py in the same directory.
  NPROC=32 NSEED=2 nohup python checkerboard_fss_combined.py > fss_combined.log 2>&1 &

Cost (single-core per point from the timing test, ~N^1.5): L8 ~15, L10 ~30, L12 ~55 min.
Full 8-filling grid, NSEED=2 (~80 pts/L on 32 cores): L8 ~45m, L10 ~90m, L12 ~165m -> ~5 h total.
n=1.0 is open-shell under periodic BC (as at L=6); treat that column as the periodic reference.
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
NSEED = int(os.environ.get("NSEED", 2)); SEEDS = list(range(1, NSEED + 1))
SIZES = [int(x) for x in os.environ.get("SIZES", "8,10,12").split(",")]
DELTAS = [0.0, 0.1, 0.2, 0.3, 0.4]
N_TARGETS = [2 * k / 36 for k in [4, 6, 8, 10, 12, 14, 16, 18]]   # the L=6 grid's fillings


def run(args):
    L, nup, delta, seed = args
    K = cb.checkerboard_hopping(L, L, T0, T1, -delta)
    Fd = cb.nn_bond_factors(L, L)[1]; Fdxy = cb.diag_bond_factors(L, L)
    q = CPMC(L, L, nup, nup, U=U, dt=DT, nwalkers=NW, seed=seed, K=K, K_dn=None)
    r = cb.run_bp_chid_cb(q, {"d": Fd, "dxy": Fdxy}, nequil=NEQ, nblocks=NBLK, bp=BP)
    return (L, nup, delta, seed, r["chi_d_vertex"], r["chi_dxy_vertex"])


if __name__ == "__main__":
    print(f"combined FSS: sizes={SIZES} U={U} NSEED={NSEED} on {NPROC} cores", flush=True)
    frames = []
    for L in SIZES:
        nups = sorted({max(1, round(nt * L * L / 2)) for nt in N_TARGETS})
        jobs = [(L, nup, d, s) for nup in nups for d in DELTAS for s in SEEDS]
        print(f"[L={L}] {len(jobs)} points, nups={nups} ...", flush=True)
        t0 = time.time()
        with Pool(NPROC) as pool:
            rows = pool.map(run, jobs)
        df = pd.DataFrame(rows, columns=["L", "nup", "delta", "seed", "chi_d", "chi_dxy"])
        df["n"] = 2 * df["nup"] / (L * L)
        df.to_csv(f"fss_L{L}.csv", index=False); frames.append(df)
        print(f"[L={L}] saved fss_L{L}.csv  rows: {len(df)}  ({(time.time()-t0)/60:.1f} min)", flush=True)
    pd.concat(frames, ignore_index=True).to_csv("fss_all.csv", index=False)
    print("saved fss_all.csv  total rows:", sum(len(f) for f in frames), flush=True)
