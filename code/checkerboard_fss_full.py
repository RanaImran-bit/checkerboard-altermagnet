#!/usr/bin/env python3
"""FULL finite-size scan (#12) -- reproduce the L=6 phase-diagram scan at L=8,10,12.
Full filling range n=0.5..1.0, full U=0..8, all delta, 4 channels (on-site s, ext-s, dx2-y2, dxy),
6 seeds. Directly comparable to chi_grid_all.csv (L=6). Answers: does the whole (n,delta,U) structure
-- interaction-driven onset, channel swap -- survive to larger lattices.

Robust for a multi-day unattended run:
  - sizes SMALLEST-FIRST (8 -> 10 -> 12), U=0 first (fast, no aux fields),
  - writes fss_full_L{L}_U{U}.csv after EACH (L,U) block -> fine-grained crash safety,
  - writes fss_full_L{L}.csv per size and fss_full_all.csv at the end.

  NPROC=32 NSEED=6 nohup python checkerboard_fss_full.py > fss_full.log 2>&1 &
Env: SIZES="8,10,12"  US="0,2,4,6,8"  NSEED=6.  Self-contained (cpqmc.py + checkerboard.py local).

Cost (single-core/pt ~N^1.5: L8~15 L10~30 L12~55 min; 6 fillings x5 delta x6 seeds=180/L,U on 32 cores):
  L=8 ~7h, L=10 ~14h, L=12 ~25h  -> ~2 days total. U=0 blocks are near-free.
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
NSEED = int(os.environ.get("NSEED", 6)); SEEDS = list(range(1, NSEED + 1))
SIZES = [int(x) for x in os.environ.get("SIZES", "8,10,12").split(",")]
US = [float(x) for x in os.environ.get("US", "0,2,4,6,8").split(",")]
DELTAS = [float(x) for x in os.environ.get("DELTAS", "0,0.1,0.2,0.3,0.4").split(",")]
N_TARGETS = [2 * k / 36 for k in [9, 10, 12, 14, 16, 18]]   # n = 0.5, 0.556, 0.667, 0.778, 0.889, 1.0


def run(args):
    L, nup, delta, U, seed = args
    K = cb.checkerboard_hopping(L, L, T0, T1, -delta)
    Fs, Fd = cb.nn_bond_factors(L, L); Fdxy = cb.diag_bond_factors(L, L); Fon = np.eye(L * L)
    q = CPMC(L, L, nup, nup, U=U, dt=DT, nwalkers=NW, seed=seed, K=K, K_dn=None)
    r = cb.run_bp_chid_cb(q, {"son": Fon, "sext": Fs, "d": Fd, "dxy": Fdxy},
                          nequil=NEQ, nblocks=NBLK, bp=BP)
    return (L, nup, delta, U, seed, r["chi_son_vertex"], r["chi_sext_vertex"],
            r["chi_d_vertex"], r["chi_dxy_vertex"])


if __name__ == "__main__":
    print(f"FULL FSS: sizes={SIZES} US={US} NSEED={NSEED} on {NPROC} cores", flush=True)
    all_frames = []
    for L in SIZES:
        nups = sorted({max(1, round(nt * L * L / 2)) for nt in N_TARGETS})
        L_frames = []
        for U in US:
            jobs = [(L, nup, d, U, s) for nup in nups for d in DELTAS for s in SEEDS]
            print(f"[L={L} U={U:g}] {len(jobs)} points, nups={nups} ...", flush=True)
            t0 = time.time()
            with Pool(NPROC) as pool:
                rows = pool.map(run, jobs)
            df = pd.DataFrame(rows, columns=["L", "nup", "delta", "U", "seed",
                                             "chi_son", "chi_sext", "chi_d", "chi_dxy"])
            df["n"] = 2 * df["nup"] / (L * L)
            df.to_csv(f"fss_full_L{L}_U{U:g}.csv", index=False)      # per-(L,U) crash safety
            L_frames.append(df); all_frames.append(df)
            print(f"[L={L} U={U:g}] saved  rows: {len(df)}  ({(time.time()-t0)/60:.1f} min)", flush=True)
        pd.concat(L_frames, ignore_index=True).to_csv(f"fss_full_L{L}.csv", index=False)
        print(f"[L={L}] saved fss_full_L{L}.csv", flush=True)
    pd.concat(all_frames, ignore_index=True).to_csv("fss_full_all.csv", index=False)
    print("saved fss_full_all.csv  total rows:", sum(len(f) for f in all_frames), flush=True)
