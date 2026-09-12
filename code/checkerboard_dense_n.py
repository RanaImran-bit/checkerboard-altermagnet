#!/usr/bin/env python3
"""DENSE-FILLING pairing scan -- publication-resolution (n, delta) phase diagram.

Matches how the group's published L=14 study sampled filling: step the electron number by one
(or two) rather than taking a handful of round fractions. At L=12 the available fillings are
n = nup/72 for nup = 36..72, i.e. 37 values with Delta n = 0.014 -- against the 6 values
(Delta n = 0.10) of the coarse cube, which cannot locate a phase boundary.

Fixed U (default 4, the working point); the U axis is already mapped at coarse filling by
checkerboard_fss_full.py. All four channels, connected vertex, 6 seeds.

  NUP_MIN=36 NUP_MAX=48 L=12 U=4 NPROC=32 NSEED=6 nohup python checkerboard_dense_n.py \
      </dev/null > dense.log 2>&1 &
Writes dense_L{L}_U{U}_nup{min}-{max}.csv after each filling (crash-safe) -- one row block per nup,
so a killed job loses at most one filling.

Cost at L=12: ~55 core-min per point; 5 delta x 6 seeds = 30 points per filling ~= 27 core-hours,
so ~52 min per filling on 32 cores. 13 fillings ~= 11 h on one node.
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
NSEED = int(os.environ.get("NSEED", 6)); SEEDS = list(range(1, NSEED + 1))
L = int(os.environ.get("L", 12))
U = float(os.environ.get("U", 4.0))
STEP = int(os.environ.get("STEP", 1))            # 1 = every electron, 2 = every other
NUP_MIN = int(os.environ.get("NUP_MIN", L * L // 4))        # n = 0.5
NUP_MAX = int(os.environ.get("NUP_MAX", L * L // 2))        # n = 1.0
DELTAS = [float(x) for x in os.environ.get("DELTAS", "0,0.1,0.2,0.3,0.4").split(",")]


def run(args):
    L_, nup, delta, U_, seed = args
    K = cb.checkerboard_hopping(L_, L_, T0, T1, -delta)
    Fs, Fd = cb.nn_bond_factors(L_, L_); Fdxy = cb.diag_bond_factors(L_, L_); Fon = np.eye(L_ * L_)
    q = CPMC(L_, L_, nup, nup, U=U_, dt=DT, nwalkers=NW, seed=seed, K=K, K_dn=None)
    r = cb.run_bp_chid_cb(q, {"son": Fon, "sext": Fs, "d": Fd, "dxy": Fdxy},
                          nequil=NEQ, nblocks=NBLK, bp=BP)
    return (L_, nup, delta, U_, seed, r["chi_son_vertex"], r["chi_sext_vertex"],
            r["chi_d_vertex"], r["chi_dxy_vertex"])


COLS = ["L", "nup", "delta", "U", "seed", "chi_son", "chi_sext", "chi_d", "chi_dxy"]

if __name__ == "__main__":
    nups = list(range(NUP_MIN, NUP_MAX + 1, STEP))
    out = f"dense_L{L}_U{U:g}_nup{NUP_MIN}-{NUP_MAX}.csv"
    print(f"dense-n scan: L={L} U={U:g} nup={NUP_MIN}..{NUP_MAX} step {STEP} "
          f"({len(nups)} fillings, n={2*NUP_MIN/(L*L):.3f}..{2*NUP_MAX/(L*L):.3f}) "
          f"x {len(DELTAS)} delta x {NSEED} seeds on {NPROC} cores -> {out}", flush=True)
    frames = []
    for i, nup in enumerate(nups):
        jobs = [(L, nup, d, U, s) for d in DELTAS for s in SEEDS]
        t0 = time.time()
        with Pool(NPROC) as pool:
            rows = pool.map(run, jobs)
        df = pd.DataFrame(rows, columns=COLS); df["n"] = 2 * df["nup"] / (L * L)
        frames.append(df)
        pd.concat(frames, ignore_index=True).to_csv(out, index=False)   # rewrite after each filling
        print(f"  [{i+1}/{len(nups)}] nup={nup} n={2*nup/(L*L):.4f} done "
              f"({(time.time()-t0)/60:.1f} min)", flush=True)
    print(f"saved {out}  total rows: {sum(len(f) for f in frames)}", flush=True)
