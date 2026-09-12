#!/usr/bin/env python3
"""FULL and VERTEX pair-field susceptibility, for the comparison with the
square-lattice benchmarks (White et al. PRB 39, 839 and Huang, Lin, Gubernatis
PRB 64, 205101).

checkerboard_fss_full.py keeps only the vertex, and stores it under column
names with no suffix, so `chi_d` in the existing masters is the VERTEX despite
the name. This driver keeps BOTH pieces and names them explicitly:

  chi_{a}_full     the full pair-field susceptibility, White's P
  chi_{a}_vertex   the connected part, full minus bubble
  bubble           = full - vertex, White's P_bar (derived at analysis time)

full > bubble (i.e. vertex > 0) means the interaction is ATTRACTIVE in that
channel, which is exactly the comparison White makes panel by panel. Reporting
the full susceptibility alone is misleading, because the bubble differs between
channels for purely kinematic reasons even at U = 0.

Default DELTAS=0: the benchmarks are square-lattice t-t' models, so they only
speak to the non-anisotropic limit of the checkerboard.

  SIZES=12 US=0,2,4 DELTAS=0 NPROC=32 NSEED=6 nohup python checkerboard_pairfull.py \
      </dev/null > pairfull.log 2>&1 &
Writes pairfull_L{L}_U{U}.csv after every delta, so a node reboot costs at most
one delta's work.
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
DELTAS = [float(x) for x in os.environ.get("DELTAS", "0").split(",")]
N_TARGETS = [2 * k / 36 for k in [9, 10, 12, 14, 16, 18]]   # n = 0.5 .. 1.0

CHAN = ["son", "sext", "d", "dxy"]
COLS = (["L", "nup", "delta", "U", "seed"]
        + [f"chi_{c}_{k}" for c in CHAN for k in ("full", "full_err",
                                                  "vertex", "vertex_err")]
        + ["n"])


def run(args):
    L, nup, delta, U, seed = args
    K = cb.checkerboard_hopping(L, L, T0, T1, -delta)
    Fs, Fd = cb.nn_bond_factors(L, L); Fdxy = cb.diag_bond_factors(L, L)
    Fon = np.eye(L * L)
    q = CPMC(L, L, nup, nup, U=U, dt=DT, nwalkers=NW, seed=seed, K=K, K_dn=None)
    r = cb.run_bp_chid_cb(q, {"son": Fon, "sext": Fs, "d": Fd, "dxy": Fdxy},
                          nequil=NEQ, nblocks=NBLK, bp=BP)
    vals = []
    for c in CHAN:
        vals += [r[f"chi_{c}"], r[f"chi_{c}_err"],
                 r[f"chi_{c}_vertex"], r[f"chi_{c}_vertex_err"]]
    return (L, nup, delta, U, seed, *vals, 2 * nup / (L * L))


if __name__ == "__main__":
    print(f"pair FULL+VERTEX: sizes={SIZES} US={US} deltas={DELTAS} "
          f"NSEED={NSEED} on {NPROC} cores", flush=True)
    for L in SIZES:
        nups = sorted({max(1, round(nt * L * L / 2)) for nt in N_TARGETS})
        for U in US:
            print(f"[L={L} U={U:g}] {len(nups)*len(DELTAS)*len(SEEDS)} points, "
                  f"nups={nups} ...", flush=True)
            t0 = time.time(); rows = []
            for d in DELTAS:
                jobs = [(L, nup, d, U, s) for nup in nups for s in SEEDS]
                td = time.time()
                with Pool(NPROC) as pool:
                    rows += pool.map(run, jobs)
                pd.DataFrame(rows, columns=COLS).to_csv(
                    f"pairfull_L{L}_U{U:g}.csv", index=False)
                print(f"   delta={d:g} done ({(time.time()-td)/60:.1f} min), "
                      f"{len(rows)} rows written", flush=True)
            print(f"[L={L} U={U:g}] saved  rows: {len(rows)}  "
                  f"({(time.time()-t0)/60:.1f} min)", flush=True)
    print("all blocks done", flush=True)
