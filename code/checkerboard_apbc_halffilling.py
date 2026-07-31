#!/usr/bin/env python3
"""Half-filling (n=1) scan under ANTI-PERIODIC boundary conditions (APBC-x), which open a
gap at half-filling and give a TRIAL-CLEAN (non-degenerate) free-electron trial -- removing
the open-shell caveat that plagued periodic half-filling (meeting-point #11 follow-up).

Per (U, delta, seed) at half-filling it computes the pairing susceptibilities chi_d, chi_dxy
AND the staggered moment M = sqrt(S^z(pi,pi)). Run it BOTH ways for a direct comparison:
    NPROC=32 NSEED=3 APX=-1 python checkerboard_apbc_halffilling.py    # APBC-x (trial-clean)
    NPROC=32 NSEED=3 APX=1  python checkerboard_apbc_halffilling.py    # periodic (open-shell)
Writes apbc_half_APX{-1,1}.csv (cols: apx,U,delta,seed,chi_d,chi_dxy,M,n).
REQUIRES the updated checkerboard.py (checkerboard_hopping now takes apx, apy)."""
import os, sys
os.environ["OMP_NUM_THREADS"] = "1"; os.environ["MKL_NUM_THREADS"] = "1"
sys.path.insert(0, os.path.expanduser("~/qmc-platform-master/pyqmc"))
import numpy as np, pandas as pd
from multiprocessing import Pool
from cpqmc import CPMC
import checkerboard as cb

L, T0, T1 = 6, -1.0, 0.3
NW, NEQ, NBLK, BP, DT = 160, 60, 40, 16, 0.05
NPROC = int(os.environ.get("NPROC", min(os.cpu_count(), 32)))
APX   = int(os.environ.get("APX", -1))          # -1 = APBC-x (trial-clean), +1 = periodic
US     = [0.0, 2.0, 4.0, 6.0, 8.0]
DELTAS = [0.0, 0.1, 0.2, 0.3, 0.4]
NSEED  = int(os.environ.get("NSEED", 3)); SEEDS = list(range(1, NSEED + 1))
NUP = L * L // 2                                 # half-filling


def run_pt(args):
    U, delta, seed = args
    K = cb.checkerboard_hopping(L, L, T0, T1, -delta, apx=APX, apy=1)
    Fd = cb.nn_bond_factors(L, L)[1]; Fdxy = cb.diag_bond_factors(L, L)
    q = CPMC(L, L, NUP, NUP, U=U, dt=DT, nwalkers=NW, seed=seed, K=K, K_dn=None)
    r = cb.run_bp_chid_cb(q, {"d": Fd, "dxy": Fdxy}, nequil=NEQ, nblocks=NBLK, bp=BP)
    q2 = CPMC(L, L, NUP, NUP, U=U, dt=DT, nwalkers=NW, seed=seed, K=K, K_dn=None)
    rs = q2.run_bp_chi_spin(nequil=NEQ, nblocks=NBLK, bp=BP)
    M = float(np.sqrt(max(rs["Ctau_q"][0][L // 2][L // 2], 0.0)))
    return (APX, U, delta, seed, r["chi_d_vertex"], r["chi_dxy_vertex"], M)


if __name__ == "__main__":
    jobs = [(U, d, s) for U in US for d in DELTAS for s in SEEDS]
    print(f"APBC half-filling scan: APX={APX}  {len(jobs)} points on {NPROC} cores", flush=True)
    with Pool(NPROC) as pool:
        rows = pool.map(run_pt, jobs)
    df = pd.DataFrame(rows, columns=["apx", "U", "delta", "seed", "chi_d", "chi_dxy", "M"])
    df["n"] = 1.0
    out = f"apbc_half_APX{APX}.csv"; df.to_csv(out, index=False)
    print("saved", out, " rows:", len(df), flush=True)
