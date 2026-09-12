#!/usr/bin/env python3
"""Scan the STAGGERED MOMENT M(U, delta, n) over the same (n, delta) x U grid as the
pairing scan, so the pairing susceptibility can be re-plotted against the altermagnetic
'polarization' Delta_pol = M * delta  (the polarization x-axis).

M = sqrt(S^z(pi,pi)), where S^z(pi,pi) is the EQUAL-TIME spin structure factor
(= Ctau_q[0] at (pi,pi) from run_bp_chi_spin). CP-clean 2-point observable.

NOTE on the singlet subtlety: the CP-AFQMC ground state is a singlet, so the single-particle
polarization Delta n(k) = n_up(k)-n_dn(k) is identically ZERO. The altermagnetism lives in
the SPIN CORRELATIONS, so the order-parameter scale is the staggered moment M (from S^z(pi,pi)),
and the dxy spin-splitting amplitude is m_AM ~ M * delta. That is the x-axis.

    python checkerboard_polarization_scan.py          # full grid
    NPROC=32 NSEED=6 python checkerboard_polarization_scan.py
Writes pol_grid_U{U}.csv + pol_grid_all.csv (cols: U,nup,delta,seed,Spipi,M,n).
Pair with chi_grid_all.csv on (U,nup,delta) to plot chi vs Delta_pol = M*delta.
"""
import os, sys
os.environ["OMP_NUM_THREADS"] = "1"; os.environ["MKL_NUM_THREADS"] = "1"
sys.path.insert(0, os.path.expanduser("~/qmc-platform-master/pyqmc"))
import numpy as np, pandas as pd
from multiprocessing import Pool
from cpqmc import CPMC
import checkerboard as cb

# ---- same knobs as checkerboard_chi_grid.py so the grids pair up ----
L, T0, T1 = 6, -1.0, 0.3
NW, NEQ, NBLK, BP, DT = 160, 60, 40, 16, 0.05
NPROC = int(os.environ.get("NPROC", min(os.cpu_count(), 32)))
US     = [0.0, 2.0, 4.0, 6.0, 8.0]
NUPS   = [4, 6, 8, 10, 12, 14, 16, 18]
DELTAS = [0.0, 0.1, 0.2, 0.3, 0.4]
NSEED  = int(os.environ.get("NSEED", 3))   # moment is CP-clean -> 3 usually enough
SEEDS  = list(range(1, NSEED + 1))
if len(sys.argv) > 1:
    US = [float(sys.argv[1])]
# --------------------------------------------------------------------


def run_pol(args):
    """One point -> (U, nup, delta, seed, S^z(pi,pi)_equal-time, M=sqrt(.))."""
    U, nup, delta, seed = args
    K = cb.checkerboard_hopping(L, L, T0, T1, -delta)
    q = CPMC(L, L, nup, nup, U=U, dt=DT, nwalkers=NW, seed=seed, K=K, K_dn=None)
    r = q.run_bp_chi_spin(nequil=NEQ, nblocks=NBLK, bp=BP)
    Spipi = float(r["Ctau_q"][0][L // 2][L // 2])   # tau=0 slice at (pi,pi)
    M = float(np.sqrt(max(Spipi, 0.0)))
    return (U, nup, delta, seed, Spipi, M)


if __name__ == "__main__":
    npoints = len(NUPS) * len(DELTAS) * len(SEEDS)
    print(f"polarization scan: US={US} -> {len(US)} x {npoints} = {len(US)*npoints} points "
          f"on {NPROC} cores", flush=True)
    frames = []
    for U in US:
        jobs = [(U, nup, d, s) for nup in NUPS for d in DELTAS for s in SEEDS]
        print(f"[U={U}] running {len(jobs)} points ...", flush=True)
        with Pool(NPROC) as pool:
            rows = pool.map(run_pol, jobs)
        df = pd.DataFrame(rows, columns=["U", "nup", "delta", "seed", "Spipi", "M"])
        df["n"] = 2 * df["nup"] / (L * L)
        out = f"pol_grid_U{U:g}.csv"
        df.to_csv(out, index=False); frames.append(df)
        print(f"[U={U}] saved {out}  rows: {len(df)}", flush=True)
    pd.concat(frames, ignore_index=True).to_csv("pol_grid_all.csv", index=False)
    print("saved pol_grid_all.csv  total rows:", sum(len(f) for f in frames), flush=True)
