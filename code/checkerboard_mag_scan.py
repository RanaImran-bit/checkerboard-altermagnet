#!/usr/bin/env python3
"""MAGNETIC / altermagnetic-order-parameter scan over the (L, n, delta, U) cube.
Pairs 1:1 with checkerboard_fss_full.py (same L, n, delta, U, seeds) so the pairing
susceptibility can be re-plotted against the AM order parameter instead of against
U and delta separately -- HoKin's Delta_pol x-axis, and the PRL referee's
"is it induced by U / is it general" question.

Per point, from the EQUAL-TIME spin structure factor S^zz(q) = Ctau_q[0]:
  Spipi  = S^zz(pi,pi)                       Neel (s-wave) magnetic scale
  M      = sqrt(Spipi)                       staggered moment (the original proxy)
  Mdelta = M * delta                         the factorized Delta_pol proxy
  Psi_AM = (1/N) sum_q sin(qx) sin(qy) S(q)  GENERAL dxy-harmonic AM order parameter
           (zero for pure Neel and zero without magnetism, so it needs BOTH U and
            delta, but WITHOUT assuming the M*delta product form)
  chi_pipi = windowed chi_zz(pi,pi)          unequal-time magnetic susceptibility

Singlet subtlety: the CP ground state is a singlet, so the one-body polarization
Delta n(k) = n_up(k) - n_dn(k) is identically zero. Altermagnetism lives in the
two-point S^zz(q), so the order parameter must be a projection of the STRUCTURE
FACTOR, not a one-body spin splitting.

The FULL S(q) grid for EVERY point is saved to sq_L{L}_U{U}.npz (~200 kB per block),
so any other harmonic can be re-projected later with zero recompute, and the k-space
figures need no extra run.

  SIZES=8,10,12 US=0,2,4,6,8 NPROC=32 NSEED=6 nohup python checkerboard_mag_scan.py \
      </dev/null > mag.log 2>&1 &
Writes mag_L{L}_U{U}.csv (crash-safe) + mag_L{L}.csv + mag_all.csv + sq_L{L}_U{U}.npz.

DEPENDENCIES (all must sit in the same directory):
  cpqmc.py  checkerboard.py  unified_scan.py
(run_bp_chi_spin imports _shift_index from unified_scan -- easy to forget when copying.)
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
SIZES = [int(x) for x in os.environ.get("SIZES", "6,8,10,12").split(",")]
US = [float(x) for x in os.environ.get("US", "0,2,4,6,8").split(",")]
DELTAS = [0.0, 0.1, 0.2, 0.3, 0.4]
N_TARGETS = [2 * k / 36 for k in [9, 10, 12, 14, 16, 18]]   # n = 0.5 .. 1.0, matches fss_full


def psi_harmonic(Sq, kind="dxy"):
    """(1/N) sum_q g(q) S(q) on the (L,L) q-grid, q = 2*pi*m/L."""
    L = Sq.shape[0]
    q = 2 * np.pi * np.arange(L) / L
    qx = q[:, None]; qy = q[None, :]
    g = {"dxy":   np.sin(qx) * np.sin(qy),
         "dx2y2": np.cos(qx) - np.cos(qy),
         "s":     np.ones_like(qx * qy)}[kind]
    return float(np.sum(g * Sq) / (L * L))


def run_point(args):
    L, nup, delta, U, seed = args
    K = cb.checkerboard_hopping(L, L, T0, T1, -delta)
    q = CPMC(L, L, nup, nup, U=U, dt=DT, nwalkers=NW, seed=seed, K=K, K_dn=None)
    r = q.run_bp_chi_spin(nequil=NEQ, nblocks=NBLK, bp=BP)
    Sq = np.array(r["Ctau_q"][0])                      # equal-time S^zz(q), shape (L,L)
    Spipi = float(Sq[L // 2, L // 2])
    N = L * L
    # S(pi,pi) is EXTENSIVE (~N m^2) once there is long-range order, so the
    # size-comparable order parameter is m = sqrt(S/N).  M_raw = sqrt(S) keeps the
    # old pol_grid_all.csv convention so the L=6 run can be cross-checked against it.
    m = float(np.sqrt(max(Spipi, 0.0) / N))
    M_raw = float(np.sqrt(max(Spipi, 0.0)))
    row = (L, nup, delta, U, seed, Spipi, m, M_raw, m * delta,
           psi_harmonic(Sq, "dxy"), psi_harmonic(Sq, "dx2y2"),
           float(r["chi_pipi"]), float(r["chi_pipi_err"]), 2 * nup / (L * L))
    return row, Sq


COLS = ["L", "nup", "delta", "U", "seed", "Spipi", "m", "M_raw", "mdelta",
        "Psi_dxy", "Psi_dx2y2", "chi_pipi", "chi_pipi_err", "n"]

if __name__ == "__main__":
    print(f"magnetic scan: sizes={SIZES} US={US} NSEED={NSEED} on {NPROC} cores", flush=True)
    all_frames = []
    for L in SIZES:
        nups = sorted({max(1, round(nt * L * L / 2)) for nt in N_TARGETS})
        L_frames = []
        for U in US:
            jobs = [(L, nup, d, U, s) for nup in nups for d in DELTAS for s in SEEDS]
            print(f"[L={L} U={U:g}] {len(jobs)} points, nups={nups} ...", flush=True)
            t0 = time.time()
            with Pool(NPROC) as pool:
                out = pool.map(run_point, jobs)
            rows = [o[0] for o in out]
            df = pd.DataFrame(rows, columns=COLS)
            df.to_csv(f"mag_L{L}_U{U:g}.csv", index=False)          # per-(L,U) crash safety
            # full S(q) for every point -> any harmonic can be re-projected later
            np.savez_compressed(f"sq_L{L}_U{U:g}.npz",
                                Sq=np.array([o[1] for o in out]),
                                meta=np.array(rows, dtype=object), cols=np.array(COLS))
            L_frames.append(df); all_frames.append(df)
            print(f"[L={L} U={U:g}] saved  rows: {len(df)}  ({(time.time()-t0)/60:.1f} min)", flush=True)
        pd.concat(L_frames, ignore_index=True).to_csv(f"mag_L{L}.csv", index=False)
        print(f"[L={L}] saved mag_L{L}.csv", flush=True)
    pd.concat(all_frames, ignore_index=True).to_csv("mag_all.csv", index=False)
    print("saved mag_all.csv  total rows:", sum(len(f) for f in all_frames), flush=True)
