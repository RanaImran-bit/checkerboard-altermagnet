#!/usr/bin/env python3
"""Twist-averaged pairing vertex: the fix for open-shell contamination.

Under periodic boundary conditions EVERY filling in our grid is open shell at
L = 12, at every anisotropy, so the trial state sits on a degenerate level. That
is visible in the data: at L = 8, n = 0.781 the extended-s vertex jumps sixteen
times its typical step exactly where the shell character flips with anisotropy.

A twist is a boundary phase, psi(r + L xhat) = exp(i theta_x) psi(r), which shifts
the k grid from 2 pi m / L to (2 pi m + theta) / L. Shell effects come from the
discrete grid landing awkwardly on the Fermi surface, and averaging over theta
fills the zone continuously and washes them out.

Only theta = 0 and pi keep the hopping matrix REAL, which matters because the
engine is real-only (cpqmc.py stores K as float and uses eigh and plain
transposes). That leaves four twists, of which two are exactly degenerate by the
x <-> y symmetry of the lattice (eigenvalues agree to 1e-14), so three distinct
runs suffice:

    (apx, apy) = (+1, +1)   weight 1   periodic, already computed
                 (-1, +1)   weight 2   half twist, stands for (+1, -1) too
                 (-1, -1)   weight 1   both antiperiodic

Average the observables with those weights. This is PBC/APBC averaging, a coarse
four-point quadrature over the twist zone. It removes the worst shell artefacts
but is not full twist averaging, which needs a continuous theta grid, hence
complex hoppings, hence a complex engine.

  APX=-1 APY=1 SIZES=12 US=4 NUPS=72 NPROC=32 NSEED=6 nohup python \
      checkerboard_twist.py </dev/null > twist.log 2>&1 &

NUPS selects fillings directly (comma separated); leave unset for the usual six.
Writes twist_L{L}_U{U}_apx{apx}apy{apy}_nups{...}.csv after every anisotropy. BOTH the
twist AND the fillings are in the filename: a run covering different fillings in the same
directory would otherwise overwrite the first, which is exactly what happened on 12 Aug.
"""
import os, sys, time, csv
os.environ["OMP_NUM_THREADS"] = "1"; os.environ["MKL_NUM_THREADS"] = "1"
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from multiprocessing import Pool
from cpqmc import CPMC
import checkerboard as cb

# Deliberately no pandas: node 257 has numpy and scipy but not pandas, and this
# driver only ever writes a flat table. The csv module is in the standard library
# everywhere, so the same file runs on every node without installing anything.


def write_csv(path, cols, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(cols); w.writerows(rows)

T0, T1 = -1.0, 0.3
NW, NEQ, NBLK, BP, DT = 160, 60, 40, 16, 0.05
NPROC = int(os.environ.get("NPROC", min(os.cpu_count(), 32)))
NSEED = int(os.environ.get("NSEED", 6)); SEEDS = list(range(1, NSEED + 1))
SIZES = [int(x) for x in os.environ.get("SIZES", "12").split(",")]
US = [float(x) for x in os.environ.get("US", "4").split(",")]
DELTAS = [float(x) for x in os.environ.get("DELTAS", "0,0.1,0.2,0.3,0.4,0.5,0.6,0.7").split(",")]
APX = int(os.environ.get("APX", 1)); APY = int(os.environ.get("APY", 1))
NUPS_ENV = os.environ.get("NUPS", "")
N_TARGETS = [2 * k / 36 for k in [9, 10, 12, 14, 16, 18]]

CHAN = ["son", "sext", "d", "dxy"]
COLS = (["L", "nup", "delta", "U", "seed", "apx", "apy"]
        + [f"chi_{c}_{k}" for c in CHAN for k in ("full", "vertex")] + ["n"])


def run(args):
    L, nup, delta, U, seed = args
    K = cb.checkerboard_hopping(L, L, T0, T1, -delta, apx=APX, apy=APY)
    Fs, Fd = cb.nn_bond_factors(L, L); Fdxy = cb.diag_bond_factors(L, L)
    Fon = np.eye(L * L)
    q = CPMC(L, L, nup, nup, U=U, dt=DT, nwalkers=NW, seed=seed, K=K, K_dn=None)
    r = cb.run_bp_chid_cb(q, {"son": Fon, "sext": Fs, "d": Fd, "dxy": Fdxy},
                          nequil=NEQ, nblocks=NBLK, bp=BP)
    vals = []
    for c in CHAN:
        vals += [r[f"chi_{c}"], r[f"chi_{c}_vertex"]]
    return (L, nup, delta, U, seed, APX, APY, *vals, 2 * nup / (L * L))


if __name__ == "__main__":
    assert APX in (-1, 1) and APY in (-1, 1), "APX/APY must be +1 or -1 to stay real"
    print(f"twist run: apx={APX:+d} apy={APY:+d} sizes={SIZES} US={US} "
          f"deltas={DELTAS} NSEED={NSEED} on {NPROC} cores", flush=True)
    for L in SIZES:
        if NUPS_ENV:
            nups = [int(x) for x in NUPS_ENV.split(",")]
        else:
            nups = sorted({max(1, round(nt * L * L / 2)) for nt in N_TARGETS})
        # Report which fillings this twist actually closes, so the run is
        # interpretable before any of the physics comes back.
        for U in US:
            K0 = cb.checkerboard_hopping(L, L, T0, T1, 0.0, apx=APX, apy=APY)
            gap = np.diff(np.sort(np.linalg.eigvalsh(K0)))
            shells = {n: ("closed" if gap[n - 1] > 1e-6 else "open") for n in nups}
            print(f"[L={L} U={U:g}] nups={nups}  shell at delta=0: {shells}", flush=True)
            t0 = time.time(); rows = []
            for d in DELTAS:
                jobs = [(L, nup, d, U, s) for nup in nups for s in SEEDS]
                td = time.time()
                with Pool(NPROC) as pool:
                    rows += pool.map(run, jobs)
                write_csv(f"twist_L{L}_U{U:g}_apx{APX}apy{APY}_nups{'-'.join(str(x) for x in nups)}.csv", COLS, rows)
                print(f"   delta={d:g} done ({(time.time()-td)/60:.1f} min), "
                      f"{len(rows)} rows", flush=True)
            print(f"[L={L} U={U:g}] saved {len(rows)} rows "
                  f"({(time.time()-t0)/60:.1f} min)", flush=True)
    print("all done", flush=True)
