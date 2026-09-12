#!/usr/bin/env python3
"""Grid scan of the connected-vertex pairing susceptibility over the (n, delta) plane,
optionally swept over U (the interaction-driven test).

Produces one CSV per U value:  chi_grid_U{U}.csv  (columns: U,nup,delta,seed,chi_son,chi_sext,chi_d,chi_dxy,n)
and a combined chi_grid_all.csv. Each is valid input to plot_cb_fig14_chi_phasediagram.py.
Seeds are independent statistical samples that get averaged per (n, delta) at plot time.

Interpretation of the U-scan: if AM/pairing are interaction-driven (not a band-structure
artifact of delta), the vertex susceptibilities should vanish at U=0 and grow with U.

Run on 251 (or any node with the platform):
    python checkerboard_chi_grid.py                     # US x 8 fillings x 5 delta x 3 seeds
    NPROC=48 python checkerboard_chi_grid.py            # more cores
    python checkerboard_chi_grid.py 4.0                 # single U=4 only (one CSV)

Wall time ~= (#US x 120 / NPROC) x per-point cost; ~23 min per U at 6x6 on 32 cores.
"""
import os, sys
os.environ["OMP_NUM_THREADS"] = "1"; os.environ["MKL_NUM_THREADS"] = "1"
sys.path.insert(0, os.path.expanduser("~/qmc-platform-master/pyqmc"))
import numpy as np, pandas as pd
from multiprocessing import Pool
from cpqmc import CPMC
import checkerboard as cb

# ---- EDIT HERE: params (match the notebook SETUP cell) ----
L, T0, T1 = 6, -1.0, 0.3                   # lattice L (6x6), t0=-1, t1=+0.3 (=-t')
NW, NEQ, NBLK, BP, DT = 160, 60, 40, 16, 0.05
NPROC = int(os.environ.get("NPROC", min(os.cpu_count(), 32)))

US     = [0.0, 2.0, 4.0, 6.0, 8.0]         # <-- the U-scan (U=0 is the interaction-off control)
# fillings: survey grid, or closed-shell preset per L (see docs/closed_shell_fillings.md).
# FILLINGS=survey (default) -> broad (n,delta) map; FILLINGS=closed -> non-degenerate trial only.
SURVEY_NUPS  = [4, 6, 8, 10, 12, 14, 16, 18]                # nup=ndn; n = 2*nup/(L*L)
CLOSED_SHELL = {6: [13], 8: [25], 10: [37], 14: [57, 61, 73]}  # gapped Fermi level at all delta
FILLINGS_MODE = os.environ.get("FILLINGS", "survey")
NUPS   = CLOSED_SHELL.get(L, SURVEY_NUPS) if FILLINGS_MODE == "closed" else SURVEY_NUPS
DELTAS = [0.0, 0.1, 0.2, 0.3, 0.4]
NSEED  = int(os.environ.get("NSEED", 3))   # seeds per point; quick=3, production=6+
SEEDS  = list(range(1, NSEED + 1))         # independent samples -> mean + error bar
# A single command-line arg overrides US with one value: `python ... 4.0`
if len(sys.argv) > 1:
    US = [float(sys.argv[1])]
# -----------------------------------------------------------


def run_pair(args):
    """One CP-AFQMC point -> connected-vertex susceptibility in 4 channels:
    on-site s (son), extended-s (sext), dx2-y2 (d), dxy."""
    U, nup, delta, seed = args
    K = cb.checkerboard_hopping(L, L, T0, T1, -delta)
    Fs, Fd = cb.nn_bond_factors(L, L)          # extended-s (NN, all +1), dx2-y2
    Fdxy = cb.diag_bond_factors(L, L)          # dxy (sin kx sin ky)
    Fon = np.eye(L * L)                         # on-site s (local pair)
    q = CPMC(L, L, nup, nup, U=U, dt=DT, nwalkers=NW, seed=seed, K=K, K_dn=None)
    r = cb.run_bp_chid_cb(q, {"son": Fon, "sext": Fs, "d": Fd, "dxy": Fdxy},
                          nequil=NEQ, nblocks=NBLK, bp=BP)
    return (U, nup, delta, seed,
            r["chi_son_vertex"], r["chi_sext_vertex"],
            r["chi_d_vertex"], r["chi_dxy_vertex"])


if __name__ == "__main__":
    npoints = len(NUPS) * len(DELTAS) * len(SEEDS)
    print(f"U-scan: US={US}  ->  {len(US)} x {npoints} = {len(US) * npoints} points "
          f"on {NPROC} cores", flush=True)
    all_frames = []
    for U in US:
        jobs = [(U, nup, d, s) for nup in NUPS for d in DELTAS for s in SEEDS]
        print(f"[U={U}] running {len(jobs)} points ...", flush=True)
        with Pool(NPROC) as pool:
            rows = pool.map(run_pair, jobs)
        df = pd.DataFrame(rows, columns=["U", "nup", "delta", "seed",
                                         "chi_son", "chi_sext", "chi_d", "chi_dxy"])
        df["n"] = 2 * df["nup"] / (L * L)
        out = f"chi_grid_U{U:g}.csv"
        df.to_csv(out, index=False)          # write per-U so partial progress survives
        all_frames.append(df)
        print(f"[U={U}] saved {out}  rows: {len(df)}", flush=True)
    pd.concat(all_frames, ignore_index=True).to_csv("chi_grid_all.csv", index=False)
    print("saved chi_grid_all.csv  total rows:", sum(len(f) for f in all_frames), flush=True)
