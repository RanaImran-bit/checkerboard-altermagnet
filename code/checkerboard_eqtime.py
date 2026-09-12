"""Equal-time AND imaginary-time-integrated pairing, saved together.

run_bp_chid_cb already computes the whole curve C_a(tau) and returns it as
"C_{tag}_tau", but the production scan drivers wrote only the integrated chi and
threw the curve away. The equal-time value is simply C_a(tau=0).

Saving both lets us do two things at once:
  1. answer the equal-time vs unequal-time comparison,
  2. validate our pairing measurement against the group's Fortran CPQMC, which
     computes ONLY equal-time correlations.

Fortran channel correspondence, read off mc2duph.f90 lines 1353-1355 and 1386:
    sowave   = gp_up(m,n)*gp_dn(m,n)                      -> on-site s
    swave    = sf(i)*sf(j)*p_sdp,  sf  = (1,1,1,1)  NN    -> extended s
    dwave    = df(i)*df(j)*p_sdp,  df  = (1,1,-1,-1) NN   -> d_x2-y2
    dd12wave = dd1+dd2, ddf = (1,1,-1,-1) on DIAGONALS    -> d_xy
(sbwave/dbwave/pbwave use pdsfb, a bond-singlet four-fermion object, not ours.)

  L=10 NUPS=25 US=4 DELTAS=0,0.1,0.2,0.3,0.4 NSEED=6 NPROC=30 \
      python checkerboard_eqtime.py
"""
import os, sys, time, csv
os.environ["OMP_NUM_THREADS"] = "1"; os.environ["MKL_NUM_THREADS"] = "1"
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from multiprocessing import Pool
from cpqmc import CPMC
import checkerboard as cb

T0, T1 = -1.0, 0.3
# Hyperparameters, now settable from the environment so the CPQMC settings can
# be matched to the Fortran production runs (Delta_tau = 0.01, beta = 32,
# N_walkers = 1000, growth-control E_T = -50). NEQ is a STEP count, so the
# projection length before measurement is beta = NEQ * DT.
NW   = int(os.environ.get("NW", 160))
NEQ  = int(os.environ.get("NEQ", 60))        # beta = NEQ * DT
NBLK = int(os.environ.get("NBLK", 40))
BP   = int(os.environ.get("BP", 16))         # tau window for chi = BP * DT
DT   = float(os.environ.get("DT", 0.05))
NPROC = int(os.environ.get("NPROC", min(os.cpu_count(), 30)))
NSEED = int(os.environ.get("NSEED", 6))
L = int(os.environ.get("L", 10))
US = [float(x) for x in os.environ.get("US", "4").split(",")]
DELTAS = [float(x) for x in os.environ.get("DELTAS", "0,0.1,0.2,0.3,0.4").split(",")]
NUPS = [int(x) for x in os.environ.get("NUPS", "25").split(",")]
CHAN = ["son", "sext", "d", "dxy"]


def run_point(args):
    L_, nup, delta, U, seed = args
    n = L_ * L_
    K = cb.checkerboard_hopping(L_, L_, T0, T1, -delta)
    Fs, Fd = cb.nn_bond_factors(L_, L_); Fdxy = cb.diag_bond_factors(L_, L_)
    Ffac = {"son": np.eye(n), "sext": Fs, "d": Fd, "dxy": Fdxy}
    q = CPMC(L_, L_, nup, nup, U=U, dt=DT, nwalkers=NW, seed=seed, K=K, K_dn=None)
    r = cb.run_bp_chid_cb(q, Ffac, nequil=NEQ, nblocks=NBLK, bp=BP)
    row = dict(L=L_, nup=nup, delta=delta, U=U, seed=seed, n=2 * nup / n)
    for c in CHAN:
        # tau=0 is the EQUAL-TIME value; the rest of the curve integrates to chi
        row[f"eq_{c}_full"] = float(r[f"C_{c}_tau"][0])
        row[f"eq_{c}_vertex"] = float(r[f"C_{c}_tau_vertex"][0])
        row[f"chi_{c}_full"] = float(r[f"chi_{c}"])
        row[f"chi_{c}_vertex"] = float(r[f"chi_{c}_vertex"])
    # KEEP the tau-resolved curve. Every earlier driver computed C_a(tau) and
    # then stored only its integral, which is exactly why the truncation of
    # chi at tau_max could not be checked afterwards. Stored as a JSON-ish
    # semicolon list so the row stays a flat CSV record.
    row["taus"] = ";".join(f"{t:g}" for t in r["taus"])
    for c in CHAN:
        row[f"Ctau_{c}_vertex"] = ";".join(f"{v:.6g}" for v in r[f"C_{c}_tau_vertex"])
    row["beta_proj"] = NEQ * DT
    row["tau_max"] = BP * DT
    row["nw"] = NW
    return row


if __name__ == "__main__":
    jobs = [(L, nu, d, U, s) for nu in NUPS for d in DELTAS
            for U in US for s in range(1, NSEED + 1)]
    print(f"equal-time + integrated pairing: L={L} nups={NUPS} US={US} "
          f"deltas={DELTAS} x {NSEED} seeds = {len(jobs)} jobs on {NPROC} cores",
          flush=True)
    t0 = time.time()
    with Pool(NPROC) as p:
        rows = p.map(run_point, jobs)
    cols = list(rows[0].keys())
    out = f"eqtime_L{L}_U{'-'.join(f'{u:g}' for u in US)}.csv"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols); w.writeheader(); w.writerows(rows)
    print(f"wrote {out}  ({len(rows)} rows, {(time.time()-t0)/60:.1f} min)", flush=True)
