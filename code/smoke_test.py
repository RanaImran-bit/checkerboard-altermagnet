#!/usr/bin/env python3
"""2-minute smoke test for a compute node before committing a long run.

Exercises exactly the code path the real driver uses -- same imports, same
multiprocessing.Pool fan-out at full width, same CPMC + checkerboard calls, same
CSV write -- but with tiny QMC settings so it finishes in ~1-2 min.

Run it the SAME way the real job is launched (detached, nohup, full python path)
so the test also covers ssh detachment and conda resolution:

  ssh -f 252 'cd ~/dense_run && NPROC=32 nohup /opt/anaconda3/bin/python smoke_test.py \
      </dev/null > smoke.log 2>&1 &'

PASS looks like:  "SMOKE TEST PASSED" plus a written smoke_out.csv.
Any SIGTERM/hang shows as a missing final line or a load-average near zero with
live processes (the wedge pattern).
"""
import os, sys, time, platform
os.environ["OMP_NUM_THREADS"] = "1"; os.environ["MKL_NUM_THREADS"] = "1"
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
from multiprocessing import Pool
from cpqmc import CPMC
import checkerboard as cb

NPROC = int(os.environ.get("NPROC", min(os.cpu_count(), 32)))
NJOBS = int(os.environ.get("NJOBS", 2 * NPROC))      # 2 full waves
L, U = 6, 4.0
T0, T1 = -1.0, 0.3


def run(args):
    seed, delta = args
    t0 = time.time()
    K = cb.checkerboard_hopping(L, L, T0, T1, -delta)
    Fs, Fd = cb.nn_bond_factors(L, L); Fdxy = cb.diag_bond_factors(L, L); Fon = np.eye(L * L)
    q = CPMC(L, L, 14, 14, U=U, dt=0.05, nwalkers=60, seed=seed, K=K, K_dn=None)
    r = cb.run_bp_chid_cb(q, {"son": Fon, "sext": Fs, "d": Fd, "dxy": Fdxy},
                          nequil=15, nblocks=8, bp=10)
    return (seed, delta, r["chi_d_vertex"], r["chi_dxy_vertex"], time.time() - t0)


if __name__ == "__main__":
    print(f"host={platform.node()}  python={sys.version.split()[0]}  numpy={np.__version__}", flush=True)
    print(f"os.cpu_count()={os.cpu_count()}  NPROC={NPROC}  NJOBS={NJOBS}", flush=True)
    print(f"loadavg before: {os.getloadavg()}", flush=True)
    t0 = time.time()
    jobs = [(s % 6 + 1, [0.0, 0.2, 0.4][s % 3]) for s in range(NJOBS)]
    with Pool(NPROC) as pool:
        rows = pool.map(run, jobs)
    wall = time.time() - t0
    df = pd.DataFrame(rows, columns=["seed", "delta", "chi_d", "chi_dxy", "sec"])
    df.to_csv("smoke_out.csv", index=False)

    cpu_sec = df.sec.sum()
    speedup = cpu_sec / wall
    print(f"loadavg after:  {os.getloadavg()}", flush=True)
    print(f"wall={wall:.1f}s  cpu={cpu_sec:.1f}s  effective parallelism={speedup:.1f}x "
          f"(ideal {NPROC})", flush=True)
    print(f"rows={len(df)}  finite chi_d={np.isfinite(df.chi_d).all()}  "
          f"chi_d mean={df.chi_d.mean():+.4f}  chi_dxy mean={df.chi_dxy.mean():+.4f}", flush=True)
    ok = (len(df) == NJOBS and np.isfinite(df[["chi_d", "chi_dxy"]].values).all()
          and speedup > 0.5 * NPROC)
    if not ok:
        print(f"SMOKE TEST FAILED  (rows {len(df)}/{NJOBS}, parallelism {speedup:.1f}/{NPROC})", flush=True)
        sys.exit(1)
    print("SMOKE TEST PASSED", flush=True)
