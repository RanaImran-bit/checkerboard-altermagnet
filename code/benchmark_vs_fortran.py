"""Cross-code check: Python CPQMC against the group's validated Fortran CPQMC.

The Fortran code is the one used in the group's published work. The Python engine
was written to add the unequal-time (imaginary-time integrated) susceptibility,
which the Fortran does not compute. Everything the Python code has been checked
against so far is INTERNAL (U=0 null test, P(R) sum rule, twist null test) or an
ED benchmark of ENERGY and S_zz only. The pairing channels have never been
compared against an independent code.

Energy is the right cross-check: it exercises the hopping matrix, trial state,
propagation, constraint and back-propagation in one number.

Fortran reference, L10n0.500u4.0tA-0.2000tt0.3N25:
    total energy      -94.4962  +/- 0.0716   (E/N = -0.944962)
    magnetic moment    0.459808 +/- 0.000157
    dt = 0.01, back-propagation 40 steps

Parallel over (dt, seed); prints unbuffered so progress is visible.
"""
import os, sys, time
os.environ["OMP_NUM_THREADS"] = "1"; os.environ["MKL_NUM_THREADS"] = "1"
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from multiprocessing import Pool
from cpqmc import CPMC
import checkerboard as cb

L, NUP, U, DELTA = 10, 25, 4.0, 0.2
T0, T1 = -1.0, 0.3
REF_E, REF_E_ERR = -94.4962481357321, 0.0716285506711404
SEEDS = [1, 2, 3, 4, 5, 6]
DTS = [0.01, 0.05]           # 0.01 matches Fortran; 0.05 is our production step
K = cb.checkerboard_hopping(L, L, T0, T1, -DELTA)


def one(args):
    dt, seed = args
    q = CPMC(L, L, NUP, NUP, U=U, dt=dt, nwalkers=200, seed=seed, K=K, K_dn=None)
    e, _ = q.run_bp(nequil=150, nblocks=40, bp=15)
    return dt, seed, e


if __name__ == "__main__":
    print(f"L={L} nup=ndn={NUP} U={U} delta={DELTA} t={T0} t1={T1}", flush=True)
    print(f"Fortran: E = {REF_E:.4f} +/- {REF_E_ERR:.4f}  E/N = {REF_E/(L*L):.6f}\n",
          flush=True)
    jobs = [(dt, s) for dt in DTS for s in SEEDS]
    t0 = time.time()
    with Pool(len(jobs)) as p:
        out = p.map(one, jobs)
    print(f"({time.time()-t0:.0f}s)\n", flush=True)
    for dt in DTS:
        es = np.array([e for d, s, e in out if d == dt and np.isfinite(e)])
        m, sd = es.mean(), es.std(ddof=1) / np.sqrt(len(es))
        dev = (m - REF_E) / np.hypot(sd, REF_E_ERR)
        print(f"dt={dt:.2f}  Python E = {m:10.4f} +/- {sd:6.4f}   "
              f"E/N = {m/(L*L):.6f}   diff = {m-REF_E:+8.4f}  ({dev:+6.1f} sigma)",
              flush=True)
