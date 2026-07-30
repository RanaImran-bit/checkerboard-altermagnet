#!/usr/bin/env python3
"""Finite-size TIMING test for the dynamic pairing susceptibility.
Runs ONE point (fixed n, delta, U, single seed) at L = 6, 8, 10 and reports the
wall-clock per point, so we can extrapolate the cost to 12x12 / 14x14 and decide
the largest size feasible for the finite-size-scaling panel (meeting point #12).

    python checkerboard_timing.py
    python checkerboard_timing.py 6 8 10 12      # custom size list

NOTE: this is only a cost probe -- the filling is approximate (not closed-shell).
For the real scaling run we will pick closed-shell fillings per L.
"""
import os, sys, time
os.environ["OMP_NUM_THREADS"] = "1"; os.environ["MKL_NUM_THREADS"] = "1"
sys.path.insert(0, os.path.expanduser("~/qmc-platform-master/pyqmc"))
import numpy as np
from cpqmc import CPMC
import checkerboard as cb

# ---- same production knobs as the U-scan so the timing is representative ----
N_TARGET, DELTA, U = 0.78, 0.4, 4.0
T0, T1 = -1.0, 0.3
NW, NEQ, NBLK, BP, DT = 160, 60, 40, 16, 0.05
LS = [int(x) for x in sys.argv[1:]] or [6, 8, 10]
# ---------------------------------------------------------------------------

print(f"timing dynamic chi:  n~{N_TARGET}, delta={DELTA}, U={U}, 1 seed", flush=True)
print(f"{'L':>3} {'sites':>6} {'nup':>4} {'n':>6} {'minutes':>9}   chi_d / chi_dxy", flush=True)
rows = []
for L in LS:
    nup = round(N_TARGET * L * L / 2)
    K = cb.checkerboard_hopping(L, L, T0, T1, -DELTA)
    Fd = cb.nn_bond_factors(L, L)[1]; Fdxy = cb.diag_bond_factors(L, L)
    t0 = time.perf_counter()
    q = CPMC(L, L, nup, nup, U=U, dt=DT, nwalkers=NW, seed=1, K=K, K_dn=None)
    r = cb.run_bp_chid_cb(q, {"d": Fd, "dxy": Fdxy}, nequil=NEQ, nblocks=NBLK, bp=BP)
    mins = (time.perf_counter() - t0) / 60
    rows.append((L, mins))
    print(f"{L:>3} {L*L:>6} {nup:>4} {2*nup/(L*L):>6.3f} {mins:>9.2f}   "
          f"{r['chi_d_vertex']:+.3f} / {r['chi_dxy_vertex']:+.3f}", flush=True)

# crude N^3 extrapolation to larger L from the largest measured point
if rows:
    L0, m0 = rows[-1]
    print("\nextrapolation (assuming cost ~ N_sites^3):", flush=True)
    for L in [12, 14]:
        if L not in [r[0] for r in rows]:
            est = m0 * ((L * L) / (L0 * L0)) ** 3
            print(f"  L={L} ({L*L} sites):  ~{est:>7.1f} min/point  (~{est/60:.1f} h)", flush=True)
