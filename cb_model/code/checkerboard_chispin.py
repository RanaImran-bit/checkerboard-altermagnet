#!/usr/bin/env python3
"""Momentum-resolved MAGNETIC spin susceptibility chi_zz(q) for the checkerboard,
via the platform's run_bp_chi_spin (single-band, ED-validated, works for any K).
Saves the full BZ q-grid to .npy and prints the (pi,pi) / max values.

  python pyqmc/checkerboard_chispin.py --lx 6 --ly 6 --nup 11 --ndn 11 --U 4 --delta 0.4
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("MKL_NUM_THREADS", "1")
import argparse, sys
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from cpqmc import CPMC
import checkerboard as cb


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lx", type=int, default=6); ap.add_argument("--ly", type=int, default=6)
    ap.add_argument("--nup", type=int, default=11); ap.add_argument("--ndn", type=int, default=11)
    ap.add_argument("--U", type=float, default=4.0)
    ap.add_argument("--t0", type=float, default=-1.0); ap.add_argument("--t1", type=float, default=0.3)
    ap.add_argument("--delta", type=float, default=0.4)
    ap.add_argument("--dt", type=float, default=0.05); ap.add_argument("--nw", type=int, default=160)
    ap.add_argument("--nequil", type=int, default=60); ap.add_argument("--nblocks", type=int, default=40)
    ap.add_argument("--bp", type=int, default=16); ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--outdir", default=".")
    a = ap.parse_args()
    n = 2 * a.nup / (a.lx * a.ly)
    K = cb.checkerboard_hopping(a.lx, a.ly, a.t0, a.t1, -a.delta)
    q = CPMC(a.lx, a.ly, a.nup, a.ndn, U=a.U, dt=a.dt, nwalkers=a.nw, seed=a.seed, K=K, K_dn=None)
    r = q.run_bp_chi_spin(nequil=a.nequil, nblocks=a.nblocks, bp=a.bp)
    chi_q = np.array(r["chi_q"])              # (lx, ly) static chi_zz(q), per site
    fn = os.path.join(a.outdir, f"chizz_L{a.lx}_n{n:.3f}_d{a.delta}_s{a.seed}.npy")
    np.save(fn, chi_q)
    imax = np.unravel_index(chi_q.argmax(), chi_q.shape)
    print(f"L={a.lx} n={n:.3f} delta={a.delta} | chi_zz(pi,pi)={r['chi_pipi']:.4f}+/-{r['chi_pipi_err']:.4f}"
          f"  chi_zz(0,0)={r['chi_q0']:.4f}  max={chi_q.max():.4f} at grid{imax}  -> {fn}")


if __name__ == "__main__":
    main()
