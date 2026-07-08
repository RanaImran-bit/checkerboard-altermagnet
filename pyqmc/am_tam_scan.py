#!/usr/bin/env python3
"""EXACT (ED) test of the manuscript's central claim on the SINGLE-BAND spin-
dependent altermagnet: does increasing the anisotropy tam ENHANCE the equal-time
CONNECTED (vertex) d-wave pairing correlation while SUPPRESSING the AFM structure
factor S(pi,pi)?  Scans tam at (near) half filling on small clusters. d-wave is
constrained-path-biased in QMC, so the trustworthy trend is ED.

    source tools/env.sh
    python pyqmc/am_tam_scan.py --lx 4 --ly 2 --nup 4 --ndn 4 --U 6
"""
from __future__ import annotations
import os
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")
import argparse, sys
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from validate_am_single import ed_am_single


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lx", type=int, default=4); ap.add_argument("--ly", type=int, default=2)
    ap.add_argument("--nup", type=int, default=4); ap.add_argument("--ndn", type=int, default=4)
    ap.add_argument("--t0", type=float, default=1.0); ap.add_argument("--U", type=float, default=6.0)
    ap.add_argument("--tams", type=float, nargs="+",
                    default=[0.0, 0.1, 0.2, 0.3, 0.4, 0.5])
    a = ap.parse_args()
    n = a.lx * a.ly
    print(f"# single-band altermagnet {a.lx}x{a.ly} ({n} sites) nup={a.nup} ndn={a.ndn} "
          f"U={a.U} t0={a.t0}  filling={ (a.nup+a.ndn)/n :.3f}")
    print(f"# EXACT ED. claim: tam^ enhances d-wave VERTEX, suppresses S(pi,pi).\n")
    print(f"{'tam':>5} {'E0':>10} {'d_full':>9} {'d_VERTEX':>9} {'s_full':>9} "
          f"{'s_VERTEX':>9} {'S(pi,pi)':>9}")
    for tam in a.tams:
        r = ed_am_single(a.lx, a.ly, a.nup, a.ndn, a.t0, tam, a.U)
        print(f"{tam:>5.2f} {r['E0']:>10.4f} {r['full_d']:>9.4f} {r['vert_d']:>9.4f} "
              f"{r['full_s']:>9.4f} {r['vert_s']:>9.4f} {r['Sq']:>9.5f}")


if __name__ == "__main__":
    main()
