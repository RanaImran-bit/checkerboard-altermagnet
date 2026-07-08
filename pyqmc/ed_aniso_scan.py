#!/usr/bin/env python3
"""EXACT (ED) scan of the two-orbital altermagnet d-wave / s-wave pairing structure
factor vs the altermagnetic anisotropy (t1 fixed, t2 varied; isotropic at t2=t1).

Pairing is biased low in constrained-path QMC (validated: ~5x for d-wave), so the
trustworthy trend comes from ED on a small cluster. We also report the nearest-
neighbour vs longest-range d-wave pair correlation to separate the PAIRING STRENGTH
(total structure factor) from the LONG-RANGE pair coherence.

    source tools/env.sh
    python pyqmc/ed_aniso_scan.py --lx 2 --ly 2 --nup 2 --ndn 2 --U 4
"""
from __future__ import annotations
import os
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")
import argparse, sys
import numpy as np
sys.path.insert(0, os.path.dirname(__file__)); sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "ed"))


def am_pair_terms(lx, ly):
    """Intra-orbital NN d-wave pair operator terms: list of (m, j, f_d) where
    Delta^dag adds a pair c^+_{m up} c^+_{j dn}; d-wave f = +1 (x), -1 (y)."""
    lxy = lx * ly
    def idx(x, y, orb): return orb * lxy + (x % lx) * ly + (y % ly)
    out = []
    for orb in (0, 1):
        for x in range(lx):
            for y in range(ly):
                m = idx(x, y, orb)
                for (dx, dy, fd) in ((1, 0, +1.0), (-1, 0, +1.0), (0, 1, -1.0), (0, -1, -1.0)):
                    out.append((m, idx(x + dx, y + dy, orb), fd, x, y, dx, dy, orb))
    return out


def scan_point(lx, ly, nup, ndn, t1, t2, t3, t4, U):
    from quspin.operators import hamiltonian
    from altermagnet_ed import _build_static
    static, basis, nsites = _build_static(lx, ly, nup, ndn, t1, t2, t3, t4, U, 0.0, 0.0)
    nc = dict(check_pcon=False, check_symm=False, check_herm=False)
    H = hamiltonian(static, [], basis=basis, dtype=np.float64, **nc)
    w, V = np.linalg.eigh(H.toarray()); psi0 = V[:, 0]; E0 = float(w[0])
    bonds = am_pair_terms(lx, ly)
    lxy = lx * ly
    def xy(site):
        loc = site % lxy; return loc // ly, loc % ly
    def dist(m, n):
        x0, y0 = xy(m); x1, y1 = xy(n)
        return min((x1 - x0) % lx, (x0 - x1) % lx) + min((y1 - y0) % ly, (y0 - y1) % ly)
    # structure factor S_a = <O_a^dag O_a> and its decomposition by the distance R
    # between the two pair CENTRES (site m of each bond): S_a = sum_R P_a(R).
    s_t, d_by_R = [], {}
    for (m, j, fd, *_ ) in bonds:
        for (nn, k, fdn, *_2) in bonds:
            s_t.append([1.0, m, nn, j, k])
            d_by_R.setdefault(dist(m, nn), []).append([fd * fdn, m, nn, j, k])
    OOs = hamiltonian([["+-|+-", s_t]], [], basis=basis, dtype=np.float64, **nc)
    Ss = float(psi0 @ OOs.dot(psi0))
    Pd_R = {}
    for R, terms in d_by_R.items():
        Op = hamiltonian([["+-|+-", terms]], [], basis=basis, dtype=np.float64, **nc)
        Pd_R[R] = float(psi0 @ Op.dot(psi0))
    Sd = float(sum(Pd_R.values()))
    return E0, Sd, Ss, Pd_R


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lx", type=int, default=2); ap.add_argument("--ly", type=int, default=2)
    ap.add_argument("--nup", type=int, default=2); ap.add_argument("--ndn", type=int, default=2)
    ap.add_argument("--t1", type=float, default=-1.0); ap.add_argument("--t3", type=float, default=0.0)
    ap.add_argument("--t4", type=float, default=0.0); ap.add_argument("--U", type=float, default=4.0)
    ap.add_argument("--t2list", type=float, nargs="+",
                    default=[-1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2])
    a = ap.parse_args()
    print(f"# altermagnet {a.lx}x{a.ly} nup={a.nup} ndn={a.ndn} U={a.U} t1={a.t1} (t3={a.t3} t4={a.t4})")
    print(f"# anisotropy alpha = 1 - t2/t1  (alpha=0 isotropic)")
    Rs = None
    rows = []
    for t2 in a.t2list:
        E0, Sd, Ss, PdR = scan_point(a.lx, a.ly, a.nup, a.ndn, a.t1, t2, a.t3, a.t4, a.U)
        if Rs is None:
            Rs = sorted(PdR)
        rows.append((t2, E0, Sd, Ss, PdR))
    hdr = f"{'t2':>6} {'alpha':>6} {'E0':>9} {'S_d':>9} {'S_s':>9} " + " ".join(f"Pd(R={R})".rjust(11) for R in Rs)
    print(hdr)
    for (t2, E0, Sd, Ss, PdR) in rows:
        alpha = 1 - t2 / a.t1
        line = f"{t2:>6.2f} {alpha:>6.2f} {E0:>9.4f} {Sd:>9.4f} {Ss:>9.4f} " + \
               " ".join(f"{PdR.get(R, float('nan')):>11.5f}" for R in Rs)
        print(line)


if __name__ == "__main__":
    main()
