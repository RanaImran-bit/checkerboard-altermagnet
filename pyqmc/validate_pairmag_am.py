#!/usr/bin/env python3
"""Validate the two-orbital ALTERMAGNET intra-orbital d-wave / s-wave pairing
structure factor (pyqmc run_bp_pairmag) against exact ED.

Pair operator (intra-orbital, nearest-neighbour):
  O_a = sum_orb sum_{(m,j) NN bond} f_a(bond) c_{m,up} c_{j,dn},
  f_d = +1 on x bonds, -1 on y bonds (d_{x^2-y^2});  f_s = +1.
Structure factor S_a = <O_a^dag O_a>, built as the number-conserving product
"+-|+-" on the altermagnet ED basis. The d-wave form factor is geometric, so the
anisotropy (t1 != t2) enters only through the ground state, not the operator.

    source tools/env.sh
    python pyqmc/validate_pairmag_am.py --lx 2 --ly 2 --nup 4 --ndn 4 --U 0
"""
from __future__ import annotations
import os
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")
import argparse, sys
import numpy as np
sys.path.insert(0, os.path.dirname(__file__)); sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "ed"))


def am_bonds(lx, ly):
    """(m, j, f_s, f_d) intra-orbital NN bonds; index orb-blocked orb in {0,1}:
    idx = orb*lxy + x*ly + y (matches build_hopping / pyqmc _bond_factors)."""
    lxy = lx * ly
    def idx(x, y, orb): return orb * lxy + (x % lx) * ly + (y % ly)
    out = []
    for orb in (0, 1):
        for x in range(lx):
            for y in range(ly):
                m = idx(x, y, orb)
                for (dx, dy, fd) in ((1, 0, +1.0), (-1, 0, +1.0), (0, 1, -1.0), (0, -1, -1.0)):
                    out.append((m, idx(x + dx, y + dy, orb), 1.0, fd))
    return out


def ed_pairmag_am(lx, ly, nup, ndn, t1, t2, t3, t4, U):
    from quspin.operators import hamiltonian
    from altermagnet_ed import _build_static
    static, basis, nsites = _build_static(lx, ly, nup, ndn, t1, t2, t3, t4, U, 0.0, 0.0)
    nc = dict(check_pcon=False, check_symm=False, check_herm=False)
    H = hamiltonian(static, [], basis=basis, dtype=np.float64, **nc)
    w, V = np.linalg.eigh(H.toarray()); psi0 = V[:, 0]
    mb = am_bonds(lx, ly)
    s_terms, d_terms = [], []
    for (m, j, fsm, fdm) in mb:
        for (nn, k, fsn, fdn) in mb:
            s_terms.append([fsm * fsn, m, nn, j, k])
            d_terms.append([fdm * fdn, m, nn, j, k])
    OOs = hamiltonian([["+-|+-", s_terms]], [], basis=basis, dtype=np.float64, **nc)
    OOd = hamiltonian([["+-|+-", d_terms]], [], basis=basis, dtype=np.float64, **nc)
    return float(w[0]), float(psi0 @ OOs.dot(psi0)), float(psi0 @ OOd.dot(psi0))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lx", type=int, default=2); ap.add_argument("--ly", type=int, default=2)
    ap.add_argument("--nup", type=int, default=4); ap.add_argument("--ndn", type=int, default=4)
    ap.add_argument("--t1", type=float, default=-1.0); ap.add_argument("--t2", type=float, default=-1.0)
    ap.add_argument("--t3", type=float, default=0.0); ap.add_argument("--t4", type=float, default=0.0)
    ap.add_argument("--U", type=float, default=0.0); ap.add_argument("--dt", type=float, default=0.01)
    ap.add_argument("--nw", type=int, default=300); ap.add_argument("--nequil", type=int, default=150)
    ap.add_argument("--nblocks", type=int, default=40); ap.add_argument("--bp", type=int, default=12)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()

    from cpqmc import CPMC
    from altermagnet_ed import build_hopping
    K = build_hopping(a.lx, a.ly, a.t1, a.t2, a.t3, a.t4)
    q = CPMC(a.lx, a.ly, a.nup, a.ndn, U=a.U, dt=a.dt, nwalkers=a.nw, seed=a.seed, K=K)
    qm = q.run_bp_pairmag(nequil=a.nequil, nblocks=a.nblocks, bp=a.bp)
    E0, Ss, Sd = ed_pairmag_am(a.lx, a.ly, a.nup, a.ndn, a.t1, a.t2, a.t3, a.t4, a.U)

    print(f"# altermagnet {a.lx}x{a.ly}, nup={a.nup} ndn={a.ndn}, t1={a.t1} t2={a.t2} U={a.U}; ED E0={E0:.5f}")
    print(f"{'observable':18s} {'ED':>11s} {'pyqmc':>11s} {'+/-err':>9s} {'z':>6s}")
    for key, label, ed in (("pair_dwave", "d-wave pairing SF", Sd),
                           ("pair_swave", "s-wave pairing SF", Ss)):
        p = qm[key]["value"]; err = qm[key]["error"]
        z = abs(p - ed) / err if err > 0 else float("nan")
        print(f"{label:18s} {ed:>11.5f} {p:>11.5f} {err:>9.5f} {z:>6.2f}")


if __name__ == "__main__":
    main()
