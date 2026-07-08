#!/usr/bin/env python3
"""Validate the SINGLE-BAND spin-dependent anisotropic-hopping altermagnet (the
manuscript's primary model) in pyqmc against exact ED: FULL and CONNECTED/VERTEX
d-wave & s-wave pairing structure factors and the AFM spin structure factor
S(pi,pi). The anisotropy is tam (mc2duph.f90): up-spin hops t0-tam on x and t0+tam
on y; down-spin rotated 90 deg.  tam=0 -> ordinary single-band Hubbard.

  full P_a   = <O_a^dag O_a>,  O_a^dag = sum_m sum_delta f_a(delta) c^+_{m up} c^+_{m+d dn}
  disc P_a   = sum_{m,n} G_up[m,n] (F_a G_dn F_a^T)[m,n]   (from the 1-particle G)
  vertex P_a = full - disc                                  (the manuscript quantity)

    source tools/env.sh
    python pyqmc/validate_am_single.py --lx 2 --ly 2 --nup 2 --ndn 2 --U 4 --tam 0.3
"""
from __future__ import annotations
import os
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")
import argparse, sys
import numpy as np
sys.path.insert(0, os.path.dirname(__file__)); sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "ed"))


def bonds_single(lx, ly):
    """(m, j, f_s, f_d) NN bonds, single-band index m = x*ly + y (pyqmc convention),
    f_d = +1 on x bonds, -1 on y bonds."""
    def idx(x, y): return (x % lx) * ly + (y % ly)
    out = []
    for x in range(lx):
        for y in range(ly):
            m = idx(x, y)
            for (dx, dy, fd) in ((1, 0, +1.0), (-1, 0, +1.0), (0, 1, -1.0), (0, -1, -1.0)):
                out.append((m, idx(x + dx, y + dy), 1.0, fd))
    return out


def ed_am_single(lx, ly, nup, ndn, t0, tam, U, t1=0.0):
    from quspin.basis import spinful_fermion_basis_general
    from quspin.operators import hamiltonian
    sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
    from cpqmc import am_hopping
    n = lx * ly
    Ku, Kd = am_hopping(lx, ly, t0, tam, t1)
    basis = spinful_fermion_basis_general(n, Nf=(nup, ndn))
    nc = dict(check_pcon=False, check_symm=False, check_herm=False)
    up_hop = [[Ku[i, j], i, j] for i in range(n) for j in range(n) if abs(Ku[i, j]) > 1e-14]
    dn_hop = [[Kd[i, j], i, j] for i in range(n) for j in range(n) if abs(Kd[i, j]) > 1e-14]
    inter = [[U, i, i] for i in range(n)]
    static = [["+-|", up_hop], ["|+-", dn_hop], ["n|n", inter]]
    H = hamiltonian(static, [], basis=basis, dtype=np.float64, **nc)
    w, V = np.linalg.eigh(H.toarray()); psi0 = V[:, 0]; E0 = float(w[0])

    bonds = bonds_single(lx, ly)
    Fs = np.zeros((n, n)); Fd = np.zeros((n, n))
    for (m, j, fs, fd) in bonds:
        Fs[m, j] += fs; Fd[m, j] += fd
    # FULL pairing via "+-|+-"
    s_t, d_t = [], []
    for (m, j, fsm, fdm) in bonds:
        for (nn, k, fsn, fdn) in bonds:
            s_t.append([fsm * fsn, m, nn, j, k]); d_t.append([fdm * fdn, m, nn, j, k])
    full_s = float(psi0 @ hamiltonian([["+-|+-", s_t]], [], basis=basis, dtype=np.float64, **nc).dot(psi0))
    full_d = float(psi0 @ hamiltonian([["+-|+-", d_t]], [], basis=basis, dtype=np.float64, **nc).dot(psi0))
    # single-particle Green's G_up[i,j]=<c^+_i c_j>, G_dn
    Gu = np.zeros((n, n)); Gd = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            Gu[i, j] = float(psi0 @ hamiltonian([["+-|", [[1.0, i, j]]]], [], basis=basis, dtype=np.float64, **nc).dot(psi0))
            Gd[i, j] = float(psi0 @ hamiltonian([["|+-", [[1.0, i, j]]]], [], basis=basis, dtype=np.float64, **nc).dot(psi0))
    disc_s = float((Gu * (Fs @ Gd @ Fs.T)).sum())
    disc_d = float((Gu * (Fd @ Gd @ Fd.T)).sum())
    # AFM S(pi,pi) = <Sstag^2>/N
    stag = [[0.5 * (-1.0) ** ((i // ly) + (i % ly)), i] for i in range(n)]
    Sst = hamiltonian([["n|", stag], ["|n", [[-c, i] for (c, i) in stag]]], [], basis=basis,
                      dtype=np.float64, **nc)
    Sp = Sst.dot(psi0); Sq = float(Sp @ Sp) / n
    return dict(E0=E0, full_s=full_s, full_d=full_d, vert_s=full_s - disc_s,
                vert_d=full_d - disc_d, Sq=Sq)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lx", type=int, default=2); ap.add_argument("--ly", type=int, default=2)
    ap.add_argument("--nup", type=int, default=2); ap.add_argument("--ndn", type=int, default=2)
    ap.add_argument("--t0", type=float, default=1.0); ap.add_argument("--tam", type=float, default=0.0)
    ap.add_argument("--t1", type=float, default=0.0)
    ap.add_argument("--U", type=float, default=4.0); ap.add_argument("--dt", type=float, default=0.01)
    ap.add_argument("--nw", type=int, default=400); ap.add_argument("--nequil", type=int, default=200)
    ap.add_argument("--nblocks", type=int, default=50); ap.add_argument("--bp", type=int, default=14)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()

    from cpqmc import CPMC, am_hopping
    Ku, Kd = am_hopping(a.lx, a.ly, a.t0, a.tam, a.t1)
    q = CPMC(a.lx, a.ly, a.nup, a.ndn, U=a.U, dt=a.dt, nwalkers=a.nw, seed=a.seed, K=Ku, K_dn=Kd)
    qm = q.run_bp_pairmag(nequil=a.nequil, nblocks=a.nblocks, bp=a.bp)
    ed = ed_am_single(a.lx, a.ly, a.nup, a.ndn, a.t0, a.tam, a.U, a.t1)

    print(f"# single-band altermagnet {a.lx}x{a.ly} nup={a.nup} ndn={a.ndn} U={a.U} tam={a.tam}; ED E0={ed['E0']:.5f}")
    print(f"{'observable':22s} {'ED':>11s} {'pyqmc':>11s} {'+/-err':>9s} {'z':>6s}")
    rows = [("d-wave FULL", "pair_dwave", ed["full_d"]),
            ("d-wave VERTEX", "pair_dwave_vertex", ed["vert_d"]),
            ("s-wave FULL", "pair_swave", ed["full_s"]),
            ("s-wave VERTEX", "pair_swave_vertex", ed["vert_s"]),
            ("S(pi,pi)", "Sq_pipi", ed["Sq"])]
    for label, key, edv in rows:
        p = qm[key]["value"]; err = qm[key]["error"]
        z = abs(p - edv) / err if err > 0 else float("nan")
        print(f"{label:22s} {edv:>11.5f} {p:>11.5f} {err:>9.5f} {z:>6.2f}")


if __name__ == "__main__":
    main()
