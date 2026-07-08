#!/usr/bin/env python3
"""Phase 6 WS1b/1c validation: singlet pairing (s + d_{x^2-y^2}) structure factors
and the antiferromagnetic spin structure factor S(pi,pi), pyqmc CPMC vs exact ED,
on the single-band square-lattice Hubbard model.

Convention-invariant scalars are used (total pairing structure factor summed over
all site pairs; S(pi,pi)), so the differing ED vs pyqmc site-index conventions do
not matter -- each code builds the operators with its own (x,y) map.

    source tools/env.sh
    python pyqmc/validate_pairmag.py --lx 4 --ly 4 --nup 2 --ndn 2 --U 4

ED operators:
  O_alpha = sum_m sum_{delta} f_alpha(delta) c_{m,up} c_{m+delta,dn}   (singlet pair)
  pairing structure factor  S_alpha = <O_alpha^dagger O_alpha> = || O_alpha|psi> ||^2
  staggered spin  Sstag = sum_i (-1)^{x+y} S^z_i ,  S(pi,pi) = <Sstag^2>/N
"""
from __future__ import annotations
import argparse, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "ed"))


def ed_pairmag(lx, ly, nup, ndn, t, U):
    from quspin.operators import hamiltonian
    from hubbard_ed import build
    basis, H, K, V = build(lx, ly, nup, ndn, t, U)
    no_check = dict(check_pcon=False, check_symm=False, check_herm=False)
    w, Vv = np.linalg.eigh(H.toarray())
    psi = Vv[:, 0]
    N = lx * ly

    def site(x, y): return (x % lx) + lx * (y % ly)   # hubbard_ed convention
    # Singlet pair operator O_a = sum_m sum_delta f_a(delta) c_{m,up} c_{m+delta,dn}.
    # O_a changes particle number, so in a fixed-Nf basis we build the number-
    # conserving product O_a^dagger O_a directly:
    #   O_a^dagger O_a = sum f_a(m,j) f_a(n,k) [c^+_{m up} c_{n up}][c^+_{j dn} c_{k dn}]
    # which is the QuSpin spinful string "+-|+-" with indices (m, n, j, k).
    bonds = [(+1, 0, +1.0), (-1, 0, +1.0), (0, +1, -1.0), (0, -1, -1.0)]   # (dx,dy,fd)
    # (site m, neighbour j, form factors fs=1, fd) for every bond
    mb = []   # list of (m, j, fs, fd)
    for x in range(lx):
        for y in range(ly):
            m = site(x, y)
            for (dx, dy, fd) in bonds:
                mb.append((m, site(x + dx, y + dy), 1.0, fd))
    s_terms, d_terms = [], []
    for (m, j, fsm, fdm) in mb:
        for (nn, k, fsn, fdn) in mb:
            s_terms.append([fsm * fsn, m, nn, j, k])
            d_terms.append([fdm * fdn, m, nn, j, k])
    OOs = hamiltonian([["+-|+-", s_terms]], [], basis=basis, dtype=np.float64, **no_check)
    OOd = hamiltonian([["+-|+-", d_terms]], [], basis=basis, dtype=np.float64, **no_check)
    pair_s = float(psi @ OOs.dot(psi)); pair_d = float(psi @ OOd.dot(psi))

    # staggered S^z = 0.5 sum_i (-1)^{x+y} (n_iu - n_id)
    stag = []
    for x in range(lx):
        for y in range(ly):
            stag.append([0.5 * (-1.0) ** (x + y), site(x, y)])
    stag_dn = [[-c, i] for (c, i) in stag]
    Sst = hamiltonian([["n|", stag], ["|n", stag_dn]], [], basis=basis,
                      dtype=np.float64, **no_check)
    Sp = Sst.dot(psi)
    Sq = float(Sp @ Sp) / N

    # local moment vector for moment^2
    du = [np.real(hamiltonian([["n|", [[1.0, i]]]], [], basis=basis,
                              dtype=np.float64, **no_check).diagonal()) for i in range(N)]
    dd = [np.real(hamiltonian([["|n", [[1.0, i]]]], [], basis=basis,
                              dtype=np.float64, **no_check).diagonal()) for i in range(N)]
    p2 = psi ** 2
    mom = np.array([float((p2 * (du[i] - dd[i])).sum()) for i in range(N)])
    return {"e0": float(w[0]), "pair_swave": pair_s, "pair_dwave": pair_d,
            "Sq_pipi": Sq, "moment2": float(np.sum(mom ** 2))}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lx", type=int, default=4); ap.add_argument("--ly", type=int, default=4)
    ap.add_argument("--nup", type=int, default=2); ap.add_argument("--ndn", type=int, default=2)
    ap.add_argument("--U", type=float, default=4.0); ap.add_argument("--t", type=float, default=1.0)
    ap.add_argument("--nw", type=int, default=300); ap.add_argument("--dt", type=float, default=0.01)
    ap.add_argument("--nequil", type=int, default=200); ap.add_argument("--nblocks", type=int, default=40)
    ap.add_argument("--bp", type=int, default=12); ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()

    print(f"# single-band {a.lx}x{a.ly}, nup={a.nup} ndn={a.ndn}, U={a.U}, t={a.t}")
    ed = ed_pairmag(a.lx, a.ly, a.nup, a.ndn, a.t, a.U)
    print(f"ED ground state E0 = {ed['e0']:.6f}")

    from cpqmc import CPMC
    q = CPMC(a.lx, a.ly, a.nup, a.ndn, t=a.t, U=a.U, dt=a.dt, nwalkers=a.nw, seed=a.seed)
    qm = q.run_bp_pairmag(nequil=a.nequil, nblocks=a.nblocks, bp=a.bp)

    keys = [("pair_swave", "s-wave pairing SF"), ("pair_dwave", "d-wave pairing SF"),
            ("Sq_pipi", "spin S(pi,pi)"), ("moment2", "sum_i <m_i^z>^2")]
    print(f"\n{'observable':22s} {'ED':>11s} {'pyqmc':>11s} {'+/-err':>9s} "
          f"{'z':>6s}  result")
    for k, label in keys:
        e = ed[k]; p = qm[k]["value"]; err = qm[k]["error"]
        z = abs(p - e) / err if err > 0 else float("nan")
        tag = "PASS" if (err > 0 and z < 2.5) else ("~0" if abs(e) < 1e-9 else "check")
        print(f"{label:22s} {e:>11.5f} {p:>11.5f} {err:>9.5f} {z:>6.2f}  {tag}")


if __name__ == "__main__":
    main()
