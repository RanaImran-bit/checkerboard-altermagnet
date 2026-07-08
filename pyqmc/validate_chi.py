#!/usr/bin/env python3
"""Phase 6 WS2 validation: imaginary-time-displaced STAGGERED spin correlation
C(tau) = <O(tau) O(0)>, O = sum_i (-1)^{x+y} S^z_i, and the static staggered
susceptibility chi_s = integral_0^inf C(tau) dtau, pyqmc CPMC vs exact ED
(Lehmann representation), single-band Hubbard.

    source tools/env.sh
    python pyqmc/validate_chi.py --lx 4 --ly 4 --nup 1 --ndn 1 --U 0   # exact gate
    python pyqmc/validate_chi.py --lx 4 --ly 4 --nup 1 --ndn 1 --U 3

ED (Lehmann): C(tau) = sum_n |<n|O|0>|^2 exp(-(E_n-E_0) tau),
              chi_s   = sum_{n>0} |<n|O|0>|^2 / (E_n - E_0).
C(tau=0) = <O^2> = N * S(pi,pi) (equal-time staggered structure factor).
Scalars are convention-invariant, so ED vs pyqmc site-index orders do not matter.
"""
from __future__ import annotations
import argparse, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "ed"))


def ed_chi(lx, ly, nup, ndn, t, U, taus):
    from quspin.operators import hamiltonian
    from hubbard_ed import build
    basis, H, _, _ = build(lx, ly, nup, ndn, t, U)
    no_check = dict(check_pcon=False, check_symm=False, check_herm=False)
    E, V = np.linalg.eigh(H.toarray())
    psi0 = V[:, 0]; E0 = E[0]
    # staggered O = sum_i 0.5 (-1)^{x+y} (n_iu - n_id)
    up, dn = [], []
    for x in range(lx):
        for y in range(ly):
            i = (x % lx) + lx * (y % ly)
            p = 0.5 * (-1.0) ** (x + y)
            up.append([p, i]); dn.append([-p, i])
    O = hamiltonian([["n|", up], ["|n", dn]], [], basis=basis, dtype=np.float64, **no_check)
    Opsi = O.dot(psi0)
    a = V.T @ Opsi                       # a_n = <n|O|0>
    dE = E - E0
    C = np.array([float(np.sum(a ** 2 * np.exp(-d * dE))) for d in taus])
    mask = dE > 1e-9
    chi = float(np.sum(a[mask] ** 2 / dE[mask]))
    return C, chi


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lx", type=int, default=4); ap.add_argument("--ly", type=int, default=4)
    ap.add_argument("--nup", type=int, default=1); ap.add_argument("--ndn", type=int, default=1)
    ap.add_argument("--U", type=float, default=0.0); ap.add_argument("--t", type=float, default=1.0)
    ap.add_argument("--dt", type=float, default=0.01); ap.add_argument("--nw", type=int, default=300)
    ap.add_argument("--nequil", type=int, default=200); ap.add_argument("--nblocks", type=int, default=40)
    ap.add_argument("--bp", type=int, default=20); ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()

    from cpqmc import CPMC
    q = CPMC(a.lx, a.ly, a.nup, a.ndn, t=a.t, U=a.U, dt=a.dt, nwalkers=a.nw, seed=a.seed)
    qm = q.run_bp_chi(nequil=a.nequil, nblocks=a.nblocks, bp=a.bp)
    taus = np.array(qm["taus"])
    Cq = np.array(qm["Ctau"]); Cqe = np.array(qm["Cerr"])
    Ce, chi_ed = ed_chi(a.lx, a.ly, a.nup, a.ndn, a.t, a.U, taus)

    print(f"# single-band {a.lx}x{a.ly}, nup={a.nup} ndn={a.ndn}, U={a.U}, t={a.t}")
    print(f"# staggered spin C(tau)=<O(tau)O(0)>, tau=l*dt (dt={a.dt}), bp={a.bp}\n")
    print(f"{'tau':>6} {'ED C':>10} {'pyqmc C':>10} {'+/-':>8} {'z':>6}")
    for l in range(0, len(taus), max(1, len(taus) // 10)):
        z = abs(Cq[l] - Ce[l]) / Cqe[l] if Cqe[l] > 0 else float("nan")
        print(f"{taus[l]:>6.3f} {Ce[l]:>10.4f} {Cq[l]:>10.4f} {Cqe[l]:>8.4f} {z:>6.2f}")
    # chi from the same trapezoidal rule on the ED curve, over the measured window
    chi_ed_win = a.dt * (Ce[1:-1].sum() + 0.5 * (Ce[0] + Ce[-1]))
    print(f"\nchi_stag (window 0..{taus[-1]:.2f}):  ED {chi_ed_win:.4f}   "
          f"pyqmc {qm['chi_stag']:.4f} +/- {qm['chi_stag_err']:.4f}")
    print(f"chi_stag (ED full integral, all tau): {chi_ed:.4f}")


if __name__ == "__main__":
    main()
