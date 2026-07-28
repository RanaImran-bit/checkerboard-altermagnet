#!/usr/bin/env python3
"""ED GATE for the checkerboard DYNAMIC pairing susceptibility VERTEX (s, dx2-y2, dxy).

Adapts validate_chid.py to (1) the checkerboard hopping and (2) the dxy diagonal-bond
channel, and compares the connected-vertex chi (ED vs CP-AFQMC). THE POINT: confirm the
SIGN of the dxy connected vertex, which the equal-time vertex and the small-cluster
susceptibility cannot be trusted on without ED. If ED and CP-AFQMC AGREE on the dxy
vertex sign -> the machinery is sign-faithful (the dxy suppression, if seen, is real).
If they DISAGREE -> constrained-path bias flips it (as in the appeal's equal-time case).

C_a(tau) = <Delta_a(tau) Delta_a^dag(0)>,  Delta_a^dag = sum_m sum_delta f_a(delta) c^+_{m up} c^+_{m+delta dn}
chi_a = int_0^tmax C_a(tau) dtau  (same window as the CP-AFQMC back-prop).
Full-Fock dense ED => keep to <= 6 sites (2x3 / 3x2). Needs QuSpin: run on the
workstation/cluster (QuSpin import hangs on the laptop).

  python pyqmc/validate_checkerboard_chid.py --lx 3 --ly 2 --nup 2 --ndn 2 --U 4 --delta 0.2
"""
from __future__ import annotations
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")
import argparse, sys
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from cpqmc import CPMC
import checkerboard as cb


def ed_cb_chid(lx, ly, nup, ndn, t0, t1, t2, U, taus, Ffac):
    """Exact windowed C_a(tau) full/bubble/vertex for the checkerboard, per channel in
    Ffac (form-factor matrices in the i = x*ly + y convention of checkerboard.py)."""
    from quspin.basis import spinful_fermion_basis_general
    from quspin.operators import hamiltonian
    N = lx * ly
    basis = spinful_fermion_basis_general(N)                 # full Fock space
    nc = dict(check_pcon=False, check_symm=False, check_herm=False)
    tk = cb.checkerboard_hopping(lx, ly, t0, t1, t2)         # same K as the CP-AFQMC run
    hop = [[tk[i, j], i, j] for i in range(N) for j in range(N) if abs(tk[i, j]) > 1e-14]
    inter = [[U, i, i] for i in range(N)]
    static = [["+-|", hop], ["|+-", hop], ["n|n", inter]]
    Hf = hamiltonian(static, [], basis=basis, dtype=np.float64, **nc).toarray()
    fu = hamiltonian([["n|", [[1.0, i] for i in range(N)]]], [], basis=basis, dtype=np.float64, **nc).diagonal().real
    fd = hamiltonian([["|n", [[1.0, i] for i in range(N)]]], [], basis=basis, dtype=np.float64, **nc).diagonal().real
    m0  = np.where((np.rint(fu) == nup)     & (np.rint(fd) == ndn))[0]
    m2  = np.where((np.rint(fu) == nup + 1) & (np.rint(fd) == ndn + 1))[0]
    m1u = np.where((np.rint(fu) == nup + 1) & (np.rint(fd) == ndn))[0]
    m1d = np.where((np.rint(fu) == nup)     & (np.rint(fd) == ndn + 1))[0]
    E0s, V0s = np.linalg.eigh(Hf[np.ix_(m0, m0)])
    psi0 = np.zeros(basis.Ns); psi0[m0] = V0s[:, 0]; E0 = float(E0s[0])
    E2s, V2s = np.linalg.eigh(Hf[np.ix_(m2, m2)])            # +pair sector (Lehmann full)
    E1u, V1u = np.linalg.eigh(Hf[np.ix_(m1u, m1u)]); dE1u = E1u - E0
    E1d, V1d = np.linalg.eigh(Hf[np.ix_(m1d, m1d)]); dE1d = E1d - E0
    Mu = np.zeros((N, len(m1u))); Md = np.zeros((N, len(m1d)))
    for a in range(N):
        cu = hamiltonian([["+|", [[1.0, a]]]], [], basis=basis, dtype=np.float64, **nc).toarray() @ psi0
        cd = hamiltonian([["|+", [[1.0, a]]]], [], basis=basis, dtype=np.float64, **nc).toarray() @ psi0
        Mu[a, :] = V1u.T @ cu[m1u]; Md[a, :] = V1d.T @ cd[m1d]
    res = {}
    for tag, F in Ffac.items():
        terms = [[F[m, j], m, j] for m in range(N) for j in range(N) if abs(F[m, j]) > 1e-12]
        Op = hamiltonian([["+|+", terms]], [], basis=basis, dtype=np.float64, **nc).toarray()
        a = V2s.T @ (Op @ psi0)[m2]; dE = E2s - E0
        Cfull = np.array([float(np.sum(a ** 2 * np.exp(-d * dE))) for d in taus])
        Cbub = np.empty(len(taus))
        for k, tau in enumerate(taus):
            Gu = (Mu * np.exp(-dE1u * tau)) @ Mu.T
            Gd = (Md * np.exp(-dE1d * tau)) @ Md.T
            Cbub[k] = float(np.sum(Gu * (F @ Gd @ F.T)))
        res[tag] = (Cfull, Cbub, Cfull - Cbub)               # full, bubble, vertex
    return E0, res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lx", type=int, default=3); ap.add_argument("--ly", type=int, default=2)
    ap.add_argument("--nup", type=int, default=2); ap.add_argument("--ndn", type=int, default=2)
    ap.add_argument("--U", type=float, default=4.0)
    ap.add_argument("--t0", type=float, default=-1.0); ap.add_argument("--t1", type=float, default=0.3)
    ap.add_argument("--delta", type=float, default=0.2)
    ap.add_argument("--dt", type=float, default=0.05); ap.add_argument("--nw", type=int, default=240)
    ap.add_argument("--nequil", type=int, default=60); ap.add_argument("--nblocks", type=int, default=80)
    ap.add_argument("--bp", type=int, default=20)
    a = ap.parse_args()
    if a.lx * a.ly > 6:
        print(f"WARNING: {a.lx}x{a.ly} = {a.lx*a.ly} sites; full-Fock dense ED needs <=6 sites (32GB at 8).")
    t2 = -a.delta
    K = cb.checkerboard_hopping(a.lx, a.ly, a.t0, a.t1, t2)
    Fs, Fd = cb.nn_bond_factors(a.lx, a.ly); Fdxy = cb.diag_bond_factors(a.lx, a.ly)
    Ffac = {"s": Fs, "d": Fd, "dxy": Fdxy}

    q = CPMC(a.lx, a.ly, a.nup, a.ndn, U=a.U, dt=a.dt, nwalkers=a.nw, seed=1, K=K, K_dn=None)
    qm = cb.run_bp_chid_cb(q, Ffac, nequil=a.nequil, nblocks=a.nblocks, bp=a.bp)
    taus = np.array(qm["taus"])
    E0, ed = ed_cb_chid(a.lx, a.ly, a.nup, a.ndn, a.t0, a.t1, t2, a.U, taus, Ffac)

    def winint(C): return a.dt * (C[1:-1].sum() + 0.5 * (C[0] + C[-1]))
    print(f"# checkerboard chi ED gate  L={a.lx}x{a.ly}  n={2*a.nup/(a.lx*a.ly):.3f}  U={a.U}  delta={a.delta}")
    print(f"# windowed chi (0..{taus[-1]:.2f})   E0(ED)={E0:.5f}")
    print(f"{'ch':5} {'chiF ED':>10} {'chiF QMC':>10} | {'chiV ED':>11} {'chiV QMC':>11} {'+/-':>8}  sign")
    allok = True
    for tag in ("s", "d", "dxy"):
        Cf, Cb, Cv = ed[tag]
        edF = winint(Cf); edV = winint(Cv)
        qF = qm[f"chi_{tag}"]; qV = qm[f"chi_{tag}_vertex"]; qVe = qm[f"chi_{tag}_vertex_err"]
        sgn = "OK" if (np.sign(edV) == np.sign(qV) or abs(edV) < 2e-3) else "FLIP"
        allok = allok and (sgn == "OK")
        print(f"{tag:5} {edF:>10.4f} {qF:>10.4f} | {edV:>+11.4f} {qV:>+11.4f} {qVe:>8.4f}  {sgn}")
    print("\nVERDICT:", "vertex signs MATCH ED (machinery sign-faithful)" if allok
          else "some vertex sign FLIPPED vs ED (constrained-path bias) -- do not trust that channel")


if __name__ == "__main__":
    main()
