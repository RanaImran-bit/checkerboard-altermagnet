#!/usr/bin/env python3
"""Phase 6: validate the unequal-time singlet PAIRING susceptibility (s-wave and
d_{x^2-y^2}) of pyqmc CPMC against exact ED (Lehmann), single-band Hubbard.

C_a(tau) = <Delta_a(tau) Delta_a^dag(0)>,
  Delta_a^dag = sum_m sum_delta f_a(delta) c^+_{m up} c^+_{m+delta dn},
chi_a = integral_0^inf C_a(tau) dtau.  The pair operator changes particle number,
so ED uses the FULL Fock space (small clusters only) and the (N+2)-particle
spectrum:  C_a(tau) = sum_n |<n|Delta_a^dag|0>|^2 exp(-(E_n-E0) tau).
C_a(tau=0) = <Delta_a Delta_a^dag> is the equal-time pair structure factor.

    source tools/env.sh
    python pyqmc/validate_chid.py --lx 2 --ly 2 --nup 2 --ndn 2 --U 0   # exact gate
    python pyqmc/validate_chid.py --lx 2 --ly 2 --nup 2 --ndn 2 --U 4
"""
from __future__ import annotations
import os
# must precede numpy/quspin import: conda openblas (libgomp) + quspin (llvm libomp)
# clash; force single-thread and allow the duplicate OpenMP runtime to load.
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("NUMBA_NUM_THREADS", "1")
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")
import argparse, sys
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "ed"))


def bonds_xy(lx, ly):
    """(m, j, f_s, f_d) for the 4 NN bonds of each site; site index i = x + lx*y;
    d-wave f_d = +1 on x bonds, -1 on y bonds; s-wave f_s = +1."""
    def site(x, y): return (x % lx) + lx * (y % ly)
    out = []
    for x in range(lx):
        for y in range(ly):
            m = site(x, y)
            for (dx, dy, fd) in ((1, 0, +1.0), (-1, 0, +1.0), (0, 1, -1.0), (0, -1, -1.0)):
                out.append((m, site(x + dx, y + dy), 1.0, fd))
    return out


def ed_chid(lx, ly, nup, ndn, t, U, taus):
    from quspin.basis import spinful_fermion_basis_general
    from quspin.operators import hamiltonian
    N = lx * ly
    basis = spinful_fermion_basis_general(N)          # full Fock space
    nc = dict(check_pcon=False, check_symm=False, check_herm=False)
    # Hubbard H on the full Fock space
    hop = []
    def site(x, y): return (x % lx) + lx * (y % ly)
    for x in range(lx):
        for y in range(ly):
            i = site(x, y)
            for (jx, jy) in (((x + 1), y), (x, (y + 1))):
                if (jx % lx, jy % ly) == (x, y):
                    continue
                hop.append((i, site(jx, jy)))
    pm = [[-t, i, j] for (i, j) in hop]; mp = [[+t, i, j] for (i, j) in hop]
    inter = [[U, i, i] for i in range(N)]
    static = [["+-|", pm], ["-+|", mp], ["|+-", pm], ["|-+", mp], ["n|n", inter]]
    Hf = hamiltonian(static, [], basis=basis, dtype=np.float64, **nc).toarray()
    # per-basis-state filling (diagonal number operators) -> exact integer sectors
    fu = hamiltonian([["n|", [[1.0, i] for i in range(N)]]], [], basis=basis,
                     dtype=np.float64, **nc).diagonal().real
    fd = hamiltonian([["|n", [[1.0, i] for i in range(N)]]], [], basis=basis,
                     dtype=np.float64, **nc).diagonal().real
    m0 = np.where((np.rint(fu) == nup) & (np.rint(fd) == ndn))[0]      # (nup,ndn)
    m2 = np.where((np.rint(fu) == nup + 1) & (np.rint(fd) == ndn + 1))[0]  # +pair
    m1u = np.where((np.rint(fu) == nup + 1) & (np.rint(fd) == ndn))[0]      # +1 up
    m1d = np.where((np.rint(fu) == nup) & (np.rint(fd) == ndn + 1))[0]      # +1 dn
    # ground state in the (nup,ndn) sector
    E0s, V0s = np.linalg.eigh(Hf[np.ix_(m0, m0)])
    psi0 = np.zeros(basis.Ns); psi0[m0] = V0s[:, 0]; E0 = float(E0s[0])
    # (N+2)-particle spectrum for the full-susceptibility Lehmann sum
    E2s, V2s = np.linalg.eigh(Hf[np.ix_(m2, m2)])
    # (N+-1) spectra for the single-particle time-displaced GFs -> the BUBBLE
    E1u, V1u = np.linalg.eigh(Hf[np.ix_(m1u, m1u)]); dE1u = E1u - E0
    E1d, V1d = np.linalg.eigh(Hf[np.ix_(m1d, m1d)]); dE1d = E1d - E0
    # M_u[a,n] = <n^{N+1up}| c^+_{a up} |0>;  G_up(tau)_ab = sum_n M[a,n]M[b,n]e^{-dE_n tau}
    Mu = np.zeros((N, len(m1u))); Md = np.zeros((N, len(m1d)))
    for a in range(N):
        cu = hamiltonian([["+|", [[1.0, a]]]], [], basis=basis, dtype=np.float64, **nc).toarray() @ psi0
        cd = hamiltonian([["|+", [[1.0, a]]]], [], basis=basis, dtype=np.float64, **nc).toarray() @ psi0
        Mu[a, :] = V1u.T @ cu[m1u]; Md[a, :] = V1d.T @ cd[m1d]
    # bond form-factor matrices (same bonds_xy convention as the full operator)
    Fmat = {"s": np.zeros((N, N)), "d": np.zeros((N, N))}
    for (m, j, fs, fd) in bonds_xy(lx, ly):
        Fmat["s"][m, j] += fs; Fmat["d"][m, j] += fd
    res = {}
    for tag, fkey in (("s", 2), ("d", 3)):
        terms = [[b[fkey], b[0], b[1]] for b in bonds_xy(lx, ly)]
        Op = hamiltonian([["+|+", terms]], [], basis=basis, dtype=np.float64, **nc).toarray()
        a = V2s.T @ (Op @ psi0)[m2]; dE = E2s - E0
        Cfull = np.array([float(np.sum(a ** 2 * np.exp(-d * dE))) for d in taus])
        F = Fmat[tag]
        Cbub = np.empty(len(taus))
        for k, tau in enumerate(taus):
            Gu = (Mu * np.exp(-dE1u * tau)) @ Mu.T   # G_up(tau) = <c_a(t) c^+_b>
            Gd = (Md * np.exp(-dE1d * tau)) @ Md.T
            Cbub[k] = float(np.sum(Gu * (F @ Gd @ F.T)))
        res[tag] = (Cfull, Cbub, Cfull - Cbub)       # full, bubble, vertex
    return E0, res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lx", type=int, default=2); ap.add_argument("--ly", type=int, default=2)
    ap.add_argument("--nup", type=int, default=2); ap.add_argument("--ndn", type=int, default=2)
    ap.add_argument("--U", type=float, default=4.0); ap.add_argument("--t", type=float, default=1.0)
    ap.add_argument("--dt", type=float, default=0.01); ap.add_argument("--nw", type=int, default=400)
    ap.add_argument("--nequil", type=int, default=200); ap.add_argument("--nblocks", type=int, default=40)
    ap.add_argument("--bp", type=int, default=20); ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()

    from cpqmc import CPMC
    q = CPMC(a.lx, a.ly, a.nup, a.ndn, t=a.t, U=a.U, dt=a.dt, nwalkers=a.nw, seed=a.seed)
    qm = q.run_bp_chid(nequil=a.nequil, nblocks=a.nblocks, bp=a.bp)
    taus = np.array(qm["taus"])
    E0, ed = ed_chid(a.lx, a.ly, a.nup, a.ndn, a.t, a.U, taus)

    def winint(C):
        return a.dt * (C[1:-1].sum() + 0.5 * (C[0] + C[-1]))
    print(f"# single-band {a.lx}x{a.ly}, nup={a.nup} ndn={a.ndn}, U={a.U}; ED E0={E0:.6f}, bp={a.bp}")
    print(f"# windowed chi (0..{taus[-1]:.2f}): FULL and VERTEX (=full-bubble) susceptibility")
    print(f"{'channel':10} {'chi FULL ED':>12} {'pyqmc':>10} | {'chi VERTEX ED':>13} {'pyqmc':>10} {'+/-':>8}")
    for tag, name in (("d", "d-wave"), ("s", "s-wave")):
        Cfull, Cbub, Cvtx = ed[tag]
        ed_full = winint(Cfull); ed_vtx = winint(Cvtx)
        q_full = qm[f"chi_{tag}"]; q_vtx = qm[f"chi_{tag}_vertex"]; q_vtx_e = qm[f"chi_{tag}_vertex_err"]
        print(f"{name:10} {ed_full:>12.4f} {q_full:>10.4f} | {ed_vtx:>13.4f} {q_vtx:>10.4f} {q_vtx_e:>8.4f}")
    print("\n# per-tau VERTEX C(tau) (d-wave):")
    Cvtx_ed = ed["d"][2]; Cvtx_q = np.array(qm["Cdtau_vertex"])
    print(f"{'tau':>6} {'ED vtx':>10} {'pyqmc vtx':>11}")
    for l in range(0, len(taus), max(1, len(taus) // 8)):
        print(f"{taus[l]:>6.3f} {Cvtx_ed[l]:>10.4f} {Cvtx_q[l]:>11.4f}")


if __name__ == "__main__":
    main()
