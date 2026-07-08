#!/usr/bin/env python3
"""Momentum-resolved MAGNETIC (spin) susceptibility for the combined altermagnet model
(spin-dependent NN anisotropy tam=tA + spin-dependent NNN t1=t'), finite-T DQMC.

    chi_zz(q) = (1/N) sum_{ij} e^{iq(ri-rj)} int_0^beta dtau <S^z_i(tau) S^z_j(0)>,
    S^z_i = (n_iu - n_id)/2.

This is the tau-integrated (thermodynamic/Kubo) extension of the equal-time spin
structure factor S^z(k) used in the two altermagnet papers (PRB 113,134443 (2026) and
the tA+t' manuscript). Per HS configuration Wick gives (same-spin exchange only;
cross-spin terms factorize):

    <Sz_i(tau)Sz_j(0)>_cfg = 1/4 [ m_i(tau) m_j(0)
                                   - Gu(0,tau)_{ji} Gu(tau,0)_{ij}
                                   - Gd(0,tau)_{ji} Gd(tau,0)_{ij} ],
    m_i(tau) = n_iu(tau) - n_id(tau)  (nonzero pattern per config; the DISCONNECTED
    density part survives config-averaging and carries q-structure -- do NOT drop it).

Time-displaced Green's functions (BSS, stable big/small UDV split, docs/STABILIZATION.md):
    G(tau_l,0) =  (B(l,0)^{-1} + B(NT,l))^{-1}          [= dqmc._green_tau]
    G(0,tau_l) = -(B(NT,l)^{-1} + B(l,0))^{-1}          [same routine, chains swapped]
    G(l,l)     =  equal-time at slice l                  [= dqmc.green(s,l)]
between stabilizations propagated per slice:
    G(l+1,0) = B_l G(l,0);  G(0,l+1) = G(0,l) B_l^{-1};  G(l+1,l+1) = B_l G(l,l) B_l^{-1}.

Validation: self-contained numpy ED (Jordan-Wigner, full 4^n Fock space) computing the
same chi_zz(q) via the Lehmann/Kubo sum -- run with --ed on 2x2 (dim 256).

    python code/dqmc_py/spin_susc.py --lx 2 --ly 2 --U 4 --mu 2 --beta 2 \
        --tam 0.3 --t1 0.2 --nmeas 2000 --ed
"""
from __future__ import annotations
import os, sys, argparse
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from dqmc import DQMC, _shift_index, reduce_mat


# ----------------------------------------------------------------------------- DQMC --
class SpinDQMC(DQMC):
    """DQMC + tau-integrated spin-spin (S^z S^z) susceptibility matrix."""

    def _expKinv(self, s):
        if not hasattr(self, "_eKi"):
            self._eKi = [np.linalg.inv(self.expK[0]), np.linalg.inv(self.expK[1])]
        return self._eKi[s]

    def _green_0tau(self, s, l):
        """Stable G(0,tau_l) = -(B(NT,l)^{-1} + B(l,0))^{-1}: the _green_tau big/small
        formula with the two UDV chains swapped."""
        if l == 0:
            return self.green(s, 0) - np.eye(self.n)
        U1, D1, V1 = self._udv_chain(s, l, self.NT)    # plays the "B1" (inverted) role
        U2, D2, V2 = self._udv_chain(s, 0, l)          # plays the "B2" role
        D1b = np.where(np.abs(D1) > 1.0, D1, 1.0); D1s = np.where(np.abs(D1) > 1.0, 1.0, D1)
        D2b = np.where(np.abs(D2) > 1.0, D2, 1.0); D2s = np.where(np.abs(D2) > 1.0, 1.0, D2)
        V2i = np.linalg.inv(V2)
        INNER = (U1.T @ V2i) / D1b[:, None] / D2b[None, :] + D1s[:, None] * (V1 @ U2) * D2s[None, :]
        X = np.linalg.solve(INNER, D1s[:, None] * V1)
        return -(V2i @ (X / D2b[:, None]))

    def chi_spin(self):
        """One-configuration tau-integrated spin correlation MATRIX
        Mzz[i,j] = dt sum_l <Sz_i(tau_l) Sz_j(0)>_cfg  (rectangle rule, matches chi_pair)
        plus the equal-time matrix Szz[i,j] = <Sz_i Sz_j>_cfg."""
        n = self.n
        G00 = [self.green(0, 0), self.green(1, 0)]
        Gl0 = [G00[0].copy(), G00[1].copy()]
        G0l = [G00[0] - np.eye(n), G00[1] - np.eye(n)]
        Gll = [G00[0].copy(), G00[1].copy()]
        m0 = (1.0 - np.diag(G00[0])) - (1.0 - np.diag(G00[1]))   # m_j(0)
        Mzz = np.zeros((n, n)); Szz = None
        for l in range(self.NT):
            if l > 0:
                if l % self.nstab == 0:                    # restabilize all three
                    Gl0 = [self._green_tau(0, l), self._green_tau(1, l)]
                    G0l = [self._green_0tau(0, l), self._green_0tau(1, l)]
                    Gll = [self.green(0, l), self.green(1, l)]
                else:                                      # cheap slice propagation
                    for s in (0, 1):
                        Bv = self._Bvec(s, l - 1)
                        B = self.expK[s] * Bv[None, :]
                        Bi = (1.0 / Bv)[:, None] * self._expKinv(s)
                        Gl0[s] = B @ Gl0[s]
                        G0l[s] = G0l[s] @ Bi
                        Gll[s] = B @ Gll[s] @ Bi
            ml = (1.0 - np.diag(Gll[0])) - (1.0 - np.diag(Gll[1]))
            M = 0.25 * (np.outer(ml, m0)
                        - G0l[0].T * Gl0[0]               # [i,j] = G(0,l)[j,i] G(l,0)[i,j]
                        - G0l[1].T * Gl0[1])
            Mzz += M
            if l == 0:
                Szz = M.copy()
        return self.dt * Mzz, Szz

    def run_spin(self, nwarm=200, nmeas=400, nsub=2):
        """Sign-weighted config average of the chi_zz and equal-time Szz matrices;
        returns q-grids chi_q = P(q)/N and S_q = P(q)/N via the reduce_mat FFT."""
        for _ in range(nwarm):
            self.sweep()
        n = self.n
        Mx = np.zeros((n, n)); Mc = np.zeros((n, n))
        sw = 0.0; saw = 0.0; dens = 0.0
        for _ in range(nmeas):
            for _ in range(nsub):
                s = self.sweep()
            xm, cm = self.chi_spin()
            Mx += s * xm; Mc += s * cm
            g = self.green(0, 0); gd = self.green(1, 0)
            dens += s * float((n - np.trace(g) + n - np.trace(gd)) / n)
            sw += s; saw += abs(s)
        Mx /= sw; Mc /= sw
        shift = _shift_index(self.lx, self.ly)
        rx = reduce_mat(Mx, shift); rc = reduce_mat(Mc, shift)
        return dict(chi_q=rx["Pq"] / n, S_q=rc["Pq"] / n,      # per-site chi(q), S(q)
                    chi_max=float(rx["Pq"].max() / n), S_max=float(rc["Pq"].max() / n),
                    chi_q0=float(rx["Pq"][0, 0] / n),          # uniform susceptibility
                    Mzz=Mx, Szz=Mc, dens=dens / sw, sign=sw / saw)


# ------------------------------------------------------------------------------- ED --
def _fock_ops(nso):
    """Jordan-Wigner fermion annihilation matrices on the 2^nso Fock space (numpy only)."""
    dim = 2 ** nso
    ops = []
    occ = ((np.arange(dim)[:, None] >> np.arange(nso)[None, :]) & 1).astype(np.int8)
    for a in range(nso):
        rows = np.where(occ[:, a] == 1)[0]
        cols = rows ^ (1 << a)                      # state after removing particle a
        parity = (-1.0) ** occ[rows, :a].sum(axis=1)
        c = np.zeros((dim, dim))
        c[cols, rows] = parity                      # c|rows> = parity |cols>
        ops.append(c)
    return ops


def ed_chi_spin(lx, ly, U, mu, beta, tam=0.0, t1=0.0, tp=0.0):
    """Exact chi_zz(q) (Kubo/Lehmann) + equal-time S^z(q) on the full Fock space.
    Returns (chi_q, S_q, dens): (lx,ly) q-grids matching reduce_mat's FFT convention."""
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "pyqmc"))
    from cpqmc import am_hopping
    n = lx * ly; nso = 2 * n
    Ku, Kd = am_hopping(lx, ly, 1.0, tam, t1, tp)
    c = _fock_ops(nso)                              # up: 0..n-1, dn: n..2n-1
    cd = [m.T for m in c]
    H = np.zeros((2 ** nso, 2 ** nso))
    for i in range(n):
        for j in range(n):
            if abs(Ku[i, j]) > 1e-15: H += Ku[i, j] * (cd[i] @ c[j])
            if abs(Kd[i, j]) > 1e-15: H += Kd[i, j] * (cd[n + i] @ c[n + j])
    Nop = np.zeros_like(H); Sz = []
    for i in range(n):
        nu = cd[i] @ c[i]; nd = cd[n + i] @ c[n + i]
        H += U * (nu @ nd) - mu * (nu + nd)
        Nop += nu + nd
        Sz.append(0.5 * (nu - nd))
    w, V = np.linalg.eigh(H)
    bw = np.exp(-beta * (w - w.min())); Z = bw.sum()
    dens = float(np.sum(bw * np.einsum("ik,ij,jk->k", V, Nop, V)) / Z) / n
    # Lehmann kernel K(Ea,Eb) = (e^{-bEb}-e^{-bEa})/(Ea-Eb) -> beta e^{-bEa} at Ea=Eb
    diff = w[:, None] - w[None, :]
    small = np.abs(diff) < 1e-9
    Kk = np.where(small, beta * bw[:, None], (bw[None, :] - bw[:, None]) / np.where(small, 1.0, diff))
    xm = np.arange(n) // ly; ym = np.arange(n) % ly
    chi_q = np.zeros((lx, ly)); S_q = np.zeros((lx, ly))
    for kx in range(lx):
        for ky in range(ly):
            phase = np.exp(-2j * np.pi * (kx * xm / lx + ky * ym / ly))
            Sq = sum(phase[i] * Sz[i] for i in range(n))
            Mt = V.T @ Sq @ V                       # real V; Sq complex
            chi_q[kx, ky] = float(np.real((np.abs(Mt).T ** 2 * Kk).sum() / Z)) / n
            SS = Sq @ Sq.conj().T
            S_q[kx, ky] = float(np.real(np.sum(bw * np.einsum("ik,ij,jk->k", V, SS, V)) / Z)) / n
    return chi_q, S_q, dens


# ------------------------------------------------------------------------------ CLI --
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lx", type=int, default=2); ap.add_argument("--ly", type=int, default=2)
    ap.add_argument("--U", type=float, default=4.0); ap.add_argument("--mu", type=float, default=2.0)
    ap.add_argument("--beta", type=float, default=2.0); ap.add_argument("--dt", type=float, default=0.0625)
    ap.add_argument("--tam", type=float, default=0.0); ap.add_argument("--t1", type=float, default=0.0)
    ap.add_argument("--tp", type=float, default=0.0)
    ap.add_argument("--nwarm", type=int, default=300); ap.add_argument("--nmeas", type=int, default=1500)
    ap.add_argument("--seed", type=int, default=1); ap.add_argument("--ed", action="store_true")
    a = ap.parse_args()
    q = SpinDQMC(a.lx, a.ly, a.U, a.mu, a.beta, a.dt, a.tam, a.t1, a.seed, tp=a.tp)
    r = q.run_spin(a.nwarm, a.nmeas)
    print(f"# spin_susc {a.lx}x{a.ly} U={a.U} mu={a.mu} beta={a.beta} (NT={q.NT}) "
          f"tam={a.tam} t1={a.t1} tp={a.tp}")
    print(f"  <sign> = {r['sign']:.4f}   density = {r['dens']:.5f}")
    print(f"  DQMC chi_zz(q) grid (rows kx=0..):\n{np.array2string(r['chi_q'], precision=4)}")
    print(f"  DQMC S^z(q)   grid:\n{np.array2string(r['S_q'], precision=4)}")
    print(f"  chi_max={r['chi_max']:.4f}  chi(q=0)={r['chi_q0']:.4f}  S_max={r['S_max']:.4f}")
    if a.ed:
        chi_q, S_q, dens = ed_chi_spin(a.lx, a.ly, a.U, a.mu, a.beta, a.tam, a.t1, a.tp)
        print(f"  ED  density = {dens:.5f}")
        print(f"  ED  chi_zz(q) grid:\n{np.array2string(chi_q, precision=4)}")
        print(f"  ED  S^z(q)   grid:\n{np.array2string(S_q, precision=4)}")


if __name__ == "__main__":
    main()
