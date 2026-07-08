#!/usr/bin/env python3
"""Pairing-matrix eigenvalue (BSE-flavored, equal-time, q=0 COM) for the DQMC -- a DIRECT
pairing measure that does NOT project onto a fixed d-wave form factor and is a robust
extremal quantity (leading eigenvalue) rather than a fragile scalar projection.

Per config with Green's cu=<c^+_i c_j>_up, cd=<c^+_i c_j>_dn:
  pair matrix  P(k,k') = <b^+_k b_{k'}>,  b^+_k = c^+_{k up} c^+_{-k dn}
  Wick:  P(k,k') = (1/N^2) cu(k,k') * cd(-k,-k')   [direct term; anomalous=0]
  with cu(k,k') = sum_{il} e^{-i k r_i} cu[i,l] e^{+i k' r_l}  (= W cu W^dag, W=DFT).
CONNECTED Pc = <P> - <cu(k,k')>*<cd(-k,-k')>(bubble, from config-avg Green's).
Leading eigenvalue lambda_max of Pc (Hermitian) + eigenvector phi(k); classify phi by its
overlap with the d_{x2-y2} (cos kx - cos ky) and ext-s (cos kx + cos ky) form factors.
"""
import os, sys
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from dqmc import DQMC


def _dft(lx, ly):
    n = lx * ly
    xs = np.arange(n) // ly; ys = np.arange(n) % ly
    kx = 2 * np.pi * (np.arange(lx)) / lx; ky = 2 * np.pi * (np.arange(ly)) / ly
    KX, KY = np.meshgrid(kx, ky, indexing="ij")
    kxf = KX.ravel(); kyf = KY.ravel()             # n k-points
    W = np.exp(-1j * (np.outer(kxf, xs) + np.outer(kyf, ys))) / np.sqrt(n)  # W[k,i]
    return W, kxf, kyf


def pair_eig(lx, ly, U, mu, beta, dt=0.0625, tam=0.0, t1=0.0, tp=0.0,
             nwarm=150, nmeas=600, seed=1):
    q = DQMC(lx, ly, U, mu, beta, dt, tam=tam, t1=t1, seed=seed, tp=tp)
    n = q.n
    W, kxf, kyf = _dft(lx, ly)
    for _ in range(nwarm):
        q.sweep()
    Pfull = np.zeros((n, n), dtype=complex)        # <P(k,k')>
    Acu = np.zeros((n, n), dtype=complex); Acd = np.zeros((n, n), dtype=complex)
    sw = 0.0
    for _ in range(nmeas):
        s = q.sweep()
        G = [q.green(0), q.green(1)]
        cu = np.eye(n) - G[0].T; cd = np.eye(n) - G[1].T          # <c^+_i c_j>
        cuk = W @ cu @ W.conj().T                                # cu(k,k')
        cdk = W @ cd @ W.conj().T                                # cd(k,k')
        # cd(-k,-k'): negate momenta -> index map
        idxm = _neg_k_index(lx, ly)
        cdk_m = cdk[np.ix_(idxm, idxm)]
        Pfull += s * (cuk * cdk_m) / n
        Acu += s * cuk; Acd += s * cdk_m
        sw += s
    Pfull /= sw; Acu /= sw; Acd /= sw
    Pc = Pfull - (Acu * Acd) / n                                 # connected pairing matrix
    Pc = 0.5 * (Pc + Pc.conj().T)                                # Hermitize
    w, v = np.linalg.eigh(Pc)
    lam = w[-1]; phi = v[:, -1]                                  # leading eigval + eigvec
    # classify eigenvector by form-factor overlap
    fd = np.cos(kxf) - np.cos(kyf); fs = np.cos(kxf) + np.cos(kyf)
    od = abs(np.vdot(fd / np.linalg.norm(fd), phi))
    os_ = abs(np.vdot(fs / np.linalg.norm(fs), phi))
    return dict(lam=float(lam), d_overlap=float(od), s_overlap=float(os_),
                lam2=float(w[-2]), sign=sw / nmeas)


def _neg_k_index(lx, ly):
    n = lx * ly; out = np.empty(n, dtype=int)
    for a in range(lx):
        for b in range(ly):
            out[a * ly + b] = ((-a) % lx) * ly + ((-b) % ly)
    return out


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--lx", type=int, default=4); ap.add_argument("--ly", type=int, default=4)
    ap.add_argument("--U", type=float, default=4.0); ap.add_argument("--mu", type=float, default=2.0)
    ap.add_argument("--beta", type=float, default=4.0); ap.add_argument("--dt", type=float, default=0.0625)
    ap.add_argument("--tam", type=float, default=0.0); ap.add_argument("--t1", type=float, default=0.0)
    ap.add_argument("--nmeas", type=int, default=600); ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    r = pair_eig(a.lx, a.ly, a.U, a.mu, a.beta, a.dt, a.tam, a.t1, nmeas=a.nmeas, seed=a.seed)
    print(f"# pair-eig {a.lx}x{a.ly} U={a.U} mu={a.mu} beta={a.beta} tam={a.tam} t1={a.t1}")
    print(f"  lambda_max = {r['lam']:.4f}  (2nd = {r['lam2']:.4f})")
    print(f"  eigvec overlap:  d-wave = {r['d_overlap']:.3f}   ext-s = {r['s_overlap']:.3f}")
    print(f"  <sign> = {r['sign']:.4f}")


def pair_eig_tau(lx, ly, U, mu, beta, dt=0.0625, tam=0.0, t1=0.0, tp=0.0,
                 nwarm=120, nmeas=300, seed=1, stride=2):
    """TAU-INTEGRATED pairing-matrix eigenvalue (the SC-relevant instability indicator):
    P(k,k') = int_0^beta dtau Gu(k,k';tau) Gd(-k,-k';tau), connected (minus bubble of the
    config-averaged time-displaced Green's), leading eigenvalue + d/ext-s overlap. Uses the
    STABLE _green_tau. tau-stride to cut the number of stabilized solves."""
    q = DQMC(lx, ly, U, mu, beta, dt, tam=tam, t1=t1, seed=seed, tp=tp)
    n, NT = q.n, q.NT
    W, kxf, kyf = _dft(lx, ly); idxm = _neg_k_index(lx, ly)
    for _ in range(nwarm):
        q.sweep()
    ls = list(range(0, NT, stride)); wtau = stride * dt
    Pf = np.zeros((n, n), complex)
    AGu = {l: np.zeros((n, n), complex) for l in ls}
    AGd = {l: np.zeros((n, n), complex) for l in ls}
    sw = 0.0
    for _ in range(nmeas):
        s = q.sweep()
        for l in ls:
            Gu = W @ q._green_tau(0, l) @ W.conj().T
            Gd = W @ q._green_tau(1, l) @ W.conj().T
            Gdm = Gd[np.ix_(idxm, idxm)]
            Pf += s * wtau * (Gu * Gdm)
            AGu[l] += s * Gu; AGd[l] += s * Gdm
        sw += s
    Pf /= sw
    Pbub = sum(wtau * (AGu[l] / sw) * (AGd[l] / sw) for l in ls)
    Pc = Pf - Pbub; Pc = 0.5 * (Pc + Pc.conj().T)
    w, v = np.linalg.eigh(Pc); lam = w[-1]; phi = v[:, -1]
    fd = np.cos(kxf) - np.cos(kyf); fs = np.cos(kxf) + np.cos(kyf)
    od = abs(np.vdot(fd / np.linalg.norm(fd), phi)); os_ = abs(np.vdot(fs / np.linalg.norm(fs), phi))
    return dict(lam=float(lam), lam2=float(w[-2]), d_overlap=float(od), s_overlap=float(os_),
                sign=sw / nmeas)
